%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory


%% Try to synchronize using strobed values
% OpenEphys board timing
dataPath = '/mnt/MGS/Ephys/Raw/taco_20240712/Record Node 101/experiment1/recording1/events/Acquisition_Board-114.Rhythm Data/TTL/';
timestampsOE = readNPY(fullfile(dataPath,'timestamps.npy'));
statesOE      =readNPY(fullfile(dataPath,'states.npy'));
full_wordsOE  =readNPY(fullfile(dataPath,'full_words.npy'));
sample_numbersOE =readNPY(fullfile(dataPath,'sample_numbers.npy'));

%% NPX imec timing
dataPath = '/mnt/MGS/Ephys/Raw/taco_20240712/Record Node 101/experiment1/recording1/events/Neuropix-PXI-100.ProbeA-AP/TTL/';
timestampsNPX = readNPY(fullfile(dataPath,'timestamps.npy'));
statesNPX      =readNPY(fullfile(dataPath,'states.npy'));
full_wordsNPX  =readNPY(fullfile(dataPath,'full_words.npy'));
sample_numbersNPX =readNPY(fullfile(dataPath,'sample_numbers.npy'));

%% First sync imec to OE
riseSent=sample_numbersNPX(statesNPX==1);
fallSent=sample_numbersNPX(statesNPX==-1);

%% Analog signals
% %Continous data on OE board saved as 32channels of
% %headstage and 8 channels of ADC for 40 channels x two bytes per channel,
% %80 bytes per sample. 30000samples/sec,2.4e6bytes per sec
% OEdaqdata=memmapfile('/mnt/MGS/Ephys/Raw/taco_20240712/Record Node 101/experiment1/recording1/continuous/Acquisition_Board-114.Rhythm Data/continuous.dat');
% %% test plots of channels to find the 1Hz imec signal
% hold off
% for ii=33:40
%     start=2.4e6*3600; %1hr into recording
%     Analogtest=double(nidaqdata.Data(start+((2*ii-1):80:2.4e7)))*1+double(nidaqdata.Data(start+((2*ii):80:2.4e7)))*2^8;plot(Analogtest)
%     title(num2str(ii))
%     pause(1)
%     hold on
% end
% %not seeing it...
% hold off
%%
riseRec=sample_numbersOE(statesOE==8);
fallRec=sample_numbersOE(statesOE==-8);

plot(riseSent,riseRec,'.')
%% Remove signals from before the start of recording, (only needed if using time in seconds rather than samples)
riseSent=riseSent(riseRec~=-1);
riseRec=riseRec(riseRec~=-1);

fallSent=fallSent(fallRec~=-1);
fallRec=fallRec(fallRec~=-1);

plot(riseSent,riseRec,'.')
%maximum drift 
max(double(abs(riseSent-riseRec)))/30000

%% Convert OE timestamps into matched ephys time
sample_numbersOE_on_NPXclock = interp1(double(riseRec),double(riseSent),double(sample_numbersOE));

%% Putting it into equivalent .kwe format
eventId = double(sign(statesOE)==1);
time_samples = sample_numbersOE_on_NPXclock; %this might be expecting samples not time
event_channels = double(abs(statesOE)-1);

%% remove last channel (7) that was used for 1Hz syncing
eventId=eventId(event_channels~=7);
time_samples=time_samples(event_channels~=7);
event_channels=event_channels(event_channels~=7);

%% save out digital events
%load events only
timestamps = time_samples;%hdf5read(filename, '/event_types/TTL/events/time_samples');
highlow = eventId;%hdf5read(filename, '/event_types/TTL/events/user_data/eventID');
bitNumber = event_channels;%hdf5read(filename, '/event_types/TTL/events/user_data/event_channels');

%no longer have a strobe signal, so need to fake it
strobeSet=find(highlow==1); %any bit high
strobeUnset=find(highlow==0); %any bit goes low
strobeUnset=[1; strobeUnset];

value=nan(size(strobeSet));
for iStrobe=1:length(strobeSet)
     ts=timestamps <= timestamps(strobeSet(iStrobe)) & timestamps >= timestamps(strobeUnset(iStrobe)) & bitNumber~=7;
     value(iStrobe)=sum(2.^bitNumber(ts) .* highlow(ts));
end

eventTimes=double(timestamps(strobeSet));
eventValues = value;
flagBits = nan(size(value));
flagData = value;
invertedBits = false;
switchbits = false;

%eventSamples depends on recording number.
% location='/recordings';
% info=h5info([inDir filesep name],location);  %where file was filename=[inDir filesep name];
% nRecs=length(info.Groups);
% 
% st_index        = double(sample_numbersNPX(1));%strcmp('start_time',{info.Groups(iRec).Attributes.Name});
% recStartTime    = double(timestampsNPX(1));%double(info.Groups(iRec).Attributes(st_index).Value);
% eventSamples    = eventTimes-recStartTime;

%save([outDir filesep 'events_' exptDate '.mat'], 'eventTimes', 'eventSamples', 'eventValues', 'flagBits', 'flagData', 'invertedBits', 'switchbits','recStartTime');



%% Loadup PTB signals sent through datapixx
load('/home/huklab/Documents/MGS/pdsdata')
PTBts=cellfun(@(X) X.datapixx.unique_number_time(:,1), PDS.data,'uniformOutput',false);
PTBword=cellfun(@(X) X.unique_number,PDS.data,'UniformOutput',false);
for ii=1:length(PTBword)
PTBword{ii}(1)=PTBword{ii}(1);
end
ptbwords=cell2mat(PTBword);
ptbt=cell2mat(PTBts')';
%%
binarisedOEwords=dec2bin(full_wordsOE);
redone_full_wordsOE=bin2dec(binarisedOEwords(:,2:end));

%%
subplot(211)
plot(ptbt,ptbwords,'.')
title('Unique signals sent by PTB')
subplot(212)
plot(timestampsOE,redone_full_wordsOE,'.')
title('Signals recorded by OE')
%%
% size(event_channels)
% size(timestamps)

tmin=2179;
tmax=2501;

keep=(timestampsOE>tmin)&(timestampsOE<tmax);
RecWords=redone_full_wordsOE(keep);
RecWords_ts=timestampsOE(keep);

ptbt_=ptbt;
%%
%If you drop a bit or think a bit is flipped this is where you can try to
%debug it by passing through dec2bin 
sentbits=dec2bin(ptbwords);

%%

%I think all bits are sent to high at some point during the presentation
%but its not recorded as part of the words so we leave them out
% strangehighs=(RecWords==62|RecWords==63);
% RecWords=RecWords(~strangehighs);
% RecWords_ts=RecWords_ts(~strangehighs);

subplot(211)
scatter(RecWords_ts,RecWords,'.')

%xlim([6280 6680]); ylim([0 64]);
subplot(212)
plot(ptbt_,ptbwords,'.')

ptbtime0=min(ptbt_);
ptbtime1=max(ptbt_);
%Cropping can be usefull
% xlim([ptbtime0 ptbtime0]); ylim([0 64]);

%%
Recbits=dec2bin(RecWords);

%% Bits 4/2 and 3/5 are mixed up: HOW!?! new board?
switchbits=Recbits;
switchbits(:,7)=Recbits(:,7);
switchbits(:,6)=Recbits(:,6);
switchbits(:,5)=Recbits(:,5); %%!!!!
switchbits(:,4)=Recbits(:,4); %%!!!!
switchbits(:,3)=Recbits(:,3); %%%!!!
switchbits(:,2)=Recbits(:,2); %%!!!!
switchbits(:,1)=Recbits(:,1);

subplot(211)
scatter(RecWords_ts,bin2dec(switchbits(:,[1 2 3 4 5 6 7])),'.')
title('Signals recorded by OE')
xlim([2150 2550]); %
ylim([0 60]);
subplot(212)

plot(ptbt_,bin2dec(sentbits(:,[5 6 7  8 9 10 11])),'.')
title('Unique signals sent by PTB')
ylim([0 60]);
%xlim([ptbtime0 ptbtime1]); ylim([0 64]);
%%    STOPPED HERE, SIGNALS STILL AREN'T RIGHT
% %% if nec, update
% Recbits=switchbits;
% RecWords=bin2dec(Recbits);
% 
% %% 
% strangehighs=(RecWords==62|RecWords==63);
% RecWords=RecWords(~strangehighs);
% RecWords_ts=RecWords_ts(~strangehighs);
% 
% %% Rough check of range, crop time to after start of recording
% 
% ptbtime0=min(ptbt_)+6000;
% ptbtime1=max(ptbt_);
% 
% 
% subplot(211)
% scatter(RecWords_ts,RecWords,'.')
% ylim([0 64]);
% subplot(212)
% plot(ptbt_,ptbwords,'.')
% ylim([0 64]);
% 
% xlim([ptbtime0 ptbtime1]); % ylim([0 64]);
% 

% 
% %% Remove unrecorded times for syncing
% 
% ptbt=ptbt_(ptbt_>ptbtime0);
% ptbwords=ptbwords(ptbt_>ptbtime0);
% 
% %%
% figure(2);clf
% subplot(211)
% scatter(RecWords_ts,RecWords,'.')
% ylim([0 64]);
% subplot(212)
% plot(ptbt,ptbwords,'.')
% ylim([0 64]);
% 
% %% if you could find the same number of strobes you'd be right, 
% % but unlikely with different sampling rates
% %potential_matches = find(RecWords==45 & ptbwords==45); etc
% 
% %% From here could take the mode within time bins and xcor?
%     % first coarsely,10ms, then refined to 30kHz with less samples, make
%     % sure they start close enough to be in the sliding window
% offset=0;
% lowbound=min(RecWords_ts);
% highbound=(max(RecWords_ts)+500);
% middletime=(highbound-lowbound)/2;
% for res=[0.01 0.001 0.0001]
% %
%     alltime_0=ptbt-ptbt(1)+offset;
%     
%     edges = middletime-(500000*res):res:middletime+(500000*res);
%     md_sg= zeros(1,length(edges)-1);
%     for bin = 1:length(edges)-1
%     ed1=edges(bin);
%     ed2=edges(bin+1);
%     idx=RecWords_ts>ed1&RecWords_ts<ed2;
%     md_sg(bin)=mode(RecWords(idx));
%     end
%     %
%     md_sg_cl= zeros(1,length(edges)-1);
%     
%     for bin = 1:length(edges)-1
%     ed1=edges(bin);
%     ed2=edges(bin+1);
%     idx=alltime_0>ed1&alltime_0<ed2;
%     md_sg_cl(bin)=mode(ptbwords(idx));
%     end
%     
%     %
%     figure(2);clf
%     middles=edges(1:end-1)+res/2;
%     subplot(211)
%     plot(middles,md_sg,'.')
%     subplot(212)
%     title('sent')
%     plot(middles,md_sg_cl,'.')
%     %
%     md_sg(isnan(md_sg))=0;
%     md_sg_cl(isnan(md_sg_cl))=0;
%     
%     % Get best lag to match time series
%     [out, lags] = xcorr(md_sg, md_sg_cl, 'None');
%     [~, id] = max(out);
%     %
%     bin_size=mode(diff(edges));
%     figure(4); clf; set(gcf, 'Color', 'w')
%     plot(lags*bin_size, out); hold on
%     xlabel('Lag (seconds)')
%     ylabel('Cross Correlation')
%     plot(lags(id)*bin_size, out(id), 'o')
%     legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')
%     add_offset = lags(id)*bin_size;
%     offset=offset+add_offset;
%     fprintf('Best guess for initial offset: %02.2f (s)\n', offset)
%     
%     
%     
%     %
%     subplot(111)
%     plot(alltime_0+offset,ptbwords,'.'); hold on
%     plot(RecWords_ts,RecWords,'.'); hold off
%     xlim([0 12000])
% 
%     %end
%    
% 
% end
% 
% %% Total offset
% toffset=offset-ptbt(1);
% sprintf('%15.15f',toffset)
% 
% w=[1 -toffset];
% @(t)(t-w(2))/w(1)
% myfun=@(t)(t-w(2))/w(1);
% 
% %%
% plot(ptbst+myfun(0),(ptbsw),'.')
% hold on
% pause(0.1)
% scatter(RecWords_ts,RecWords,'.')
% hold off
% 
%% check
% clf
offset=2179.12-1720821827.75;
ptbtimeend=max(RecWords_ts-offset);
lowbound=2150;
highbound=2550; 
x1=lowbound; x2=highbound;
hold off
for ll=1:1
    figure(3);subplot(111)
    plot(ptbt,(ptbwords),'.')
    hold on
    xlim([x1-offset x2-offset]); ylim([0 64]);
    pause(0.2)
    scatter(RecWords_ts-offset,RecWords,'o')
    plot([ptbtimeend ptbtimeend],[0 64])
    hold off
    xlim([x1-offset x2-offset]); ylim([0 64]);
    pause(0.2)
end
% 
% 
% 
%% Remove unrecorded times for syncing
ptbt=ptbt(ptbt<ptbtimeend);
ptbwords=ptbwords(ptbt_<ptbtimeend);
%% Correcting for drift, best to take matches

wordstocheck=15;%unique(ptbwords);
clear tdiff rec_match ptb_match
for ii=1:length(wordstocheck)
    word=wordstocheck(ii);
    ptb_ind=find(ptbwords==word);
    rec_ind=find(RecWords==word);
    if nnz(ptb_ind)>nnz(rec_ind) % more ptb, loop through rec
        for jj=1: nnz(rec_ind)
            [tdiff{ii}(jj),ptb_match{ii}(jj)]=min(abs(ptbt+offset-RecWords_ts(rec_ind(jj))));
            rec_match{ii}(jj)=rec_ind(jj);
        end
    else % more rec, loop through ptb
        for jj=1: nnz(ptb_ind)
            [tdiff{ii}(jj),rec_match{ii}(jj)]=min(abs(ptbt(ptb_ind(jj))+offset-RecWords_ts));
            ptb_match{ii}(jj)=ptb_ind(jj);
        end
    end
    max(tdiff{ii})

end
pmatches=cell2mat(ptb_match);
ematches=cell2mat(rec_match);

fun = synchtime.align_clocks(ptbt(pmatches)', RecWords_ts(ematches));
plot(ptbt(pmatches), RecWords_ts(ematches), 'o'); hold on
plot(xlim, fun(xlim), 'k')
xlabel('PTB clock')
ylabel('Ephys clock')
%%
clf
x1=lowbound; x2=highbound;
for ll=1:10
figure(2);subplot(111)
plot(fun(ptbt),(ptbwords),'.')
hold on
xlim([x1 x2]); ylim([0 64]);
pause(0.2)
scatter(RecWords_ts,RecWords,'o')
hold off
xlim([x1 x2]); ylim([0 64]);
pause(0.2)

end

%% If happy lock it in!!!
Exp.ptb2Ephys = fun;

%%
