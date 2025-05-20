%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory
cd('/home/huklab/Documents/NPX_pilot')

%% 
ephysPath='/mnt/NPX/Rocky/20230921/Rocky20230921_V2_g0';
ephysPrefix='Rocky20230921_V2_g0';

 %% Get timing strobes from spikeGLX setup
% Nidq.bin has data stream sampled on nidaq. We can take the start and end
% times from the meta file or we can treat it as its own stream starting at
% t=0s since we need to register and shift other clocks to this anyway. We
% have 4 clocks that need to be alligned (excluding eyetracking): 

%   imec card on which the ephys is recorded, let's make this the master clock that everything is synced to
%   Nidaq (gets datapix,PD,imec timing triggers), for timing
%   PTB - stimulus machine clock
%   datapixx - connected to PTB and not expected to introduce large lags,
%   this final clock can be replaced with some other output hardware 



 %timingfilename='/mnt/NPX/Rocky/20230921/Rocky20230921_V2_g0/Rocky20230921_V2_g0_t0.nidq.bin';
timingfilename=[ephysPath '/' ephysPrefix '_t0.nidq.bin'];
nchannels=4;
samples=[];
dataformat=[];
memorymap=1;
NiDaQout = spGLX_load_binary_file(timingfilename,nchannels,samples,dataformat,memorymap);

%Recorded digital lines are stored in: 
ttl=dec2bin(NiDaQout.Data.data(4,:));
datapixx=ttl(:,1:7)=='1'; % datapixx strobing uses 7 bits, 

%Trying to recreate 'event' recordings ala plexon, a strobe bit gates
%recording of events. Here we just want to collect any highs on these wires
Event.samplenumber=find(sum(datapixx,2)>0);
Event.Ni_ts=(Event.samplenumber)/30002.845833;% from nidaq metafile
Event.signal=ttl(Event.samplenumber,1:7);

%Final signal is recorded on first bit, least significant digit so last in
%byte order
imecSync=ttl(:,8)=='1';
riseRec=find(diff(imecSync)>0); %Recorded rising edges in samples on nidaq

%This is typically a 1Hz signal from the imec card that signals the start
%of recording on the first rising edge and maintains a pulse to correct for
%drifting of clocks between the imec and nidaq cards

%Analog signals are also recorded on the nidaq:
PDanalog=NiDaQout.Data.data(1,:);
PD=PDanalog>1.2*10^4;
risePD=find(diff(PD)>0);% Rising edges
risePD_Ni_ts=(risePD)/30002.845833;
fallPD=find(diff(PD)<0);%Falling edges
fallPD_Ni_ts=(fallPD)/30002.845833;% 



% All signals are now in nidaq samples, and need to be converted to imec
% timing, using the 1Hz sync signal as reference


%% Received signals from imec card for syncing 
%To load the sent imec signals, need to load up the entire raw spiking
%file, this may actually be unnecessary if we are very confident in a lack
%of drift between the imec-nidaq clocks and take start times from metafiles

spikefilename=[ephysPath '/' ephysPrefix '_imec0/' ephysPrefix '_t0.imec0.ap.bin'];
nchannels=385;
samples=[];
dataformat=[];
memorymap=1;
Tout = spGLX_load_binary_file(spikefilename,nchannels,samples,dataformat,memorymap);

%Signals that were sent from imec card recorded on the ephys clock
imecSent=Tout.Data.data(385,:); %this takes forever
riseSent=find(diff(imecSent)>40); %in samples on AP



%%
%plot(riseRec,riseSent')
%save('/home/huklab/Documents/SpikeSorting/Output/2023-09-21/Rocky20230921_V2_g0_t0.nidq_imecAp.mat','riseSent','riseRec','imecSent','imecSync','risePD','fallPD')
% clear Tout
 load('/home/huklab/Documents/SpikeSorting/Output/2023-09-21/Rocky20230921_V2_g0_t0.nidq_imecAp.mat')

%% Converting signals from nidaq to ephys(imec) time
% Datapixx 'events'
Event.ClosestAPSamples = interp1(riseRec,riseSent,Event.samplenumber);
Event.Ap_ts=Event.ClosestAPSamples/30000.074166666665; % from AP metafile

%Photodiode rise times
%risephD=find(diff(phDiode)>0); %in samples on nidaq
Event.risePhDiodeClosestAPSamples = interp1(riseRec,riseSent,risePD);
Event.risePD_Ap_ts=Event.risePhDiodeClosestAPSamples/30000.074166666665; % from AP metafile
Event.fallPhDiodeClosestAPSamples = interp1(riseRec,riseSent,fallPD);
Event.fallPD_Ap_ts=Event.fallPhDiodeClosestAPSamples/30000.074166666665; % from AP metafile

% Saving everything into an Events structure for saving / ease of loading


%% Events recorded on NiDaq are now in terms of seconds on Ephys(Imec) clock from t=0s at start of NPX recording

%save(['/home/huklab/Documents/SpikeSorting/Output/2023-09-21/' 'Rocky20230921_V2_t0.nidq' 'dp_ev_ts2.mat'],'Event')
load(['/home/huklab/Documents/SpikeSorting/Output/2023-09-21/' 'Rocky20230921_V2_t0.nidq' 'dp_ev_ts2.mat'],'Event')


%% Syncing between PTB and imec (using signals recorded on nidaq)
% Strobe signal goes high on sent strobes, this gates the recording of 
Strobebit=str2num(Event.signal(:,7));
RecStrobes=[diff(Strobebit); 0]>0; %just take rising edges to prevent duplicates

%%
RecStrobes=(Event.signal(:,7)=='1');
strobeIdx=find(RecStrobes);
Recbits=(Event.signal(strobeIdx,1:6));
RecWords=bin2dec(Event.signal(strobeIdx,1:6));
RecWords_ts=Event.Ap_ts(strobeIdx);

%%

nStrobes=sum(RecStrobes)
t=Event.Ap_ts(RecStrobes);
%%
% Recbits=Event.signal(:,1:6);
% RecWords=Recbits(sum(Recbits,2)>0,:);
% RecWords=bin2dec(num2str(RecWords));
% RecWords_ts=Event.Ap_ts(sum(Recbits,2)>0);
%%
AllWords=bin2dec(Event.signal(:,1:6));
scatter(Event.Ap_ts,AllWords,'.')
hold on
scatter(RecWords_ts,RecWords,'.')
hold off
title('All rec. words vs Strobed words')

%% We should now have a time and a word that we can match to the sent messages from ptb

%% If any stim files crashed out or didn't save to *z.mat fix them here
% zdata=load('/mnt/NPX/Rocky/20230921/Rocky 2023-09-21/FaceCal_Rocky_210923_04.mat');
% S=zdata.S;
% D=cell(1,1);
% ND=length(fields(zdata));
% for k=1:(ND-1)
% Dstring=sprintf('D%d',k);
% D{k,1}= zdata.(Dstring);
% end
% 
% %save locally then move to mnt
% save('FaceCal_Rocky_210923_04z.mat','S','D')


%% Load stimuli data from once data was recording, move calibration/setup etc to subfolder
StimPath='/mnt/NPX/Rocky/20230921/Rocky 2023-09-21/';
Exp = io.basic_marmoview_import(StimPath);

%% Remove weird trial(s)
jj=0;
for ii=1:length(Exp.D)
    jj=jj+1;
    if(isfield(Exp.D{ii},'STARTCLOCK'))
        Exp2.D{jj,1}=Exp.D{ii};
    else
        jj=jj-1; %skip
    end
end

%% just one
Exp.D=Exp2.D;
clear Exp2

%% Check number of timestamps
ptbsw = cell2mat(cellfun(@(x) x.STARTCLOCK(:)', Exp.D, 'uni', 0));
ptbew = cell2mat(cellfun(@(x) x.ENDCLOCK(:)', Exp.D, 'uni', 0));

ptbwords = reshape([ptbsw'; ptbew'], [], 1); % all strobe words

fprintf('PTB sent %d strobed words during this session\n', numel(ptbwords))
fprintf('Ephys recovered %d strobes (flip on channel 3)\n', nStrobes)

%% Prepare to synch using xcorr
ptbst = cellfun(@(x) x.STARTCLOCKTIME, Exp.D); % ptb trial starts
ptbet = cellfun(@(x) x.ENDCLOCKTIME, Exp.D); % ptb trial ends
ptbt = reshape([ptbst'; ptbet'], [], 1); % all strobe times from ptb

% so we have to spoof timings or set them to zero, they should be so fast
% that the 30kHz recording may not separate them

ptbt_=ptbt;
ptbt_=(ptbt_'+[0:0.001:0.005]');
ptbt_=ptbt_(:);

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
title('Received Signals')
%xlim([6280 6680]); ylim([0 64]);
subplot(212)
plot(ptbt_,ptbwords,'.')
%xlim([1694139100 1694139500]); ylim([0 64]);
title('Sent PTB Signals')

%% Bits 4/2 and 3/5 are mixed up: HOW!?! new board?
switchbits=Recbits;
switchbits(:,6)=Recbits(:,6);
switchbits(:,5)=Recbits(:,5); %%!!!!
switchbits(:,4)=Recbits(:,4); %%!!!!
switchbits(:,3)=Recbits(:,3); %%%!!!
switchbits(:,2)=Recbits(:,2); %%!!!!
switchbits(:,1)=Recbits(:,1);

subplot(211)
scatter(RecWords_ts,bin2dec(switchbits(:,1:6)),'.')
%xlim([6280 6680]); %ylim([0 64]);
subplot(212)




plot(ptbt_,bin2dec(sentbits(:,1:6)),'.')
%xlim([1694139100 1694139500]); % ylim([0 64]);

%% if nec, update
Recbits=switchbits;
RecWords=bin2dec(Recbits);

%% 
strangehighs=(RecWords==62|RecWords==63);
RecWords=RecWords(~strangehighs);
RecWords_ts=RecWords_ts(~strangehighs);

subplot(211)
scatter(RecWords_ts,RecWords,'.')
ylim([0 64]);
subplot(212)
plot(ptbt_,ptbwords,'.')
ylim([0 64]);

% xlim([1694132500 1694145000]); % ylim([0 64]);



%% Remove unrecorded times for syncing

ptbt=ptbt_(ptbt_>0);
ptbwords=ptbwords(ptbt_>0);

%%
figure(2);clf
subplot(211)
scatter(RecWords_ts,RecWords,'.')
ylim([0 64]);
subplot(212)
plot(ptbt,ptbwords,'.')
ylim([0 64]);

%% if you could find the same number of strobes you'd be right, 
% but unlikely with different sampling rates
%potential_matches = find(RecWords==45 & ptbwords==45); etc

%% From here could take the mode within time bins and xcor?
    % first coarsely,10ms, then refined to 30kHz with less samples, make
    % sure they start close enough to be in the sliding window
offset=0;
lowbound=min(RecWords_ts);
highbound=(max(RecWords_ts)+500);
middletime=(highbound-lowbound)/2;
for res=[0.01 0.001 0.0001]
%
    alltime_0=ptbt-ptbt(1)+offset;
    
    edges = middletime-(500000*res):res:middletime+(500000*res);
    md_sg= zeros(1,length(edges)-1);
    for bin = 1:length(edges)-1
    ed1=edges(bin);
    ed2=edges(bin+1);
    idx=RecWords_ts>ed1&RecWords_ts<ed2;
    md_sg(bin)=mode(RecWords(idx));
    end
    %
    md_sg_cl= zeros(1,length(edges)-1);
    
    for bin = 1:length(edges)-1
    ed1=edges(bin);
    ed2=edges(bin+1);
    idx=alltime_0>ed1&alltime_0<ed2;
    md_sg_cl(bin)=mode(ptbwords(idx));
    end
    
    %
    figure(2);clf
    middles=edges(1:end-1)+res/2;
    subplot(211)
    plot(middles,md_sg,'.')
    subplot(212)
    title('sent')
    plot(middles,md_sg_cl,'.')
    %
    md_sg(isnan(md_sg))=0;
    md_sg_cl(isnan(md_sg_cl))=0;
    
    % Get best lag to match time series
    [out, lags] = xcorr(md_sg, md_sg_cl, 'None');
    [~, id] = max(out);
    %
    bin_size=mode(diff(edges));
    figure(4); clf; set(gcf, 'Color', 'w')
    plot(lags*bin_size, out); hold on
    xlabel('Lag (seconds)')
    ylabel('Cross Correlation')
    plot(lags(id)*bin_size, out(id), 'o')
    legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')
    add_offset = lags(id)*bin_size;
    offset=offset+add_offset;
    fprintf('Best guess for initial offset: %02.2f (s)\n', offset)
    
    
    
    %
    subplot(111)
    plot(alltime_0+offset,ptbwords,'.'); hold on
    plot(RecWords_ts,RecWords,'.'); hold off
    xlim([0 12000])

    %end
   

end

%% Total offset
toffset=offset-ptbt(1);
sprintf('%15.15f',toffset)

w=[1 -toffset];
@(t)(t-w(2))/w(1)
myfun=@(t)(t-w(2))/w(1);

%%
plot(ptbt+myfun(0),(ptbwords),'.')
hold on
pause(0.1)
scatter(RecWords_ts,RecWords,'.')
hold off

%% check
clf
x1=lowbound; x2=highbound;
for ll=1:10
figure(2);subplot(111)
plot(ptbt+myfun(0),(ptbwords),'.')
hold off
xlim([x1 x2]); ylim([0 64]);
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2]); ylim([0 64]);
pause(0.2)

end

%% If happy lock it in!!!
Exp.ptb2Ephys = myfun;





%% Check screen lag, Photodiode times vs 'Flashtimes' recorded
%% Photodiode timing checks, 
foragetrials=find(cellfun(@(x) strcmp(x.PR.name,'ForageProceduralNoise'),Exp.D));
trial_start_times=cellfun(@(x) x.STARTCLOCKTIME,Exp.D);
trial_end_times=cellfun(@(x) x.ENDCLOCKTIME,Exp.D);

Forage_time_windows=[trial_start_times(foragetrials) trial_end_times(foragetrials)];
cellfun(@(x) x.PR.Flashtime,Exp.D(foragetrials),'UniformOutput',false)
%Does Flashtime have the screen flip times? Or is it still the flip before
PD_PTBtimes=cellfun(@(x) x.PR.Flashtime,Exp.D(foragetrials),'UniformOutput',false);
PD_PTBtimes_=cell2mat(PD_PTBtimes);
PD_PTBinAPtimes=Exp.ptb2Ephys(PD_PTBtimes_);
%%
PD_PTBtimes_sent=[];
for trial=1:length(Exp.D)
    if isfield(Exp.D{trial}.PR,'Flashtime')
        PD_PTB{trial}=Exp.D{trial}.PR.Flashtime;
        
        PD_PTBtimes_sent=[PD_PTBtimes_sent; PD_PTB{trial}];
        
    end
end

PD_PTBinAPtimes=Exp.ptb2Ephys(PD_PTBtimes_sent);
%%
hold off
plot(Event.fallPD_Ap_ts,.95*ones(1,length(Event.fallPD_Ap_ts)),'.')
hold on
plot(Event.risePD_Ap_ts,1.05*ones(1,length(Event.risePD_Ap_ts)),'.')
plot(PD_PTBinAPtimes,ones(1,length(PD_PTBinAPtimes)),'.')
xlim([5150 5200]);ylim([0.8 1.2])

%% Analog PD signal on ephys time
nidaqtime_s=1:length(PDanalog)/30002.845833;
nidaqClosestAPSamples = interp1(riseRec,riseSent,1:length(PDanalog));
nidaqsampon_Ap_ts=nidaqClosestAPSamples/30000.074166666665; % from AP metafile
%%
PD_PTBinAPtimes(PD_PTBinAPtimes<0)=[];
%NOTE; Flashtimes were saved with the 'currenttime' before the draw of the screen with
%the photodiode, need to add in an additional screenflip (1/60)

%%
figure(11);clf
hold off
plot(nidaqsampon_Ap_ts(8e7:4:end-2e8),PDanalog(8e7:4:end-2e8))
hold on 
plot(PD_PTBinAPtimes+(1/60),4000*ones(1,length(PD_PTBinAPtimes)),'.')
%xlim([3000 9000]);
xlim([5350 5360]);
% 20-30ms? 25ms of system delay to add in to correct screen throws
% add screendelay for each stim, or just subtract 0.025 from spike times


%% 
% Exp=load('Rock20230921_V2_nodpi.mat')
%% Reorder D by trial time for ease
trial_start_times=cellfun(@(x) x.STARTCLOCKTIME,Exp.D);
[sortedtime reorder]=sort(trial_start_times);

Exp.D={Exp.D{reorder}}';



%% Interjecting to look at waves trials
foragetrials=cellfun( @(x) strcmp(x.PR.name,'ForageProceduralNoise') ,Exp.D);
forindx=find(foragetrials);
sumwaves=cellfun( @(x) x.PR.noisetype==8 ,Exp.D(foragetrials));
wvnoisetrials=forindx(find(sumwaves));
Dnoise={Exp.D{wvnoisetrials}};

%%
FullNoiseHistory=[];
for trials=1:length(Dnoise)
    FullNoiseHistory=cat(1,FullNoiseHistory,Dnoise{trials}.PR.NoiseHistory);
end



wvtime=FullNoiseHistory(:,1); 
%1 time entry then 4*6 wave params, ori, sf, phase, dir-90, contrast,speed, contrast

%probably should have saved the relative phase but this is honestly a nice
%way to check absolute vs relative.
% test:
wv1.ph=FullNoiseHistory(:,4);
wv1.ori=FullNoiseHistory(:,5);
wv1.sf=FullNoiseHistory(:,3);

phmod=wv1.sf/min(unique(FullNoiseHistory(:,3)));
%phase should be 90+phmod*Rphase
Rphase=(wv1.ph-90)./phmod;
unique(Rphase)
%Looks like some in/decrements from frame to frame shift

%%
%% Load spikes
dataPath = '/home/huklab/Documents/SpikeSorting/Output/2023-09-21/';
Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp.st = Exp.osp.st - .025;
Exp.osp.st(Exp.osp.st<0)=0;
% move zero to the end
id = max(Exp.osp.clu)+1;

Exp.osp.clu(Exp.osp.clu==0) = id;
Exp.osp.cids = Exp.osp.cids2;
Exp.osp.cids(Exp.osp.cids==0) = id;

% lbledgood=find(Exp.osp.cgs==2)+1;

% Exp.osp.st = st_s.spikeTimes(:,1)+.6;
% Exp.osp.clu = st_s.spikeTimes(:,2);
% Exp.osp.depths = st_s.depths;
% Exp.osp.templates = st_s.templates;
% Exp.osp.cids = unique(Exp.osp.clu);


figure(1); clf
plot.raster(Exp.osp.st, Exp.osp.clu, 1)
ylabel('Unit ID')
xlabel('Time (seconds)')


%%
wvtime_ephys_t=Exp.ptb2Ephys(wvtime);
sp=Exp.osp;

%Every frame
t0=wvtime_ephys_t; %Should ADD LATENCY HERE, but already subtracted from ?
t1=t0+1/60;
windowsize_=t1-t0;

[~, ~, ind1] = histcounts(sp.st, t0);
[~, ~, ind2] = histcounts(sp.st, t1);

ix = ind1 ~= 0 & ind2 ~= 0;
% ix = ix & (ind2-ind1 == mode(ind2-ind1));
ix = ix & (ind2+1 == ind1);
% ix = ix & sp.clu > 0;
ix = ix & sp.clu  > 0;
Y = sparse(ind1(ix), double(sp.clu(ix)), ones(sum(ix), 1), numel(t0), max(sp.clu));

Y = full(Y); %sp count in fixes x units
Y = Y./windowsize_;

%% Assign spikes to units

cids = Exp.osp.cids;

if isempty(intersect(unique(Exp.osp.clu),  Exp.osp.cids))
    warning('provided cids does not match unit ids')
    cids = unique(Exp.osp.clu);
end

Robs = Y(:,cids);

%% Singles
% Separate according to pairs? Six per stim frame -> 6x the Noisehistory
% size
Wavepairs_t=(reshape(repmat(wvtime,1,6)',[],1));
% w1-w2
% w1-w3
% w1-w4
% w2-w3
% w2-w4
% w3-w4

wv1=FullNoiseHistory(:,2:7);
wv2=FullNoiseHistory(:,7+(1:6));
wv3=FullNoiseHistory(:,1+2*6+(1:6));
wv4=FullNoiseHistory(:,1+3*6+(1:6));

Allwaves_t=(reshape(repmat(wvtime,1,4),[],1));% just repeating time series (not interleaved)
Allwaves_R=(reshape(repmat(wvtime,1,4),[],1));% just repeating time series (not interleaved)
Allwaves=([wv1;wv2;wv3;wv4]);

Allwaves(Allwaves(:,4)>=180,4)=Allwaves((Allwaves(:,4)>=180),4)-180;


Robs_all=(repmat(Robs,4,1));
Robs_all(1:10,1)
Robs_all((1:10)+size(Robs,1),1)

[oriset, ~, oriind]=unique(Allwaves(:,4));

% Basically no preference, given how each was merged with at least tow other oris it makes sense
for ii=1:length(oriset)
    Rori(ii,:)=mean(Robs_all(Allwaves(:,4)==oriset(ii),:));
end
mean(Rori,2)
%%
% Separate according to pairs? Six per stim frame -> 6x the Noisehistory
% size
Wavepairs_t=(reshape(repmat(wvtime,1,6)',[],1));
% w1-w2
% w1-w3
% w1-w4
% w2-w3
% w2-w4
% w3-w4

wv1=FullNoiseHistory(:,2:7);
wv2=FullNoiseHistory(:,7+(1:6));
wv3=FullNoiseHistory(:,1+2*6+(1:6));
wv4=FullNoiseHistory(:,1+3*6+(1:6));

Wavepairs_t=(reshape(repmat(wvtime,1,6),[],1));% just repeating time series (not interleaved)
Wavepairs=([[wv1 wv2];[wv1 wv3];[wv1 wv4];[wv2 wv3];[wv2 wv4];[wv3 wv4]]);

Wavepairs(Wavepairs(:,4)>=180,4)=Wavepairs((Wavepairs(:,4)>=180),4)-180;
Wavepairs(Wavepairs(:,10)>=180,10)=Wavepairs((Wavepairs(:,10)>=180),10)-180;

%Find simplest case: matching ori, but don't know best ori..
orimatched=find(Wavepairs(:,4)==Wavepairs(:,10));

subset=(Wavepairs(orimatched,4)==120);

Robs_all=(repmat(Robs,6,1));
Robs_all(1:10,1)
Robs_all((1:10)+size(Robs,1),1)

%%
MatchedPairs=Wavepairs(orimatched(subset),:);
Robs_MatchedPairs=Robs_all(orimatched(subset),:);

phasediff=MatchedPairs(:,3)-MatchedPairs(:,6+3);

phmod1=MatchedPairs(:,2)/min(unique(FullNoiseHistory(:,3)));
phmod2=MatchedPairs(:,2+6)/min(unique(FullNoiseHistory(:,3)));
Rphase1=(MatchedPairs(:,3)-90)./phmod1;
Rphase2=(MatchedPairs(:,9)-90)./phmod2;

for jj=1:10
MatchedPairs(MatchedPairs(:,3)>=360,3)=MatchedPairs((MatchedPairs(:,3)>=360),3)-360;
MatchedPairs(MatchedPairs(:,9)>=360,9)=MatchedPairs((MatchedPairs(:,9)>=360),9)-360;

MatchedPairs(MatchedPairs(:,3)<=-360,3)=MatchedPairs((MatchedPairs(:,3)<=-360),3)+360;
MatchedPairs(MatchedPairs(:,9)<=-360,9)=MatchedPairs((MatchedPairs(:,9)<=-360),9)+360;
end

Rphasediff=Rphase1-Rphase2;
% If below zero, add 360 to first wave
Rphasediff(Rphasediff<=-10)=Rphasediff(Rphasediff<=-10)+360;
phasediff=MatchedPairs(:,3)-MatchedPairs(:,6+3);
uniqueshifts=unique(phasediff);
binaround=0:45:315;%
%%
for ii=1:8
    ind_bin=abs(Rphasediff-binaround(ii))<15;
    phasebin{ii}=Robs(ind_bin,:);
    meanRphasepairs(ii,:)=mean(Robs(ind_bin,:));
    stdRphasepairs(ii,:)=std(Robs(ind_bin,:));

end
%% Not great, but the SNR is max supposed to be 1/6 of what's possible
goodunits=find(mean(meanRphasepairs,1)>1 & Exp.osp.clusterDepths>3000);

for cell=1:length(goodunits)
hold off
plot(binaround,meanRphasepairs(:,goodunits(cell)))
% hold on
% plot(binaround,meanRphasepairs(:,NC)+2.*stdRphasepairs(:,NC),'r-')
% plot(binaround,meanRphasepairs(:,NC)-2.*stdRphasepairs(:,NC),'r-')
ylim([0 max(meanRphasepairs(:,goodunits(cell)))])
pause(0.5)
end





%% From the top
%ori, sf, phase, dir-90 ,speed, contrast
%just sf, phase, ori, contrast
wv1=FullNoiseHistory(:,1+([2 3 4 6]));
wv2=FullNoiseHistory(:,1+6+([2 3 4 6]));
wv3=FullNoiseHistory(:,1+6*2+([2 3 4 6]));
wv4=FullNoiseHistory(:,1+6*3+([2 3 4 6]));

phmd1=(wv1(:,1)/0.5);
phmd2=(wv2(:,1)/0.5);
phmd3=(wv3(:,1)/0.5);
phmd4=(wv4(:,1)/0.5);

Rphase1=(wv1(:,2)-90)./phmd1;
Rphase2=(wv2(:,2)-90)./phmd2;
Rphase3=(wv3(:,2)-90)./phmd3;
Rphase4=(wv4(:,2)-90)./phmd4;

%%

% Take orientations > 180 and rotate? Need to invert phases from opposing
% directions
Rphase1(wv1(:,3)>=180)=360-Rphase1(wv1(:,3)>=180);
wv1(wv1(:,3)>=180,3)=wv1(wv1(:,3)>=180,3)-180;
Rphase2(wv2(:,3)>=180)=360-Rphase2(wv2(:,3)>=180);
wv2(wv2(:,3)>=180,3)=wv2(wv2(:,3)>=180,3)-180;
Rphase3(wv3(:,3)>=180)=360-Rphase3(wv3(:,3)>=180);
wv3(wv3(:,3)>=180,3)=wv3(wv3(:,3)>=180,3)-180;
Rphase4(wv4(:,3)>=180)=360-Rphase4(wv4(:,3)>=180);
wv4(wv4(:,3)>=180,3)=wv4(wv4(:,3)>=180,3)-180;

Rphase1(Rphase1>=345)=Rphase1(Rphase1>=345)-360;
Rphase2(Rphase2>=345)=Rphase2(Rphase2>=345)-360;
Rphase3(Rphase3>=345)=Rphase3(Rphase3>=345)-360;
Rphase4(Rphase4>=345)=Rphase4(Rphase4>=345)-360;

%Bin phases back to 8 starting phases
binaround=0:45:315;%
for ii=1:8
    ind_bin1=abs(Rphase1-binaround(ii))<15;
    Rphase1(ind_bin1)=binaround(ii);
    Rphind1(ind_bin1)=ii;
    ind_bin2=abs(Rphase2-binaround(ii))<15;
    Rphase2(ind_bin2)=binaround(ii);
    Rphind2(ind_bin2)=ii;
    ind_bin3=abs(Rphase3-binaround(ii))<15;
    Rphase3(ind_bin3)=binaround(ii);
    Rphind3(ind_bin3)=ii;
    ind_bin4=abs(Rphase4-binaround(ii))<15;
    Rphase4(ind_bin4)=binaround(ii);
    Rphind4(ind_bin4)=ii;
end

%Want to index into a space matrix with params for 
%mapping  %[0 60 120] to [1 2 3]
oriind1=wv1(:,3)/60+1;
oriind2=wv2(:,3)/60+1;
oriind3=wv3(:,3)/60+1;
oriind4=wv4(:,3)/60+1;

% mapping to 1:8
sfind1=phmd1;
sfind2=phmd2;
sfind3=phmd3;
sfind4=phmd4;
%% Big bad covariance matrix

nstim=size(FullNoiseHistory,1);
%Build a stimulus space matrix, ori*sf*rphase=192 dims, for each frame
%(oriind1-1)*8*8+(sfind1-1)*8+(Rphind1)
%ori indexes 3 blocks of 64, sfind into subblocks of 8 and so on
Stimmat=zeros(192,nstim);
Stimcovmat=zeros(192,192,nstim);
% Loop through each stimulus frame
% For each wave put contrast in the correct ori,sf and phase bin
for stim=1:nstim
    %Adding contrasts, in case of doubling up in bins
    Stimmat((oriind1(stim)-1)*8*8+(sfind1(stim)-1)*8+(Rphind1(stim)),stim)=...
        Stimmat((oriind1(stim)-1)*8*8+(sfind1(stim)-1)*8+(Rphind1(stim)),stim) + wv1(stim,4);
    Stimmat((oriind2(stim)-1)*8*8+(sfind2(stim)-1)*8+(Rphind2(stim)),stim)=...
        Stimmat((oriind2(stim)-1)*8*8+(sfind2(stim)-1)*8+(Rphind2(stim)),stim) + wv2(stim,4);
    Stimmat((oriind3(stim)-1)*8*8+(sfind3(stim)-1)*8+(Rphind3(stim)),stim)=...
        Stimmat((oriind3(stim)-1)*8*8+(sfind3(stim)-1)*8+(Rphind3(stim)),stim) + wv3(stim,4);
    Stimmat((oriind4(stim)-1)*8*8+(sfind4(stim)-1)*8+(Rphind4(stim)),stim)=...
        Stimmat((oriind4(stim)-1)*8*8+(sfind4(stim)-1)*8+(Rphind4(stim)),stim) + wv4(stim,4);
    %Big bad covariance matrix, remember we only have positive numbers
    Stimcovmat(:,:,stim)=Stimmat(:,stim)*Stimmat(:,stim)';
end

%% Framewise multiplication by rate observed
STAish=(Stimmat*Robs); %Independant response per dim, no covariance of waves
[what who]=max(mean(Robs));
% This is wrong, 'largest' dim is ori: order ph,sf,ori 
ch1=reshape(STAish(:,who),8,8,3);

imagesc(squeeze(ch1(:,:,1)))%.*repmat([1:8],8,1))
set(gca, 'YTick',1:8,'YTickLabels',binaround,'XTick',1:8,'XTickLabels',0.5*(1:8))
axis equal tight xy
xlabel('SF')
ylabel('Phase')

%%
goodunits=find(mean(meanRphasepairs,1)>1 & Exp.osp.clusterDepths>3000);
for cell=1:length(goodunits)
    ch_sta=reshape(STAish(:,goodunits(cell)),8,8,3);
    for jj=1:3
imagesc(squeeze(ch_sta(:,:,jj)))
set(gca, 'YTick',1:8,'YTickLabels',binaround,'XTick',1:8,'XTickLabels',0.5*(1:8))
axis equal tight xy
xlabel('SF')
ylabel('Phase')
title(num2str(goodunits(cell)))
pause(0.05)
    end
end

%% Covariance 
STCish=reshape(Stimcovmat,192^2,nstim)*Robs;
NC=size(Robs,2);
STCish=reshape(STCish,192,192,NC);

%%
imagesc(STCish(:,:,who))
[U,S,V] = svd(STCish(:,:,who)) ;
ch_stc1=reshape(U(:,1),8,8,3);
imagesc(squeeze(ch_stc1(:,:,1)))
imagesc(squeeze(ch_stc1(:,:,2)))
imagesc(squeeze(ch_stc1(:,:,3)))
% ch_stc2=reshape(U(:,2),8,8,3);
% imagesc(squeeze(ch_stc1(:,:,2)))
set(gca, 'YTick',1:8,'YTickLabels',binaround,'XTick',1:8,'XTickLabels',0.5*(1:8))
axis equal tight xy
xlabel('SF')
ylabel('Phase')


%% loop through each decent V2 cell
goodunits=find(mean(meanRphasepairs,1)>1 & Exp.osp.clusterDepths>3000);
for cell=1:length(goodunits)
    
    [U,S,V] = svd(STCish(:,:,goodunits(cell))) ;
    ch_stc1=reshape(U(:,1),8,8,3); %PC1
    for jj=1:3
    imagesc(squeeze(ch_stc1(:,:,jj)))
    set(gca, 'YTick',1:8,'YTickLabels',binaround,'XTick',1:8,'XTickLabels',0.5*(1:8))
    axis equal tight xy
    xlabel('SF')
    ylabel('Phase')
    title(num2str(goodunits(cell)))
    pause(0.05)
    end
end



%% Simpler: plot SF vs (absolute phase) for each orientation
% We really need to account for absolute phase if we want to say anything
% definitive, need to use eyetracking

% Every time point, and every wave need to be offset by a different number
% according to their SF and angle, remembering that 0 is a vertical line
% and that positive is to the right. Ori=0 only cares about x shift, x is
% in degrees, x/cpd = xcycles


%ori, sf, phase, dir-90 ,speed, contrast
%to just 
%sf, phase, ori, contrast
wv1=FullNoiseHistory(:,1+([2 3 4 6]));
wv2=FullNoiseHistory(:,1+6+([2 3 4 6]));
wv3=FullNoiseHistory(:,1+6*2+([2 3 4 6]));
wv4=FullNoiseHistory(:,1+6*3+([2 3 4 6]));
tic
%Loops are slow but easy to read for now
for stim=1:nstim
    %Closest previous call to the eyetracker
    eyeidx=find(wvtime(stim)<Exp.vpx.smo(:,1),1); %this takes forever
    eyex(stim)=Exp.vpx.smo(eyeidx,2);
    eyey(stim)=Exp.vpx.smo(eyeidx,3);

    %Dva *Cpd to cycles, multiply by 360 for phase(?)
    phshift1=wv1(stim,1)*(eyex(stim)*cosd(wv1(stim,3))+eyey(stim)*sind(wv1(stim,3)));
    wv1(stim,5)=wv1(stim,2)+phshift1*360;

    phshift2=wv2(stim,1)*(eyex(stim)*cosd(wv2(stim,3))+eyey(stim)*sind(wv2(stim,3)));
    wv2(stim,5)=wv2(stim,2)+phshift2*360;

    phshift3=wv3(stim,1)*(eyex(stim)*cosd(wv3(stim,3))+eyey(stim)*sind(wv3(stim,3)));
    wv3(stim,5)=wv3(stim,2)+phshift3*360;

    phshift4=wv4(stim,1)*(eyex(stim)*cosd(wv4(stim,3))+eyey(stim)*sind(wv4(stim,3)));
    wv4(stim,5)=wv4(stim,2)+phshift4*360;

    if mod(stim,round(nstim/10))==0
        disp([num2str((stim/nstim)*100,'%.2f') '% complete'])
    end
end
toc
%% Redo, accidently used wvx(stim,3) for phases
for stim=1:nstim
    

    %Dva *Cpd to cycles, multiply by 360 for phase(?)
    phshift1=wv1(stim,1)*(eyex(stim)*cosd(wv1(stim,3))+eyey(stim)*sind(wv1(stim,3)));
    wv1(stim,5)=wv1(stim,2)+phshift1*360;

    phshift2=wv2(stim,1)*(eyex(stim)*cosd(wv2(stim,3))+eyey(stim)*sind(wv2(stim,3)));
    wv2(stim,5)=wv2(stim,2)+phshift2*360;

    phshift3=wv3(stim,1)*(eyex(stim)*cosd(wv3(stim,3))+eyey(stim)*sind(wv3(stim,3)));
    wv3(stim,5)=wv3(stim,2)+phshift3*360;

    phshift4=wv4(stim,1)*(eyex(stim)*cosd(wv4(stim,3))+eyey(stim)*sind(wv4(stim,3)));
    wv4(stim,5)=wv4(stim,2)+phshift4*360;

    if mod(stim,round(nstim/10))==0
        disp([num2str((stim/nstim)*100,'%.2f') '% complete'])
    end
end




Exp.wvsave.wv1=wv1;
Exp.wvsave.wv2=wv2;
Exp.wvsave.wv3=wv3;
Exp.wvsave.wv4=wv4;



%%
ph1=wv1(:,5);
ph2=wv2(:,5);
ph3=wv3(:,5);
ph4=wv4(:,5);


    ph1=rem(ph1,360);
    ph1(ph1<0)=ph1(ph1<0)+360;
    ph2=rem(ph2,360);
    ph2(ph2<0)=ph2(ph2<0)+360;
    ph3=rem(ph3,360);
    ph3(ph3<0)=ph3(ph3<0)+360;
    ph4=rem(ph4,360);
    ph4(ph4<0)=ph4(ph4<0)+360;


%%
phbins=0:15:360;% go past 360 first

ph1_bins=((15)*round(ph1/15));
ph1_bins(ph1_bins==360)=0;
ph2_bins=((15)*round(ph2/15));
ph2_bins(ph2_bins==360)=0;
ph3_bins=((15)*round(ph3/15));
ph3_bins(ph3_bins==360)=0;
ph4_bins=((15)*round(ph4/15));
ph4_bins(ph4_bins==360)=0;

%sf wv1(:,1);
%ori wv(:,3);
sfset=unique(wv1(:,1));
phset=unique(ph1_bins);
phset=phset(~isnan(phset));
oriset=unique(wv1(:,3));

%% Rate per bin, mean 
Robs_conds1=zeros(length(oriset),length(sfset),length(phset),size(Robs,2));
for ori_ind=1:length(oriset)
    for sf_ind=1:length(sfset)
        for ph_ind=1:length(phset)
            cond_inds1=(wv1(:,3)==oriset(ori_ind)) & ...
                (wv1(:,1)==sfset(sf_ind)) & ...
                (ph1_bins==phset(ph_ind));
            % Changing order from above? ori,sf,ph
             Robs_conds1(ori_ind,sf_ind,ph_ind,:)=mean(Robs(cond_inds1,:),1);
        end
    end

end


%%
%% loop through cells, this is only 1/4 of the waves at a time, no covariance
goodunits=find(mean(Robs,1)>5); %& Exp.osp.clusterDepths>3000);
for cell=1:length(goodunits)
subplot(8,8,cell)
imagesc(squeeze(Robs_conds1(2,:,:,cell))'); axis xy
end





%% Now remake the STC
%% Big bad covariance matrix again, but this time we want to use the abs phases we found in phx_bins

% %Bin phases back to 8 starting phases
% binaround=0:45:315;%
% for ii=1:8
%     ind_bin1=abs(ph1_bins-binaround(ii))<15;
%     Aphase1(ind_bin1)=binaround(ii);
%     Aphind1(ind_bin1)=ii;
%     ind_bin2=abs(ph2_bins-binaround(ii))<15;
%     Aphase2(ind_bin2)=binaround(ii);
%     Aphind2(ind_bin2)=ii;
%     ind_bin3=abs(ph3_bins-binaround(ii))<15;
%     Aphase3(ind_bin3)=binaround(ii);
%     Aphind3(ind_bin3)=ii;
%     ind_bin4=abs(ph4_bins-binaround(ii))<15;
%     Aphase4(ind_bin4)=binaround(ii);
%     Aphind4(ind_bin4)=ii;
% end

%phases are already binned into 15 degree wide bins
Aphind1=ph1_bins/15+1;
Aphind2=ph2_bins/15+1;
Aphind3=ph3_bins/15+1;
Aphind4=ph4_bins/15+1;



%%
nstim=size(FullNoiseHistory,1);
%Build a stimulus space matrix, ori*sf*Aphase=1152 dims, for each frame
%6*8*24
%(oriind1-1)*8*24+(sfind1-1)*24+(Aphind1)
%ori indexes 3 blocks of 192, sfind into subblocks of 24 and so on
Stimmat=zeros(1152,nstim);
%Stimcovmat=zeros(1152,1152,nstim);
% Loop through each stimulus frame
% For each wave put contrast in the correct ori,sf and phase bin
for stim=1:nstim
    if ~isnan(Aphind1(stim))
    %Adding contrasts, in case of doubling up in bins
    Stimmat((oriind1(stim)-1)*24*8+(sfind1(stim)-1)*24+(Aphind1(stim)),stim)=...
        Stimmat((oriind1(stim)-1)*24*8+(sfind1(stim)-1)*24+(Aphind1(stim)),stim) + 1;%wv1(stim,4);
    Stimmat((oriind2(stim)-1)*24*8+(sfind2(stim)-1)*24+(Aphind2(stim)),stim)=...
        Stimmat((oriind2(stim)-1)*24*8+(sfind2(stim)-1)*24+(Aphind2(stim)),stim) + 1;%wv2(stim,4);
    Stimmat((oriind3(stim)-1)*24*8+(sfind3(stim)-1)*24+(Aphind3(stim)),stim)=...
        Stimmat((oriind3(stim)-1)*24*8+(sfind3(stim)-1)*24+(Aphind3(stim)),stim) + 1;%wv3(stim,4);
    Stimmat((oriind4(stim)-1)*24*8+(sfind4(stim)-1)*24+(Aphind4(stim)),stim)=...
        Stimmat((oriind4(stim)-1)*24*8+(sfind4(stim)-1)*24+(Aphind4(stim)),stim) + 1;%wv4(stim,4);
    %Big bad covariance matrix, remember we only have positive numbers
    %Stimcovmat(:,:,stim)=Stimmat(:,stim)*Stimmat(:,stim)';
    end
end



%% Framewise multiplication by rate observed
STAish=(Stimmat*(Robs-mean(Robs))); %Independant response per dim, no covariance of waves
[what who]=max(mean(Robs));
% This is wrong, 'largest' dim is ori: order ph,sf,ori 
ch1=reshape(STAish(:,who),24,8,6);

imagesc(squeeze(ch1(:,:,1)))%.*repmat([1:8],8,1))
% set(gca, 'YTick',1:8,'YTickLabels',binaround,'XTick',1:8,'XTickLabels',0.5*(1:8))
% axis equal tight xy
% xlabel('SF')
% ylabel('Phase')

%%
goodunits=find(mean(Robs)>1 & Exp.osp.clusterDepths>3000);
for cell=1:length(goodunits)
    ch_sta=reshape(STAish(:,goodunits(cell)),24,8,6);
    for jj=1:3
imagesc(squeeze(ch_sta(:,:,jj)))
set(gca, 'YTick',1:24,'YTickLabels',phset,'XTick',1:8,'XTickLabels',0.5*(1:8))
axis equal tight xy
xlabel('SF')
ylabel('Phase')
title(num2str(goodunits(cell)))
pause(0.5)
    end
end

%% Covariance, unfortunately has to be the slow way

NC=size(Robs,2);
Robs2=reshape(Robs-mean(Robs),nstim,1,NC);
STCish=zeros(1152,1152,NC);
for stim=1:nstim
STCish=STCish+Stimmat(:,stim)*Stimmat(:,stim)'.*Robs2(stim,:,:);
end



%%
imagesc(STCish(:,:,who))
[U,S,V] = svd(STCish(:,:,who)) ;
ch_stc1=reshape(U(:,1),24,8,6);
imagesc(squeeze(ch_stc1(:,:,1)))
imagesc(squeeze(ch_stc1(:,:,2)))
imagesc(squeeze(ch_stc1(:,:,3)))
% ch_stc2=reshape(U(:,2),8,8,3);
% imagesc(squeeze(ch_stc1(:,:,2)))
set(gca, 'YTick',1:24,'YTickLabels',phset,'XTick',1:8,'XTickLabels',0.5*(1:8))
axis equal tight xy
xlabel('SF')
ylabel('Phase')


%% loop through each decent V2 cell
goodunits=find(mean(Robs)>1 );
for cell=1:length(goodunits)
    
    [U,S,V] = svd(STCish(:,:,goodunits(cell))) ;
    ch_stc1=reshape(U(:,1),24,8,6); %PC1
    ch_stc2=reshape(U(:,2),24,8,6); %PC1
    ch_stc3=reshape(U(:,3),24,8,6); %PC1
    for jj=1
        subplot(131)
        imagesc(squeeze(ch_stc1(:,:,jj)))
        set(gca, 'YTick',1:24,'YTickLabels',phset,'XTick',1:8,'XTickLabels',0.5*(1:8))
        axis equal tight xy
        xlabel('SF')
        ylabel('Phase')
        title([num2str(goodunits(cell)) ' PC1'])
        subplot(132)
        imagesc(squeeze(ch_stc2(:,:,jj)))
        set(gca, 'YTick',1:24,'YTickLabels',phset,'XTick',1:8,'XTickLabels',0.5*(1:8))
        axis equal tight xy
        xlabel('SF')
        ylabel('Phase')
        title([num2str(goodunits(cell)) ' PC2'])
        subplot(133)
        imagesc(squeeze(ch_stc3(:,:,jj)))
        set(gca, 'YTick',1:24,'YTickLabels',phset,'XTick',1:8,'XTickLabels',0.5*(1:8))
        axis equal tight xy
        xlabel('SF')
        ylabel('Phase')
        title([num2str(goodunits(cell)) ' PC3'])
    pause(0.5)
    end
end












%% ---------%%

%% fake strobe times (necessary to synchronize the eye tracker)
for iTrial = 1:numel(Exp.D)

    t = Exp.ptb2Ephys(Exp.D{iTrial}.STARTCLOCKTIME);
    if t > 0
        Exp.D{iTrial}.START_EPHYS = t;
    end

    t = Exp.ptb2Ephys(Exp.D{iTrial}.ENDCLOCKTIME);
    if t > 0
        Exp.D{iTrial}.END_EPHYS = t;
    end
end


%% Eyetracking calibration
% Load facecall data
dataPath1 ='/mnt/NPX/Rocky/20230921/Rocky 2023-09-21/';
%EyeExp = io.basic_marmoview_import(dataPath1);
fid = 1;
% [Exp, fig] = io.import_eye_position(Exp, dataPath1, 'fid', fid, 'zero_mean', true, 'normalize', true);
% 

[Exp, fig] = io.import_eye_position(Exp, dataPath1, 'fid', fid, 'zero_mean', true, 'normalize', true);
%[Exp, fig] = io.import_eye_position_ddpi09072023(Exp, dataPath1, 'fid', fid, 'zero_mean', true, 'normalize', true);

%%
%%
C = calibGUI(Exp,'use_smo','true');
C.cmat
%% Saving Cmat so I don't need to refine again

cmat=[0.9828    0.9222   -2.0009   -0.2875   -0.0349]; % online track
%cmat=[0.1429    0.1401   -1.9981  -39.9112   20.2618]; %Offline dpi
C.cmat=cmat;
%%
% eyePos = EyeExp.vpx.smo(:,2:3);
% eyePos =get_eyepos(C);
% fig = io.checkCalibration(EyeExp, eyePos);
% %

%% And data set
% fid=1;
% [Exp, fig] = io.import_eye_position(Exp, StimPath, 'fid', fid, 'zero_mean', true, 'normalize', true);
% eyePosRaw = Exp.vpx.smo(:,2:3);
% eyePosRaw = Exp.vpx.raw0(:,2:3);
% 
% 
%             th = cmat(3);
%             R = [cosd(th) -sind(th); sind(th) cosd(th)];
%             S = [cmat(1) 0; 0 cmat(2)];
%             A = (R*S)';
%             Ainv = pinv(A);
% 
%             eyePos = (Exp.vpx.raw0(:,2:3) - cmat(4:5))*Ainv;
%             
% 
% %Why doesn't this work???
% fig = io.checkCalibration(Exp, eyePos);
% % Now it works?

%%
eyePos =C.get_eyepos();


%% Update the eyePos in Exp (Overwriting the raw data feels like the wrong way to do this, but subsequent analysis expect it?)
Exp.vpx.smo(:,2:3)=eyePos;

%% Update to raw eyePos in Exp
C = calibGUI(Exp,'use_smo',0);
%% Raw data is way different
 cmat=[11.2678   10.7247   -3.3977  326.9572 -233.3391];%online t
%cmat=[7.0738    7.4387   -0.8886  -91.1588  -60.6488];%offline dpi
C.cmat=cmat
%%
eyePos =C.get_eyepos();
io.checkCalibration(Exp, eyePos);
%% Update
Exp.vpx.raw(:,2:3)=eyePos;
%%

% Saccade processing:
% Perform basic processing of eye movements and saccades
Exp = saccadeflag.run_saccade_detection_cloherty(Exp, ...
    'ShowTrials', false,...
    'accthresh', 2e4,...
    'velthresh', 10,...
    'velpeak', 10,...
    'isi', 0.02);
%
% track invalid sampls
Exp.vpx.Labels(isnan(Exp.vpx.smo(:,2))) = 4;

validTrials = io.getValidTrials(Exp, 'Grating');
for iTrial = validTrials(:)'
    if ~isfield(Exp.D{iTrial}.PR, 'frozenSequence')
        Exp.D{iTrial}.PR.frozenSequence = false;
    end
end
        
fprintf(fid, 'Done importing session\n');


%% Load spikes
dataPath = '/home/huklab/Documents/SpikeSorting/Output/2023-09-21/';
Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp.st = Exp.osp.st - .025;
Exp.osp.st(Exp.osp.st<0)=0;
% move zero to the end
id = max(Exp.osp.clu)+1;

Exp.osp.clu(Exp.osp.clu==0) = id;
Exp.osp.cids = Exp.osp.cids2;
Exp.osp.cids(Exp.osp.cids==0) = id;

% lbledgood=find(Exp.osp.cgs==2)+1;

% Exp.osp.st = st_s.spikeTimes(:,1)+.6;
% Exp.osp.clu = st_s.spikeTimes(:,2);
% Exp.osp.depths = st_s.depths;
% Exp.osp.templates = st_s.templates;
% Exp.osp.cids = unique(Exp.osp.clu);


figure(1); clf
plot.raster(Exp.osp.st, Exp.osp.clu, 1)
ylabel('Unit ID')
xlabel('Time (seconds)')

%% Get drifting gratings
%validTrials = io.getValidTrials(Exp, 'DriftingGrating');



%%
figure(1); clf
plot.raster(Exp.osp.st, Exp.osp.clu, 1); hold on
ylabel('Unit ID')
xlabel('Time (seconds)')
cmap = lines;
stims = {'Dots', 'DriftingGrating', 'Gabor', 'BackImage', 'FaceCal','FixRsvpStim'};

for istim = 1:numel(stims)
    validTrials = io.getValidTrials(Exp, stims{istim});
    for iTrial = validTrials(:)'
        t_ = Exp.ptb2Ephys(Exp.D{iTrial}.STARTCLOCKTIME);
        t1 = Exp.ptb2Ephys(Exp.D{iTrial}.ENDCLOCKTIME);
        yd = ylim;
        h(istim) = fill([t_ t_ t1 t1], [yd(1) yd(2) yd(2) yd(1)], 'k', 'FaceColor', cmap(istim,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
    end
end

legend(h, stims)

%% Get spatial map
dotTrials = io.getValidTrials(Exp, 'Dots');
dotTrials=dotTrials(dotTrials<1518);
Exp_subset=Exp;
%Exp_subset.D={Exp.D{dotTrials}}';
if ~isempty(dotTrials)
    
    BIGROI = [-15 -10 5 10];
    %eyePos = C.refine_calibration();
    binSize = .2;
    Frate = 60;
    [Xstim, RobsSpace, opts] = io.preprocess_spatialmapping_data(Exp_subset, ...
        'ROI', BIGROI*Exp.S.pixPerDeg, 'binSize', binSize*Exp.S.pixPerDeg, ...
        'eyePosExclusion', 2e3, ...
        'eyePos', Exp.vpx.smo(:,2:3), 'frate', Frate, 'fastBinning', true);
    
    % use indices while fixating
    ecc = hypot(opts.eyePosAtFrame(:,1), opts.eyePosAtFrame(:,2))/Exp.S.pixPerDeg;
    ix = opts.eyeLabel==1 & ecc < 15.2;
    
end


spike_rate = mean(RobsSpace)*Frate;

%% Get spatial map from only static dots 
dotTrials = io.getValidTrials(Exp, 'StaticDots');

Exp_subset=Exp;
%Exp_subset.D={Exp.D{dotTrials}}';
if sum(cellfun(@(x) strcmp(x.PR.name,'FixedProceduralNoise'),Exp.D))
    
    BIGROI = [-15 -10 5 10];
    %eyePos = C.refine_calibration();
    binSize = .2;
    Frate = 60;
    [Xstim, RobsSpace, opts] = io.preprocess_staticspatialmapping_data(Exp_subset, ...
        'ROI', BIGROI*Exp.S.pixPerDeg, 'binSize', binSize*Exp.S.pixPerDeg, ...
        'eyePosExclusion', 2e3, ...
        'eyePos',Exp.vpx.fix, 'frate', Frate, 'fastBinning', true);% Exp.vpx.smo(:,2:3)
    
    % use indices while fixating
    ecc = hypot(opts.eyePosAtFrame(:,1), opts.eyePosAtFrame(:,2))/Exp.S.pixPerDeg;
    ix = opts.eyeLabel==1 & ecc < 15.2;
    
end


spike_rate = mean(RobsSpace)*Frate;

%%


figure(2); clf; set(gcf, 'Color', 'w')
subplot(3,1,1:2)
h = [];
h(1) = stem(spike_rate, '-ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4);

ylabel('Firing Rate During Stimulus')
hold on
goodunits = find(spike_rate > .1);
h(2) = plot(goodunits, spike_rate(goodunits), 'og', 'MarkerFaceColor', 'g', 'MarkerSize', 4);
h(3) = plot(xlim, [1 1], 'r-');
legend(h, {'All units', 'Good Units', '.1 Spike/Sec'}, 'Location', 'Best')

subplot(3,1,3)
stem(spike_rate, '-ok', 'MarkerFaceColor', 'k', 'MarkerSize', 4)
ylim([0, 1])
xlabel('Unit #')
ylabel('Firing rate of bad units')


fprintf('%d / %d fire enough spikes to analyze\n', numel(goodunits), size(RobsSpace,2))
drawnow
%
win = [-3 15];
stas = forwardCorrelation(Xstim, mean(RobsSpace,2), win);
stas = stas / std(stas(:)) - mean(stas(:));
wm = [min(stas(:)) max(stas(:))];
nlags = size(stas,1);
%
figure(1); clf
for ilag = 1:nlags
    subplot(2, ceil(nlags/2), ilag)
    imagesc(opts.xax, opts.yax, reshape(stas(ilag, :), opts.dims), wm)
end
%%
figure(1); clf
for ilag = 8
    subplot(1,1,1)
    xax = opts.xax/Exp.S.pixPerDeg;
    yax = opts.yax/Exp.S.pixPerDeg;
    imagesc(xax, yax, reshape(stas(ilag, :), opts.dims), wm)
end
axis equal tight xy
hold on
        scatter(0,0,'k+');
        hold off
%% Analyze population
R = RobsSpace(:,goodunits);
%stas = forwardCorrelation(Xstim, R-mean(R), win, find(ix));
stas = forwardCorrelation(Xstim, R-mean(R), win);
NC = size(R,2);

%% plotting
sx = 5;
sy = 5; %ceil(NC/sx);
nperfig = sx*sy;
nfigs = ceil(NC / nperfig);
fprintf('%d figs\n', nfigs)

if nfigs == 1
    figure(1); clf; set(gcf, 'Color', 'w')
end

for cc = 1:NC
    if nfigs > 1
        fig = ceil(cc/nperfig);
        figure(fig); set(gcf, 'Color', 'w')
    end

    si = mod(cc, nperfig);  
    if si==0
        si = nperfig;
    end
    subplot(sy,sx, si)

    I = squeeze(stas(:,:,cc));
    I = (I - mean(I(:)) )/ std(I(:));

    [bestlag,j] = find(I==max(I(:)));
    I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
    I = mean(I);
%     I = std(I);
    xax = opts.xax/Exp.S.pixPerDeg;
    yax = opts.yax/Exp.S.pixPerDeg;
    imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2]);
    axis equal tight
%    set(gca, 'XTick',[],'YTick',[])
    colormap(plot.coolwarm)
     title(Exp.osp.cids(goodunits(cc)))
end
    

%% plot subset
subset = [33 152 156 168 171 187 214 253 256 258 259 262 264 266 318 319 ...
    327 365 383 396 400 423 448 454 460 463 471 473 493 536 688 721 722 776 ...
    875 1001 1032 1034 1035 1174 1289 1396];
nsub = numel(subset);
[~, ccid] = ismember(subset, Exp.osp.cids(goodunits));
%[~, ccid] = ismember(subset, Exp.osp.cids);

nfigs = ceil(nsub / nperfig);
fig = ceil(1); clf
figure(fig); set(gcf, 'Color', 'w')

sx = ceil(sqrt(nsub));
sy = round(sqrt(nsub));

for cc = 1:nsub
    subset(cc)
%     if nfigs > 1
%         fig = ceil(cc/nperfig);
%         figure(fig); set(gcf, 'Color', 'w')
%     end
% 
%     si = mod(cc, nperfig);  
%     if si==0
%         si = nperfig;
%     end
    si=nsub-cc+1;
    subplot(sy,sx, si)
    if ccid(cc)>0
        I = squeeze(stas(:,:,ccid(cc)));
        I = (I - mean(I(:)) )/ std(I(:));
    
        [bestlag,j] = find(I==max(I(:)));
        I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
        I = mean(I);
    %     I = std(I);
        xax = opts.xax/Exp.S.pixPerDeg;
        yax = opts.yax/Exp.S.pixPerDeg;
        imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2]); axis xy 
        hold on
        scatter(0,0,'k+');
        colormap(plot.coolwarm)
        title(Exp.osp.cids(goodunits(ccid(cc))))

        %NOTE: depth is -ve (distance up from tip), IE deepest unit would
        %be 0
        title(5000-Exp.osp.clusterDepths(goodunits(ccid(cc))))
        grid on
    end
end


%% Grating tuning prep
Exp.spikeTimes = Exp.osp.st;
Exp.spikeIds = Exp.osp.clu;

nTrials=numel(Exp.D);
for ii=1:nTrials
    Exp.D{ii}.treadmill=Exp.D{ii}.inputs{2};
end

%%

% convert to D struct format
D = io.get_drifting_grating_output(Exp);

% Note: PLDAPS sessions and MarmoV5 will load into different formats

%% Step 2.1: plot one session
% This just demonstrates how to index into individual sessions and get
% relevant data

sessionId = 1; % pick a session id

% index into grating times from that session
gratIx = D.sessNumGratings==sessionId;
gratIx = gratIx & D.GratingOnsets > min(D.spikeTimes) & D.GratingOffsets < max(D.spikeTimes);
D.sessNumGratings(~gratIx) = 10;
gratOnsets = D.GratingOnsets(gratIx);
gratOffsets = D.GratingOffsets(gratIx);
gratDir = D.GratingDirections(gratIx);




%%

% find treadmill times that correspond to this session 
treadIx = D.treadTime > gratOnsets(1) & D.treadTime < gratOffsets(end);
treadTime = D.treadTime(treadIx);

%% PLOT SPIKE TIMES
spikeIds = unique(D.spikeIds(D.sessNumSpikes==sessionId));
NC = numel(spikeIds);
spikeRate = zeros(numel(gratOnsets), NC);

%bs = diff(treadTime);
gs=diff(gratOnsets);
for cc = 1:NC
    spikeRate(:,cc) = [0 histcounts(D.spikeTimes(D.spikeIds==spikeIds(cc)), gratOnsets)./gs'];
end

spikeRate = spikeRate ./ max(spikeRate); % normalize for visualization



%%
figure(111); clf
% PLOT GRATINGS TIME VS. DIRECTION
subplot(3,1,1)
nGratings = numel(gratOnsets);
for iG = 1:nGratings
    plot([gratOnsets(iG) gratOffsets(iG)], gratDir(iG)*[1 1], 'r', 'Linewidth', 2); hold on
end
ylabel('Direction')
xlim(treadTime([1 end]))
title('Gratings')
ylim([0 360])


%%

for dir=1:numel(unique(gratDir))
spikeCond(dir,:)=mean(spikeRate(gratDir==30*dir,:));
end
for cc=1:500
plot(spikeCond(:,cc))
pause(0.1)
end

%% PLOT SPIKE TIMES
spikeIds = unique(D.spikeIds(D.sessNumSpikes==sessionId));
NC = numel(spikeIds);
spikeRate = zeros(numel(treadTime), NC);

bs = diff(treadTime);
for cc = 1:NC
    spikeRate(:,cc) = [0 histcounts(D.spikeTimes(D.spikeIds==spikeIds(cc)), treadTime)./bs'];
end

spikeRate = spikeRate ./ max(spikeRate); % normalize for visualization

treadSpeed = D.treadSpeed(treadIx);
runThresh=3;

isrunning = treadSpeed > runThresh;
onsets = find(diff(isrunning) == 1);
offsets = find(diff(isrunning) == -1);

if onsets(1) > offsets(1)
    onsets = [1; onsets];
end

if offsets(end) < onsets(end)
    offsets = [offsets; numel(treadSpeed)];
end

assert(numel(onsets)==numel(offsets), "onset offset mismatch")

subplot(3,1,2) % plot spike count
imagesc(treadTime, 1:NC, spikeRate'); hold on
for i = 1:numel(onsets)
    fill(treadTime([onsets(i) onsets(i) offsets(i) offsets(i)]), [ylim fliplr(ylim)], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
end
xlim(treadTime([1 end]))
colormap(1-gray)
title('Spikes')
ylabel('Unit #')   

% % PLOT TREADMILL RUNNING SPEED
% subplot(3,1,3) % tread speed
% plot(treadTime, treadSpeed , 'k'); hold on
% xlabel('Time (s)')
% ylabel('Speed (cm/s)')
% 
% for i = 1:numel(onsets)
%     fill(treadTime([onsets(i) onsets(i) offsets(i) offsets(i)]), [ylim fliplr(ylim)], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
% end
% xlim(treadTime([1 end]))
% 
% title('Treadmill')
% trmax = prctile(D.treadSpeed(treadIx), 99);
% 
% plot.fixfigure(gcf, 12, [8 11])
% saveas(gcf, fullfile(figDir, 'grating_sess_birdseye.pdf'))
% 
% %% play movie
% vidObj = VideoWriter(fullfile(figDir, 'session.mp4'), 'MPEG-4');
% vidObj.Quality = 100;
% vidObj.FrameRate = 5;
% open(vidObj);
% 
% 
% % play a movie of the session
% figure(gcf)
% t0 = treadTime(1);
% win = 25;
% for t = 1:500
%     xd = [t0 t0 + win];
%     for i = 1:3
%         subplot(3,1,i)
%         xlim(xd)
%     end
%     ylim([0 trmax])
%     drawnow
%     t0 = t0 + win/10;
%     if t0 > treadTime(end)
%         break
%     end
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
% end
% 
% close(vidObj);
% 
% %% Firing rates
% D.treadSpeed = max(D.treadSpeed, 0);
% D.sessNumTread = double(D.treadTime > 0);
% 
% iix = ismember(D.spikeIds, Exp.osp.cids(goodunits));
% D.spikeIds(~iix) = 9999;
% 
% %%
% sessionId = 1;
% [StimDir, spksb, runningSpeed, Dstat] = bin_population(D, sessionId);
% 
% %% PSTH
% stimids = unique(StimDir);
% nstim = numel(stimids);
% 
% Rtuning = zeros(nstim, Dstat.NLags, Dstat.NCells);
% for istim = 1:nstim
%     iix = StimDir==stimids(istim);
%     Rtuning(istim,:,:) = squeeze(mean(spksb(iix,:,:)));
% end
% 
% % running?
% isrunning = mean(runningSpeed,2) > 3;
% Rtc = zeros(nstim, 2, Dstat.NCells);
% tix = Dstat.lags > .05 & Dstat.lags < .9;
% R = squeeze(mean(spksb(:,tix,:), 2));
% for istim = 1:nstim
%     iix = StimDir==stimids(istim);
% 
%     Rtc(istim,1,:) = mean(R(iix & isrunning,:));
%     Rtc(istim,2,:) = mean(R(iix & ~isrunning,:));
% end
% 
% %% check out individual cell
% %subset = [31,32,108,144,358, 598, 648, 892, 692,905,918,934, 949, 964, 965, 980, 1174, 1183, 1055];
% subset = [6, 31,34,37,44,107,114,124,143,147,164,227,276,384, 488, 584, 630, 761, 805,819];
% 
% nsub = numel(subset);
% [~, ccid] = ismember(subset, Exp.osp.cids(goodunits));
% 
% %%
% cc = cc + 1;
%     
% I = squeeze(stas(:,:,ccid(cc)));
% I = (I - mean(I(:)) )/ std(I(:));
% 
% [bestlag,j] = find(I==max(I(:)));
% I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
% I = mean(I);
% 
% figure(1); clf
% subplot(2,2,1)
% xax = opts.xax/Exp.S.pixPerDeg;
% yax = opts.yax/Exp.S.pixPerDeg;
% imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2])
% colormap(plot.coolwarm)
% title(Exp.osp.cids(goodunits(cc)))
% 
% i = find(I == max(I(:)));
% subplot(2,2,2)
% tkernel = stas(:, i, ccid(cc));
% plot(win(1):win(2), tkernel, 'k', 'Linewidth', 2)
% 
% subplot(2,2,3)
% for istim = 1:nstim
%     plot(Dstat.lags, imboxfilt(squeeze(Rtuning(istim,:,ccid(cc))),11)*60, 'Color', cmap(istim,:)); hold on
% end
% 
% subplot(2,2,4)
% plot(stimids, Rtc(:,1,ccid(cc))*60, '-o'); hold on
% plot(stimids, Rtc(:,2,ccid(cc))*60, '-o');
% 
% 
% %%
% %% plotting
% sx = 15;
% NC = size(stas,3);
% sy = ceil(NC/sx);
% 
% figure(1); clf; set(gcf, 'Color', 'w')
% ax = plot.tight_subplot(sy, sx, 0.01, 0.01);
% 
% for cc = 1:NC
% 
%     set(gcf, 'currentaxes', ax(cc))
% 
%     I = squeeze(stas(:,:,cc));
%     I = (I - mean(I(:)) )/ std(I(:));
% 
%     [bestlag,j] = find(I==max(I(:)));
%     I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
%     I = mean(I);
%     xax = opts.xax/Exp.S.pixPerDeg;
%     yax = opts.yax/Exp.S.pixPerDeg;
%     imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2])
%     colormap(plot.coolwarm)
%     title(Exp.osp.cids(goodunits(cc)))
% end
%     
% 
% 
% %%
% figure(2); clf
% 
% cmap = plot.viridis(nstim);
% ax = plot.tight_subplot(sy, sx, 0.01, 0.01);
% 
% for cc = 1:Dstat.NCells
%     set(gcf, 'currentaxes', ax(cc))
%     for istim = 1:nstim
%         plot(Dstat.lags, imboxfilt(squeeze(Rtuning(istim,:,cc)),11), 'Color', cmap(istim,:)); hold on
%     end
%     axis off
% end
% 
% 
% figure(3); clf
% ax = plot.tight_subplot(sy, sx, 0.01, 0.01);
% for cc = 1:Dstat.NCells
%     set(gcf, 'currentaxes', ax(cc))
%     imagesc(Dstat.lags, 1:nstim, squeeze(Rtuning(:,:,cc)))
%     axis off
% end
% 
% figure(4); clf
% ax = plot.tight_subplot(sy, sx, 0.01, 0.01);
% for cc = 1:Dstat.NCells
%     set(gcf, 'currentaxes', ax(cc))
%     plot(stimids, Rtc(:,1,cc)*60, '-o'); hold on
%     plot(stimids, Rtc(:,2,cc)*60, '-o');
% end
% %%
% 
% 
% 
% %% Unit by Unit simple analysis
% thresh = 3;
% nboot = 100;
% 
% unitList = unique(D.spikeIds);
% NC = numel(unitList);
% 
% corrRho = zeros(NC,1);
% corrPval = zeros(NC,1);
% 
% frBaseR = zeros(NC,3);
% frBaseS = zeros(NC,3);
% 
% frStimR = zeros(NC,3);
% frStimS = zeros(NC,3);
% 
% for cc = 1:NC
%     unitId = unitList(cc);
% 
%     [stimDir0, robs, runSpd, opts] = bin_ssunit(D, unitId, 'win', [-.2 .1]);
%     
% %     pause
%     goodIx = getStableRange(sum(robs,2), 'plot', false);
%     
%     %stimDir0 = stimDir0(goodIx);
%     stimDir = stimDir0{1};
%     stimDir = stimDir(goodIx);
%     
%     robs = robs(goodIx,:);
%     runSpd = runSpd{1}(goodIx,:);
%     
%     iix = opts.lags < 0;
%     frbase = sum(robs(:,iix),2) / (max(opts.lags(iix)) - min(opts.lags(iix)));
%     
%     spd = mean(runSpd,2);
%     
%     [corrRho(cc), corrPval(cc)] = corr(spd, frbase, 'type', 'Spearman');
%     
%     runTrials = find(spd > thresh);
%     statTrials = find(abs(spd) < 1);
%     mixTrials = [runTrials; statTrials];
%     
%     nrun = numel(runTrials);
%     nstat = numel(statTrials);
%     
%     n = min(nrun, nstat);
%     
%     frBaseR(cc,:) = prctile(mean(frbase(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
%     frBaseS(cc,:) = prctile(mean(frbase(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
%     
%     iix = opts.lags > 0.04 & opts.lags < opts.lags(end)-.15;
%     frstim = sum(robs(:,iix),2) / (max(opts.lags(iix)) - min(opts.lags(iix)));
%     
%     frStimR(cc,:) = prctile(mean(frstim(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
%     frStimS(cc,:) = prctile(mean(frstim(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
%     
% end
% 
% %% plot some outcomes
% 
% incBaseIx = find(frBaseR(:,2) > frBaseS(:,3));
% decBaseIx = find(frBaseR(:,2) < frBaseS(:,1));
% 
% incStimIx = find(frStimR(:,2) > frStimS(:,3));
% decStimIx = find(frStimR(:,2) < frStimS(:,1));
% 
% figure(1); clf
% set(gcf, 'Color', 'w')
% ms = 4;
% cmap = lines;
% subplot(1,2,1)
% plot(frBaseS(:,[2 2])', frBaseR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
% plot(frBaseS(:,[1 3])', frBaseR(:,[2 2])', 'Color', .5*[1 1 1])
% plot(frBaseS(:,2), frBaseR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
% plot(frBaseS(incBaseIx,2), frBaseR(incBaseIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
% plot(frBaseS(decBaseIx,2), frBaseR(decBaseIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))
% 
% xlabel('Stationary')
% ylabel('Running')
% title('Baseline Firing Rate')
% 
% set(gca, 'Xscale', 'log', 'Yscale', 'log')
% plot(xlim, xlim, 'k')
% 
% subplot(1,2,2)
% plot(frStimS(:,[2 2])', frStimR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
% plot(frStimS(:,[1 3])', frStimR(:,[2 2])', 'Color', .5*[1 1 1])
% plot(frStimS(:,2), frStimR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
% plot(frStimS(incStimIx,2), frStimR(incStimIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
% plot(frStimS(decStimIx,2), frStimR(decStimIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))
% 
% xlabel('Stationary Firing Rate')
% ylabel('Running Firing Rate')
% title('Stim-driven firing rate')
% 
% set(gca, 'Xscale', 'log', 'Yscale', 'log')
% plot(xlim, xlim, 'k')
% 
% nIncBase = numel(incBaseIx);
% nDecBase = numel(decBaseIx);
% 
% nIncStim = numel(incStimIx);
% nDecStim = numel(decStimIx);
% 
% modUnits = unique([incBaseIx; decBaseIx; incStimIx; decStimIx]);
% nModUnits = numel(modUnits);
% 
% fprintf('%d/%d (%02.2f%%) increased baseline firing rate\n', nIncBase, NC, 100*nIncBase/NC)
% fprintf('%d/%d (%02.2f%%) decreased baseline firing rate\n', nDecBase, NC, 100*nDecBase/NC)
% 
% fprintf('%d/%d (%02.2f%%) increased stim firing rate\n', nIncStim, NC, 100*nIncStim/NC)
% fprintf('%d/%d (%02.2f%%) decreased stim firing rate\n', nDecStim, NC, 100*nDecStim/NC)
% 
% fprintf('%d/%d (%02.2f%%) total modulated units\n', nModUnits, NC, 100*nModUnits/NC)
% 
% [pvalStim, ~, sStim] = signrank(frStimS(:,2), frStimR(:,2));
% [pvalBase, ~, sBase] = signrank(frBaseS(:,2), frBaseR(:,2));
% 
% fprintf('Wilcoxon signed rank test:\n')
% fprintf('Baseline rates: p = %02.10f\n', pvalBase)
% fprintf('Stim-driven rates: p = %02.10f\n', pvalStim)
% 
% good = ~(frBaseR(:,2)==0 | frBaseS(:,2)==0);
% 
% m = geomean(frBaseR(good,2)./frBaseS(good,2));
% ci = bootci(nboot, @geomean, frBaseR(good,2)./frBaseS(good,2));
% 
% fprintf("geometric mean baseline firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good)) 
% 
% m = geomean(frStimR(:,2)./frStimS(:,2));
% % ci = bootci(nboot, @geomean, frStimR(:,2)./frStimS(:,2));
% ci = nan(1,2);
% fprintf("geometric mean stim-driven firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)
% 
% 
% 
% 
% %% Do direction / orientation decoding
% 
% 
% % example session to explore
% Dstat = decode_stim(D, 1, 'figDir', figDir);
% 
% %% Bootstrapped empirical analyses and tuning curve fits
% 
% fpath = getpref('FREEVIEWING', 'HUKLAB_DATASHARE');
% 
% % D.sessNumTread = ones(size(D.treadTime))*1;
% 
% D.treadSpeed(D.treadSpeed < 0) = nan;
% stat = tuning_empirical(D, 'binsize', 10e-3, ...
%     'runningthresh', 3, ...
%     'nboot', 500, ...
%     'seed', 1234);  
% 
% 
% 
% %% Plot all tuning curves
% fitS = stat.TCfitS;
% fitR = stat.TCfitR;
% NC = numel(fitS);
% 
% sx = ceil(sqrt(NC));
% sy = round(sqrt(NC));
% 
% 
% figure(1); clf
% ax = plot.tight_subplot(sx, sy, 0.02);
% for cc = 1:NC
%     if min(fitS(cc).numTrials, fitR(cc).numTrials) < 50
%         continue
%     end
% %     
% %     if fitS(cc).llrpval > 0.05 && fitR(cc).llrpval > 0.05
% %         continue
% %     end
%     fprintf("Unit %d/%d\n", cc, NC)
% 
%     set(gcf, 'currentaxes', ax(cc))
%     
%     cmap = lines;
%     % STATIONARY
%     h = errorbar(fitS(cc).thetas, fitS(cc).tuningCurve, fitS(cc).tuningCurveSE, '-o', 'Color', cmap(1,:));
%     h.MarkerSize = 2;
%     h.MarkerFaceColor = cmap(1,:);
%     h.CapSize = 0;
%     hold on
%     if fitS(cc).llrpval > 0.05
%         plot(xlim, mean(fitS(cc).tuningCurve)*[1 1], 'Color', cmap(1,:))
%     else
%         plot(linspace(0, 360, 100), fitS(cc).tuningFun(linspace(0, 360, 100)), 'Color', cmap(1,:))
%     end
%     
%     % RUNNING
%     h = errorbar(fitR(cc).thetas, fitR(cc).tuningCurve, fitR(cc).tuningCurveSE, '-o', 'Color', cmap(2,:));
%     h.MarkerSize = 2;
%     h.MarkerFaceColor = cmap(2,:);
%     h.CapSize = 0;
%     hold on
%     if fitR(cc).llrpval > 0.05
%         plot(xlim, mean(fitR(cc).tuningCurve)*[1 1], 'Color', cmap(2,:))
%     else
%         plot(linspace(0, 360, 100), fitR(cc).tuningFun(linspace(0, 360, 100)), 'Color', cmap(2,:))
%     end
% 
%     set(gca, 'XTick', [], 'YTick', [])
% %     set(gca, 'XTick', 0:180:360)
%     xlim([0 360])
%     text(10, .1*max(ylim), sprintf('%d', cc), 'fontsize', 5)
% %     axis off
%     title(sprintf('Unit: %d', Exp.osp.cids(cc)))
% end
% 
% plot.fixfigure(gcf, 10, [sx sy]*2, 'offsetAxes', false);
% saveas(gcf, fullfile(figDir, 'tuning_curves.pdf'))
% 
% %% If happy, save processe file
% % pathparts = regexp(dataPath, '/', 'split');
% % Exp.FileTag = [pathparts{end-1} '_' pathparts{end}(2:end) '.mat'];
% Exp.FileTag = 'gru/g20220412.mat'
% save(fullfile(HUKDATASHARE, Exp.FileTag), '-v7.3', '-struct', 'Exp')


%%
inds = io.getValidTrials(Exp, 'FixRsvpStim');

spikeIds=unique(Exp.osp.clu);
NC=length(spikeIds);
spikeRate = zeros(111, NC, length(inds));
most=mode(Exp.osp.clu);
hero=find(spikeIds==most);

figure(1); clf

%%
for trial = 1:length(inds)
    thisTrial=inds(trial);
    if length(Exp.D{thisTrial}.PR.NoiseHistory(:,4))>0
    tt=Exp.D{thisTrial}.PR.NoiseHistory(:,1);
    plot( tt-tt(1) ,Exp.D{thisTrial}.PR.NoiseHistory(:,4)) ; hold on
    
    %Convert start time (In PTB time) to ephys time
    t_=Exp.ptb2Ephys(Exp.D{thisTrial}.PR.NoiseHistory(1,1)); 
    
    t0=t_-10;
    t1=t_+10;
    idx=Exp.osp.st>t0&Exp.osp.st<t1;
    subset_st=Exp.osp.st(idx);
    subset_clu=Exp.osp.clu(idx);
    
    
    tv=linspace(t0,t1,100);
    tsp=mode(diff(tv));
    for cc = 1:NC
        spikeRate2(:,cc,trial) = histcounts(subset_st(subset_clu==spikeIds(cc)), tv)./tsp;
    end
    end
end



msr=mean(spikeRate2,3);
figure(1)
plot(msr(:,hero))
tv=linspace(-10,10,100);

%plot.raster(time, neuron_id)

%% Exp.slist has saccadetimes
% Exp slist(1) is saccadetimes, try saccade triggered av, then try
% histogram of lags

spikeIds=unique(Exp.osp.clu);
NC=length(spikeIds);
stimes=Exp.slist(:,1);
spikeRate2 = zeros(99, NC, length(stimes));
%
for ii = 1:10:length(stimes)
    tt=stimes(ii);
    %Convert start time (In PTB time) to ephys time
    t_=Exp.ptb2Ephys(tt); 
    
    t0=t_-.01;
    t1=t_+.4;
    idx=Exp.osp.st>t0&Exp.osp.st<t1;
    subset_st=Exp.osp.st(idx);
    subset_clu=Exp.osp.clu(idx);
    
    
    tv=linspace(t0,t1,100);
    tsp=mode(diff(tv));
    for cc = 1:NC
        spikeRate2(:,cc,ii) = histcounts(subset_st(subset_clu==spikeIds(cc)), tv)./tsp;
    end

end


%%
msr=mean(spikeRate2,3);
figure(1);clf
%
tv=linspace(-.01,.4,100);


subplot(211)
plot(tv(2:end),msr)
%
subplot(212)
plot(tv(2:end),mean(msr'))

%% Only during static trials
staticfixtrials=find(cellfun(@(x) strcmp(x.PR.name,'FixedProceduralNoise'),Exp.D));
trial_start_times=cellfun(@(x) x.STARTCLOCKTIME,Exp.D);
trial_end_times=cellfun(@(x) x.ENDCLOCKTIME,Exp.D);

staticwindows=[trial_start_times(staticfixtrials) trial_end_times(staticfixtrials)];
for sac=1:length(stimes)
insac(sac)=sum(stimes(sac)>staticwindows(:,1) & stimes(sac)<staticwindows(:,2));
end
%%
stimes3=stimes(insac>0);
spikeRate3 = zeros(99, NC, length(stimes3));
%%
for ii = 1:length(stimes3)
    tt=stimes3(ii);
    %Convert start time (In PTB time) to ephys time
    t_=Exp.ptb2Ephys(tt); 
    
    t0=t_-.1;
    t1=t_+.4;
    idx=Exp.osp.st>t0&Exp.osp.st<t1;
    subset_st=Exp.osp.st(idx);
    subset_clu=Exp.osp.clu(idx);
    
    
    tv=linspace(t0,t1,100);
    tsp=mode(diff(tv));
    for cc = 1:NC
        spikeRate3(:,cc,ii) = histcounts(subset_st(subset_clu==spikeIds(cc)), tv)./tsp;
    end

end

%
%%
msr3=mean(spikeRate3,3);
figure(2);clf
subplot(211)
tv=linspace(-.1,.4,100);
%
subplot(211)
plot(tv(2:end),msr3)
%
subplot(212)
plot(tv(2:end),mean(msr3'))


%% Definitely two peaks

goodre=mean(mean(spikeRate3,1),3)>.1;
msr3re=msr3(:,goodre,:);
subplot(211)
plot(tv(2:end),msr3re)
%
subplot(212)
plot(tv(2:end),mean(msr3re'))


%%
diffcheck=diff(stimes3);
diffcheck(diffcheck>1)=nan;
figure(4)
plot(sort(diffcheck))
%%
figure(4)
hist(diffcheck,0:0.01:1)
ylabel('count')
xlabel('time bins (s)')
title('time between saccades during static stim')
%%
% Exp.osp.st=Exp.osp.st-0.5;
% Exp.osp.st(Exp.osp.st<0)=0; %don't do this, but for now 
%%
eyePos=Exp.vpx.smo(:,2:3);


%%


