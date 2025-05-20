%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory
cd('/home/huklab/Documents/NPX_pilot')

%% 
ephysPath='/mnt/NPX/Rocky/20240427/Rocky20240427_V1V2_g0';
%These were unfortunately mislabeled
ephysPrefix='Rocky20240427_V1V2_g0';

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



 %timingfilename='/mnt/NPX/Rocky/20240427/Rocky20240427_V2_g0/Rocky20240427_V2_g0_t0.nidq.bin';
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
Event.Ni_ts=(Event.samplenumber)/30003.0003;% from nidaq metafile
Event.signal=ttl(Event.samplenumber,1:7);

%Final signal is recorded on first bit, least significant digit so last in
%byte order
imecSync=ttl(:,8)=='1';
riseRec=find(diff(imecSync)>0); %Recorded rising edges in samples on nidaq

%% This is typically a 1Hz signal from the imec card that signals the start
%of recording on the first rising edge and maintains a pulse to correct for
%drifting of clocks between the imec and nidaq cards

%Analog signals are also recorded on the nidaq:
% PDanalog=NiDaQout.Data.data(1,:);
PDanalog=NiDaQout.Data.data(3,:);

PD=PDanalog>3*10^3;
risePD=find(diff(PD)>0);% Rising edges
risePD_Ni_ts=(risePD)/30003.0003;
fallPD=find(diff(PD)<0);%Falling edges
fallPD_Ni_ts=(fallPD)/30003.0003;% 



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
if isempty(dir('/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/'))
    mkdir('/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/')
end
 save('/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/Rocky20240427_V2_g0_t0.nidq_imecAp.mat','riseSent','riseRec','imecSent','imecSync','Event')
% clear Tout
 %load('/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/Rocky20240427_V2_g0_t0.nidq_imecAp.mat')

%% Converting signals from nidaq to ephys(imec) time
% Datapixx 'events'
Event.ClosestAPSamples = interp1(riseRec,riseSent,Event.samplenumber);
Event.Ap_ts=Event.ClosestAPSamples/30000.076569037657; % from AP metafile

%Photodiode rise times
%risephD=find(diff(phDiode)>0); %in samples on nidaq
Event.risePhDiodeClosestAPSamples = interp1(riseRec,riseSent,risePD);
Event.risePD_Ap_ts=Event.risePhDiodeClosestAPSamples/30000.076569037657; % from AP metafile
Event.fallPhDiodeClosestAPSamples = interp1(riseRec,riseSent,fallPD);
Event.fallPD_Ap_ts=Event.fallPhDiodeClosestAPSamples/30000.076569037657; % from AP metafile

% Saving everything into an Events structure for saving / ease of loading


%% Events recorded on NiDaq are now in terms of seconds on Ephys(Imec) clock from t=0s at start of NPX recording

save(['/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/' 'Rocky20240427_V2_t0.nidq' 'dp_ev_ts2.mat'],'Event')
%load(['/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/' 'Rocky20240427_V2_t0.nidq' 'dp_ev_ts2.mat'],'Event')


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



%% Load stimuli data from once data was recording, move calibration/setup etc to subfolder
StimPath='/mnt/NPX/Rocky/20240427/Rocky20240427/';
Exp = io.basic_marmoview_import(StimPath);


%% Remove broken trails (failed before start of first stimulus)
keep=cellfun(@(x) isfield(x,'STARTCLOCKTIME'),Exp.D);
Exp.D=Exp.D(keep);

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
%xlim([6280 6680]); ylim([0 64]);
subplot(212)
plot(ptbt_,ptbwords,'.')
ptbtime0=min(ptbt_);
ptbtime1=max(ptbt_);
%Cropping can be usefull
% xlim([ptbtime0 ptbtime0]); ylim([0 64]);


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
%xlim([ptbtime0 ptbtime1]); ylim([0 64]);

%% if nec, update
Recbits=switchbits;
RecWords=bin2dec(Recbits);

%% 
strangehighs=(RecWords==62|RecWords==63);
RecWords=RecWords(~strangehighs);
RecWords_ts=RecWords_ts(~strangehighs);

%% Rough check of range, crop time to after start of recording

ptbtime0=min(ptbt_)+11500;
ptbtime1=max(ptbt_);


subplot(211)
scatter(RecWords_ts,RecWords,'.')
ylim([0 64]);
subplot(212)
plot(ptbt_,ptbwords,'.')
ylim([0 64]);

xlim([ptbtime0 ptbtime1]); % ylim([0 64]);

%%

%% Remove unrecorded times for syncing

ptbt=ptbt_(ptbt_>ptbtime0);
ptbwords=ptbwords(ptbt_>ptbtime0);

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
plot(ptbst+myfun(0),(ptbsw),'.')
hold on
pause(0.1)
scatter(RecWords_ts,RecWords,'.')
hold off

%% check
clf
x1=lowbound; x2=highbound;
for ll=1:10
figure(2);subplot(111)
plot(ptbst+myfun(0),(ptbsw),'.')
hold on
xlim([x1 x2]); ylim([0 64]);
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2]); ylim([0 64]);
pause(0.2)

end



%% Correcting for drift, best to take matches

wordstocheck=unique(ptbwords);
clear tdiff rec_match ptb_match
for ii=1:length(wordstocheck)
    word=wordstocheck(ii);
    ptb_ind=find(ptbwords==word);
    rec_ind=find(RecWords==word);
    if nnz(ptb_ind)>nnz(rec_ind) % more ptb, loop through rec
        for jj=1: nnz(rec_ind)
            [tdiff{ii}(jj),ptb_match{ii}(jj)]=min(abs(ptbt+myfun(0)-RecWords_ts(rec_ind(jj))));
            rec_match{ii}(jj)=rec_ind(jj);
        end
    else % more rec, loop through ptb
        for jj=1: nnz(ptb_ind)
            [tdiff{ii}(jj),rec_match{ii}(jj)]=min(abs(ptbt(ptb_ind(jj))+myfun(0)-RecWords_ts));
            ptb_match{ii}(jj)=ptb_ind(jj);
        end
    end
    max(tdiff{ii})

end
pmatches=cell2mat(ptb_match);
ematches=cell2mat(rec_match);

fun = synchtime.align_clocks(ptbt(pmatches), RecWords_ts(ematches));
plot(ptbt(pmatches), RecWords_ts(ematches), 'o'); hold on
plot(xlim, fun(xlim), 'k')
xlabel('PTB clock')
ylabel('Ephys clock')
%%
clf
x1=lowbound; x2=highbound;
for ll=1:10
figure(2);subplot(111)
plot(fun(ptbst),(ptbsw),'.')
hold on
xlim([x1 x2]); ylim([0 64]);
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2]); ylim([0 64]);
pause(0.2)

end

%% If happy lock it in!!!
Exp.ptb2Ephys = fun;




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
xlim([1000 5500]);ylim([0.8 1.2])

%~20ms

%% Analog PD signal on ephys time
nidaqtime_s=1:length(PDanalog)/30003.0003;% from nidaq metafile
nidaqClosestAPSamples = interp1(riseRec,riseSent,1:length(PDanalog));
nidaqsampon_Ap_ts=nidaqClosestAPSamples/30000.076569037657; % from AP metafile
%%
PD_PTBinAPtimes(PD_PTBinAPtimes<0)=[];
%NOTE; Flashtimes were saved with the 'currenttime' before the draw of the screen with
%the photodiode, need to add in an additional screenflip (1/60)

%%
figure(11);clf
hold off
% plot(nidaqsampon_Ap_ts(8e7:4:end-2e8),PDanalog(8e7:4:end-2e8))
plot(nidaqsampon_Ap_ts(1:4:end),PDanalog(1:4:end))
hold on 
plot(PD_PTBinAPtimes+(1/60),4000*ones(1,length(PD_PTBinAPtimes)),'.')
%xlim([3000 9000]);
%xlim([5700 5860]);
% 20-30ms? 25ms of system delay to add in to correct screen throws
% add screendelay for each stim, or just subtract 0.025 from spike times

legend('measuredPDtrace','PTBfliptimes')





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
% dataPath1 ='/mnt/NPX/Rocky/20240427/Rocky 2024-04-27/facecal';
% EyeExp = io.basic_marmoview_import(dataPath1);
% fid = 2;
%  [EyeExp, fig] = io.import_eye_position(EyeExp, dataPath1, 'fid', fid, 'zero_mean', false, 'normalize', false);
% %% 
% C0 = calibGUI(EyeExp,'use_smo','true');
% C0.cmat

%%
fid = 1;
[Exp, fig] = io.import_eye_position(Exp, StimPath, 'fid', fid, 'zero_mean', false, 'normalize', false);
%[Exp, fig] = io.import_eye_position_ddpi12092023(Exp, dataPath1, 'fid', fid, 'zero_mean', true, 'normalize', true);

%%
caltrials=find(cellfun(@(x) strcmp(x.PR.name,'FaceCal'),Exp.D));
EyeExp=Exp;
EyeExp.D=Exp.D(caltrials(1:end-100));
%%
C = calibGUI(EyeExp,'use_smo','true');
C.cmat= [   0.9169    0.9533   -3.9989    0.1143   0.2144];

%%
clear C
C = calibGUI(Exp,'use_smo','true');
C.cmat= [   0.9169    0.9533   -3.9989    0.1143   0.2144];
%1.1687    1.2130   -4.9988   -0.6326    0.0306
%% Saving Cmat so I don't need to refine again

% cmat=[1.2028    1.0752   -2.9994   -2.5274   -3.0364]; % online track
%cmat=[0.1429    0.1401   -1.9981  -39.9112   20.2618]; %Offline dpi
%cmat=[1.6017    1.1430   -5.6833   -0.4870    0.0254]
%cmat=[ 11.3761    1.0767   -6.6919   -0.4219    0.0254];
% cmat= [1.2477    1.1346   -1.9999   -4.3242   -0.7232];
% C.cmat=cmat;
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
clear C
C = calibGUI(EyeExp,'use_smo',0);
%% Raw data is way different
%  cmat=[11.2678   10.7247   -3.3977  326.9572 -233.3391];%online t
%cmat=[7.0738    7.4387   -0.8886  -91.1588  -60.6488];%offline dpi
cmat=[ 5.8750    6.3559   -1.1063  -32.3766  -65.0466];
C.cmat=cmat
%%
clear C
C = calibGUI(Exp,'use_smo',0);
C.cmat=[ 5.8750    6.3559   -1.1063  -32.3766  -65.0466];
C.cmat=cmat; %update plot
%% then get full session Eyetrack
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
Exp.vpx.Labels(isnan(Exp.vpx.raw(:,2))) = 4;

validTrials = io.getValidTrials(Exp, 'Grating');
for iTrial = validTrials(:)'
    if ~isfield(Exp.D{iTrial}.PR, 'frozenSequence')
        Exp.D{iTrial}.PR.frozenSequence = false;
    end
end
        
fprintf(fid, 'Done importing session\n');
%% save here if you aren't ready to put in spikes
Exp.FileTag = 'Rocky20240427.mat';
save(fullfile('/home/huklab/Documents/NPX_pilot/Output', Exp.FileTag), '-v7.3', '-struct', 'Exp')
%Exp=load(fullfile('/home/huklab/Documents/NPX_pilot/Output', Exp.FileTag));
%% Load spikes
dataPath = '/home/huklab/Documents/SpikeSorting/Output/Rocky-2024-04-27-imec0-redo';
Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp.st = Exp.osp.st;
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

%% Second probe
dataPath = '/home/huklab/Documents/SpikeSorting/Output/Rocky/20240427-imec1';
Exp.osp2 = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp2.st = Exp.osp2.st;
Exp.osp2.st(Exp.osp2.st<0)=0;
% move zero to the end
id = max(Exp.osp2.clu)+1;

Exp.osp2.clu(Exp.osp2.clu==0) = id;
Exp.osp2.cids = Exp.osp2.cids2;
Exp.osp2.cids(Exp.osp2.cids==0) = id;

% lbledgood=find(Exp.osp2.cgs==2)+1;

% Exp.osp2.st = st_s.spikeTimes(:,1)+.6;
% Exp.osp2.clu = st_s.spikeTimes(:,2);
% Exp.osp2.depths = st_s.depths;
% Exp.osp2.templates = st_s.templates;
% Exp.osp2.cids = unique(Exp.osp2.clu);


figure(1); clf
plot.raster(Exp.osp2.st, Exp.osp2.clu, 1)
ylabel('Unit ID')
xlabel('Time (seconds)')

%%
NC1=max(Exp.osp.clu);
Exp.osp2.clu=Exp.osp2.clu+NC1;
Exp.osp2.cids=Exp.osp2.cids+NC1;

%% Combine
Exp.ospfull.probe=[ones(size(Exp.osp.st,1),1,'int8'); 2*ones(size(Exp.osp2.st,1),1,'int8')];
Exp.ospfull.st=[Exp.osp.st; Exp.osp2.st];
Exp.ospfull.clu=[Exp.osp.clu; Exp.osp2.clu+NC1];
Exp.ospfull.cids=[Exp.osp.cids; Exp.osp2.cids+NC1];
Exp.ospfull.clusterDepths=[Exp.osp.clusterDepths Exp.osp2.clusterDepths];

%%
Exp.osp1=Exp.osp;
Exp.osp=Exp.ospfull;

%% 
% %% Retry with handful
% %% Load spikes
% dataPath = '/home/huklab/Documents/SpikeSorting/Output/2024-04-27/Retry';
% Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
% Exp.osp.st = Exp.osp.st;
% Exp.osp.st(Exp.osp.st<0)=0;
% % move zero to the end
% id = max(Exp.osp.clu)+1;
% 
% Exp.osp.clu(Exp.osp.clu==0) = id;
% Exp.osp.cids = Exp.osp.cids2;
% Exp.osp.cids(Exp.osp.cids==0) = id;
% 
% % lbledgood=find(Exp.osp.cgs==2)+1;
% 
% % Exp.osp.st = st_s.spikeTimes(:,1)+.6;
% % Exp.osp.clu = st_s.spikeTimes(:,2);
% % Exp.osp.depths = st_s.depths;
% % Exp.osp.templates = st_s.templates;
% % Exp.osp.cids = unique(Exp.osp.clu);
% 
% 
% figure(1); clf
% plot.raster(Exp.osp.st, Exp.osp.clu, 1)
% ylabel('Unit ID')
% xlabel('Time (seconds)')
% 
% 
% %% Retry with full re-sort
% 
% Exp=load('/home/huklab/Documents/NPX_pilot/Rock20240427_V1V2_onlinedpi.mat');
% resync=load('/home/huklab/Documents/NPX_pilot/Rocky20240427_resync.mat');
% Exp.ptb2Ephys=resync.fun;
% Exp.FileTag='Rocky20240427_V1V2_resorted';
% %% Load spikes
% dataPath = '/home/huklab/Documents/SpikeSorting/Output/2024-04-27_FullyResorted/Rocky_2024-04-27-prb1/';
% Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
% Exp.osp.st = Exp.osp.st;
% Exp.osp.st(Exp.osp.st<0)=0;
% % move zero to the end
% id = max(Exp.osp.clu)+1;
% 
% Exp.osp.clu(Exp.osp.clu==0) = id;
% Exp.osp.cids = Exp.osp.cids2;
% Exp.osp.cids(Exp.osp.cids==0) = id;
% 
% % lbledgood=find(Exp.osp.cgs==2)+1;
% 
% % Exp.osp.st = st_s.spikeTimes(:,1)+.6;
% % Exp.osp.clu = st_s.spikeTimes(:,2);
% % Exp.osp.depths = st_s.depths;
% % Exp.osp.templates = st_s.templates;
% % Exp.osp.cids = unique(Exp.osp.clu);
% 
% 
% figure(1); clf
% plot.raster(Exp.osp.st, Exp.osp.clu, 1)
% ylabel('Unit ID')
% xlabel('Time (seconds)')
% 
% %% Second probe
% dataPath = '/home/huklab/Documents/SpikeSorting/Output/2024-04-27_FullyResorted/Rocky_2024-04-27-prb2/';
% Exp.osp2 = load(fullfile(dataPath, 'spkilo.mat'));
% Exp.osp2.st = Exp.osp2.st;
% Exp.osp2.st(Exp.osp2.st<0)=0;
% % move zero to the end
% id = max(Exp.osp2.clu)+1;
% 
% Exp.osp2.clu(Exp.osp2.clu==0) = id;
% Exp.osp2.cids = Exp.osp2.cids2;
% Exp.osp2.cids(Exp.osp2.cids==0) = id;
% 
% % lbledgood=find(Exp.osp2.cgs==2)+1;
% 
% % Exp.osp2.st = st_s.spikeTimes(:,1)+.6;
% % Exp.osp2.clu = st_s.spikeTimes(:,2);
% % Exp.osp2.depths = st_s.depths;
% % Exp.osp2.templates = st_s.templates;
% % Exp.osp2.cids = unique(Exp.osp2.clu);
% 
% 
% figure(1); clf
% plot.raster(Exp.osp2.st, Exp.osp2.clu, 1)
% ylabel('Unit ID')
% xlabel('Time (seconds)')
% 
% %%
% NC1=max(Exp.osp.clu);
% Exp.osp2.clu=Exp.osp2.clu+NC1;
% Exp.osp2.cids=Exp.osp2.cids+NC1;
% 
% %% Combine
% Exp.ospfull.probe=[ones(size(Exp.osp.st,1),1,'int8'); 2*ones(size(Exp.osp2.st,1),1,'int8')];
% Exp.ospfull.st=[Exp.osp.st; Exp.osp2.st];
% Exp.ospfull.clu=[Exp.osp.clu; Exp.osp2.clu+NC1];
% Exp.ospfull.cids=[Exp.osp.cids; Exp.osp2.cids+NC1];
% Exp.ospfull.clusterDepths=[Exp.osp.clusterDepths Exp.osp2.clusterDepths];
% 



%% Get drifting gratings
%validTrials = io.getValidTrials(Exp, 'DriftingGrating');



%%
figure(1); clf
plot.raster(Exp.osp.st, Exp.osp.clu, 1); hold on
ylabel('Unit ID')
xlabel('Time (seconds)')
cmap = lines;
stims = {'Dots', 'DriftingGrating', 'BackImage', 'FaceCal','FixRsvpStim'};

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
%dotTrials=dotTrials(dotTrials<1518);
Exp_subset=Exp;
%Exp_subset.D={Exp.D{dotTrials}}';
if ~isempty(dotTrials)
    
    BIGROI = [-15 -15 15 15];
    %eyePos = C.refine_calibration();
    binSize = .5;
    Frate = 60;
    [Xstim, RobsSpace, opts] = io.preprocess_spatialmapping_data(Exp, ...
        'ROI', BIGROI*Exp.S.pixPerDeg, 'binSize', binSize*Exp.S.pixPerDeg, ...
        'eyePosExclusion', 2e3, ...
        'eyePos', Exp.vpx.smo(:,2:3), 'frate', Frate, 'fastBinning', true);
    
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
%%
figure(1); clf
for ilag = 1:nlags
    subplot(2, ceil(nlags/2), ilag)
    imagesc(opts.xax, opts.yax, reshape(stas(ilag, :), opts.dims), wm)
end
%%
figure(1); clf
for ilag = 6
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
%R = RobsSpace;
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
        set(gcf,'Position',[0 0 1080 1080]) 
    end

    si = mod(cc, nperfig);  
    if si==0
        si = nperfig;
    end
    subplot(sy,sx, si)

    I = squeeze(stas(:,:,cc));
    I = (I - mean(I(:)) )/ std(I(:));

    [bestlag,j] = find(I==max(I(:)));
    if ~isempty(bestlag)
    I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
    I = mean(I);
%     I = std(I);
    xax = opts.xax/Exp.S.pixPerDeg;
    yax = opts.yax/Exp.S.pixPerDeg;
    imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2]);
    hold on
        scatter(0,0,'k+');
        hold off
    axis equal tight xy
%    set(gca, 'XTick',[],'YTick',[])
    colormap(plot.coolwarm)
%      title((goodunits(cc)))
    title(((cc)))
%      title(Exp.osp.cids(goodunits(cc)))
     title(5000-Exp.osp.clusterDepths(goodunits((cc))))
    end
end
    

%% Print figures to pdf
for fig=1:nfigs
    figure(fig)
    set(gcf,'Position',[0 0 1080 1080])
saveas(figure(fig),['/home/huklab/Documents/NPX_pilot/Output/RockyV2V1_20240427/Rocky20240427_V2V1_dotRFs/Rocky20240427_V2V1_dotRFs_' num2str(nperfig*(fig-1)+1) '_' num2str(nperfig*fig) '.pdf'])
end


%% plot subset
% subset = [1 5 8 12 15 17 19 22 26 60 61 70 71 77 80 83 84 99 107 111 114 120 123 ...
%     132 140 170 178 183 193 209 218 219 220 221 226 228 232 233 234 238 239 240 241 242 243 244 246 ...
%     247 255 256 257 258 268 269 270 272 276 277 280 281 282 290 297 298 301 311 314 323 324 325 ...
%     333 339 341 345 348 349 350 351 354 358 359 360 362 364 366 375 377 384 410 411 419 430 ...
%     447 450 469 471 473 477 491 497 502 506 514 518 570 571 591 595 598 600 606 647 648 650 653 658 661 677 703 705 706 788 812 ...
%     833 834 845 847 866 872 878 881 913 955 977 982 1006 1009 1033 1037 1059 1068 1081 1082 1083 1148 1150 1154 1155 ...
%     1167 1179 1180 1181 1183 1184 1186 1188 1191];
% subset = [1 2 4 5 8  9 10 12 13 15 16 17 18 19 20 21 22 23 24 25 ...
%     26 27 28 29 30 31 32 33 34 35 41 43 44 45 46 47 48 49 50 ...
%     51 52 53 54 55 56 57 61 62 63 64 69 74 ...
%     80 81 82 83 84 85 87 89 90 93 95 97 100 ...
%     102 103 106 108 111 117 122 123 ...
%     129 132 133 134 137 138 140 147 ...
%     152 161 162 163 164 165 168 173 175 ...
%     176 177 178 180 181 187 188 191 192 195 ...    
% ];
subset=1:size(Exp.osp.cids);
nsub = numel(subset);
subset=Exp.osp.cids(subset);
%[~, ccid] = ismember(subset, Exp.osp.cids(goodunits));
[~, ccid] = ismember(subset, Exp.osp.cids);
[~, ccid] = Exp.osp.cids(goodunits);

nfigs = ceil(nsub / nperfig);
fig = ceil(1); clf
figure(fig); set(gcf, 'Color', 'w')

sx = ceil(sqrt(nsub));
sy = round(sqrt(nsub));
%[~,di]=sort(Exp.osp.clusterDepths(goodunits(ccid)));
[~,di]=sort(Exp.osp.clusterDepths((ccid)));
for ci = 1:nsub
    cc=di(ci);
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
    si=nsub-ci+1;
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
%        title(Exp.osp.cids(goodunits(ccid(cc))))

        %NOTE: depth is -ve (distance up from tip)
%          title(7000-560-2000-Exp.osp.clusterDepths(goodunits(ccid(cc))))
         title(7000-560-2000-Exp.osp.clusterDepths((ccid(cc))))
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

inds = io.getValidTrials(Exp, 'BackImage');
spikeIds=unique(Exp.osp.clu);
NC=length(spikeIds);
spikeRate = zeros(111, NC, length(inds));
most=mode(Exp.osp.clu);
hero=find(spikeIds==most);

figure(1); clf

%%
for trial = 1:length(inds)
    thisTrial=inds(trial);
   % if length(Exp.D{thisTrial}.PR.NoiseHistory(:,4))>0
%     tt=Exp.D{thisTrial}.PR.NoiseHistory(:,1);
    %plot( tt-tt(1) ,Exp.D{thisTrial}.PR.NoiseHistory(:,4)) ; hold on

    
    %Convert start time (In PTB time) to ephys time
    %t_=Exp.ptb2Ephys(Exp.D{thisTrial}.PR.NoiseHistory(1,1)); 

    tt=Exp.D{thisTrial}.STARTCLOCKTIME;
    t_=Exp.ptb2Ephys(tt); 

    t0=t_-.10;
    t1=t_+.5;
    idx=Exp.osp.st>t0&Exp.osp.st<t1;
    subset_st=Exp.osp.st(idx);
    subset_clu=Exp.osp.clu(idx);
    
    
    tv=linspace(t0,t1,100);
    tsp=mode(diff(tv));
    for cc = 1:NC
        spikeRate2(:,cc,trial) = histcounts(subset_st(subset_clu==spikeIds(cc)), tv)./tsp;
    end
  %  end
end


%%
msr=mean(spikeRate2,3);
figure(1);clf

tv=linspace(-.10,.50,100);
plot(tv(1:99),msr(:,:)')
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

%% Only during images trials
%staticfixtrials=find(cellfun(@(x) strcmp(x.PR.name,'FixedProceduralNoise'),Exp.D));
staticfixtrials=find(cellfun(@(x) strcmp(x.PR.name,'BackImage'),Exp.D));
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
tv=linspace(-.1,.4,100)-0.01;
%
subplot(211)
plot(tv(2:end),msr3)
%
subplot(212)
plot(tv(2:end),mean(msr3'))

%%
[~,order]=(max(msr3))
[~,ordered]=sort(order)
imagesc(msr3(:,ordered)')

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


