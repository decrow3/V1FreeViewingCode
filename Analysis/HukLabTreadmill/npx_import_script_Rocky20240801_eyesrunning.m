%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory
cd('/home/huklab/Documents/NPX_pilot')



%% Load stimuli data from once data was recording, move calibration/setup etc to subfolder
StimPath='/home/huklab/Documents/Eyetracking/Rocky 2024-08-01/';%/mnt/NPX/Rocky/20240801/Rocky20240801/';
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

ptbtime0=min(ptbt_)+6000;
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
% dataPath1 ='/mnt/NPX/Rocky/20240801/Rocky 2024-08-01/facecal';
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
%EyeExp.D=Exp.D(caltrials(1:end-100));
%%
C = calibGUI(EyeExp,'use_smo','true');
%C.cmat= [  0.9570    0.8249   -3.9992   -4.7732 -21.3140];

%%
clear C
C = calibGUI(Exp,'use_smo','true');
C.cmat= [  0.9570    0.8249   -3.9992   -4.7732 -21.3140];
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
% clear C
% C = calibGUI(EyeExp,'use_smo',0);
% %% Raw data is way different
% %  cmat=[11.2678   10.7247   -3.3977  326.9572 -233.3391];%online t
% %cmat=[7.0738    7.4387   -0.8886  -91.1588  -60.6488];%offline dpi
% cmat=[  5.9015    6.4353   -4.0889  -33.0357  -70.1293];
% C.cmat=cmat
%%
clear C
C = calibGUI(Exp,'use_smo',0);
C.cmat=[ 5.9015    6.4353   -4.0889  -33.0357  -70.1293];
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


%% Grating tuning prep
% Exp.spikeTimes = Exp.osp.st;
% Exp.spikeIds = Exp.osp.clu;
%%
jj=0;
nTrials=numel(Exp.D);
for ii=1:nTrials
    if length(Exp.D{ii}.inputs)>1
    Exp.D{ii}.treadmill=Exp.D{ii}.inputs{2};
    removetrial(ii)=0;
    jj=jj+1;
    Exp.D2{jj}=Exp.D{ii};
    else
    removetrial(ii)=1;
    end
end
%%
Exp.D0=Exp.D;
Exp.D=Exp.D2';

%% BREAK FOR DPI AND TREADSPEED EXPORT, NOT RUNNING get_drifting_grating_output
Exp.vpx0=Exp.vpx;
%%
DPIloadup20240801
%%
Exp.vpx.raw1=Exp.vpx.raw;
Exp.vpx.raw=Exp.vpx.dpi;
%%

Exp.vpx.onlinesmo=Exp.vpx.smo;
Exp.vpx.onlineraw=Exp.vpx.raw;
Exp.vpx.dvadpi(:,1)=Exp.vpx.dpi(:,1);
Exp.vpx.dvadpi(:,2)=Exp.vpx.dpi(:,2)/8.5 +3.5;
Exp.vpx.dvadpi(:,3)=-Exp.vpx.dpi(:,3)/8.5 + 7.5;
%%

% clear C
% C = calibGUI(Exp,'use_smo',0);
% %C.cmat=[ 5.9015    6.4353   -4.0889  -33.0357  -70.1293];
% %C.cmat=cmat; %update plot
% %% then get full session Eyetrack
% eyePos =C.get_eyepos();
% io.checkCalibration(Exp, eyePos);
%% Update
fix_mitchelllab_exports2

%%
Exp.spikeTimes = [min(Exp.vpx.smo(:,1)) max(Exp.vpx.smo(:,1))];
Exp.spikeIds = [0 0];

% convert to D struct format
D = io.get_drifting_grating_output(Exp,'truncate',false);

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


