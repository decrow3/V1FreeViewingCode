%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory
cd('/home/huklab/Documents/NPX_pilot')
 %% Get timing strobes from spikeGLX setup
% timingfilename='/mnt/NPX/Rocky/20230608/Rocky-v1-0608_g0/Rocky-v1-0608_g0_t0.nidq.bin';
% nchannels=4;
% samples=[];
% dataformat=[];
% memorymap=1;
% NiDaQout = spGLX_load_binary_file(timingfilename,nchannels,samples,dataformat,memorymap);
% ttl=dec2bin(NiDaQout.Data.data(4,:));
% 
% 
% datapixx=ttl(:,1:6)=='1';
% 
% 
% Event.samplenumber=find(sum(datapixx,2)>0);
% Event.ts=(Event.samplenumber)/30003.000300;% from nidaq metafile
% Event.signal=ttl(Event.samplenumber,1:6);
% %
% imecSync=ttl(:,8)=='1';
% %% Received signals from imec card for syncing 
% spikefilename='/mnt/NPX/Rocky/20230608/Rocky-v1-0608_g0/Rocky-v1-0608_g0_imec0/Rocky-v1-0608_g0_t0.imec0.ap.bin';
% nchannels=385;
% samples=[];
% dataformat=[];
% memorymap=1;
% Tout = spGLX_load_binary_file(spikefilename,nchannels,samples,dataformat,memorymap);
% imecSent=Tout.Data.data(385,:); %this takes forever
% %Signals that were sent from imec card
% riseRec=find(diff(imecSync)>0); %in samples on nidaq
% riseSent=find(diff(imecSent)>40); %in samples on AP
% %%
% %save('/home/huklab/Documents/SpikeSorting/Output/2023-06-08/Rocky-0608_g0_t0.nidq_imecAp.mat','riseSent','riseRec','imecSent','imecSync')
% load('/home/huklab/Documents/SpikeSorting/Output/2023-06-08/Rocky-0608_g0_t0.nidq_imecAp.mat')
% %%
% 
% Event.ClosestAPSamples = interp1(riseRec,riseSent,Event.samplenumber);
% Event.Ap_ts=Event.ClosestAPSamples/30000.076569037657; % from AP metafile
% 
% risephD=find(diff(phDiode)>0); %in samples on nidaq
% Event.PhDiodeClosestAPSamples = interp1(riseRec,riseSent,risephD);
% Event.PD_ts=Event.PhDiodeClosestAPSamples/30000.076569037657; % from AP metafile

%% Events recorded on NiDaq are now in terms of seconds from start of NPX
%recording

% save(['/home/huklab/Documents/SpikeSorting/Output/2023-06-08/' 'Rocky-0608_g0_t0.nidq' 'dp_ev_ts2.mat'],'Event')
% 
% close all
% clear Tout NiDaQout
%%
%/mnt/NPX/Rocky/20230608/Rocky-v1-0608_g0/Rocky-v1-0608_g0_t0.nidq.bin'
load(['/home/huklab/Documents/SpikeSorting/Output/2023-06-08/' 'Rocky-0608_g0_t0.nidq' 'dp_ev_ts.mat'],'Event')

%%
RecStrobes=(Event.signal(:,6));
nStrobes=sum(Event.signal(:,6))
t=Event.Ap_ts(RecStrobes);

Recbits=Event.signal(:,1:5);
RecWords=Recbits(sum(Recbits,2)>0,:);
RecWords=bin2dec(num2str(RecWords));
RecWords_ts=Event.Ap_ts(sum(Recbits,2)>0);

scatter(RecWords_ts,RecWords,'.')

%% Load stimuli data from once data was recording, move calibration/setup etc to subfolder
StimPath='/mnt/NPX/Rocky/20230608/Stim/';
Exp = io.basic_marmoview_import(StimPath);

%% Check number of timestamps
ptbsw = cell2mat(cellfun(@(x) x.STARTCLOCK(:)', Exp.D, 'uni', 0));
ptbew = cell2mat(cellfun(@(x) x.ENDCLOCK(:)', Exp.D, 'uni', 0));

ptbwords = reshape([ptbsw'; ptbew'], [], 1); % all strobe words

fprintf('PTB sent %d strobed words during this session\n', numel(ptbwords))
fprintf('Ephys recovered %d strobes (flip on channel 3)\n', nStrobes)

%% Prepare to synch using xcorr
ptbst = cellfun(@(x) x.STARTCLOCKTIME, Exp.D); % ptb trial starts
ptbet = cellfun(@(x) x.ENDCLOCKTIME, Exp.D); % ptb trial ends

%% Plot sent signals
figure(1);clf
hold on
scatter(ptbst,ptbsw(:,1))
scatter(ptbst,ptbsw(:,2))
scatter(ptbst,ptbsw(:,3))
scatter(ptbst,ptbsw(:,4))
scatter(ptbst,ptbsw(:,5))
scatter(ptbst,ptbsw(:,6))


%% Remove unrecorded times for syncing
%
ptbt = reshape([ptbst'; ptbet'], [], 1); % all strobe times from ptb

%%
figure(1); clf
plot(diff(ptbt), 'o')
hold on


%t = ts(find(diff(ts) > 2e-3) + 1);
%figure(2); clf
%T = max([t-t(1);ptbt-ptbt(1)]);
T=max([ptbt-ptbt(1)]);

% bin timestamps to find best timelag
bin_size = 1e-3;
bins = 0:bin_size:T;
ecnt = histcounts(t-t(1), bins);
pcnt = histcounts(ptbt-ptbt(1), bins);

figure(3); clf; set(gcf, 'Color', 'w')
plot(bins(1:end-1), pcnt); hold on
plot(bins(1:end-1), ecnt, '-')


legend({'PTB Strobes', 'EPHYS Strobes'})
xlabel('Time from first Strobe (seconds)')
title('Binned Strobe Times')

%% Get best lag to match time series
[out, lags] = xcorr(pcnt, ecnt, 'Coef');
[~, id] = max(out);

figure(4); clf; set(gcf, 'Color', 'w')
plot(lags*bin_size, out); hold on
xlabel('Lag (seconds)')
ylabel('Cross Correlation')
plot(lags(id)*bin_size, out(id), 'o')
legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')

offset = lags(id)*bin_size;
fprintf('Best guess for initial offset: %02.2f (s)\n', offset)

%% re-bin with lag accounted to match times
bin_size = 1e-3;
bins = 0:bin_size:T;

[ecnt, ~, eid] = histcounts(t-t(1)+offset, bins);
[pcnt, ~, pid] = histcounts(ptbt-ptbt(1), bins);

figure(4); clf; set(gcf, 'Color', 'w')
plot(bins(1:end-1), pcnt); hold on
plot(bins(1:end-1), ecnt, '.')
%xlim([3000, 7000])
legend({'PTB Strobes', 'EPHYS Strobes'})
xlabel('Time from first PTB Strobe (seconds)')
title('Binned Strobe Times')


%%
 potential_matches = find(ecnt==1 & pcnt==1);
%potential_matches = find(ecnt>0 & pcnt>0);
fprintf('found %d potential matches at lag %02.3fs\n', numel(potential_matches), lags(id)*bin_size)

ematches = ismember(eid,potential_matches);
pmatches = ismember(pid, potential_matches);

assert(sum(ematches)==sum(pmatches), 'numer of matching timestamps does not match')
clf

%%

fun = synchtime.align_clocks(ptbt(pmatches), t(ematches));
plot(ptbt(pmatches), t(ematches), 'o'); hold on
plot(xlim, fun(xlim), 'k')
xlabel('PTB clock')
ylabel('Ephys clock')


sprintf('%15.15f',fun(0))
% '-1686257265.850085020065308'

%% Plot check
figure(2)
plot(ptbst+fun(0),ptbsw,'.')
hold on
scatter(RecWords_ts,RecWords,'.')
hold off
%%
x1=-1; x2=12200;
for ll=1:10
figure(2)
plot(ptbst+fun(0),rem(ptbsw,32),'.')
hold on
xlim([x1 x2])
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2])
pause(0.2)

end
%% Not looking great

%Exp.ptb2Ephys = fun;


%% Run Fixbiter.m
toffset=-1686257209.144262790679932;
sprintf('%15.15f',toffset)
w=[1 -toffset];
@(t)(t-w(2))/w(1)
myfun=@(t)(t-w(2))/w(1);

% %% Manual sync gave offset of -1686257209.140790939331055 for ptb2Ephys
% w=[1 1686257209.14079];
% @(t)(t-w(2))/w(1)
% myfun=@(t)(t-w(2))/w(1);

%%
plot(ptbst+myfun(0),rem(ptbsw,32),'.')
scatter(RecWords_ts,RecWords,'.')
%% check
x1=4100; x2=4200;
for ll=1:10
figure(2)
plot(ptbst+myfun(0),rem(ptbsw,32),'.')
hold off
xlim([x1 x2]); ylim([0 30]);
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2]); ylim([0 30]);
pause(0.2)

end

%%
Exp.ptb2Ephys = myfun;


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
% dataPath1 = [StimPath];
% EyeExp = io.basic_marmoview_import(dataPath1);
% fid = 1;
% [EyeExp, fig] = io.import_eye_position(EyeExp, dataPath1, 'fid', fid, 'zero_mean', true, 'normalize', true);


%%

% eyex=EyeExp.D{2, 1}.eyeData(2);
% eyey=EyeExp.D{2, 1}.eyeData(3);
% 
% cx=EyeExp.D{1, 1}.C.c(1);
% cy=EyeExp.D{1, 1}.C.c(2);
% dx=EyeExp.D{1, 1}.C.dx;
% dy=EyeExp.D{1, 1}.C.dy;
% 
%     % convert to d.v.a.
%     eyex2 = (eyex- cx)/(dx * Exp.S.pixPerDeg);
%     eyey2 = (eyey - cy)/(dy * Exp.S.pixPerDeg);

%%
% C = calibGUI(EyeExp,'use_smo','true');
% C.cmat
%%
%0.0061   -0.0035    0.0008    0.0004   -0.0002
% cmat=[1.1834    1.2302    0.0004    2.4289   -1.1900];
% cmat=[1.1775    1.2302    0.0004    2.4289   -1.1900];

cmat=[1.1775    1.2302    0.0004    2.4249   -1.0550];
%C.cmat=cmat;
%%
% eyePos = EyeExp.vpx.smo(:,2:3);
% eyePos =get_eyepos(C);
% fig = io.checkCalibration(EyeExp, eyePos);
% %

%% And data set
fid=1;
[Exp, fig] = io.import_eye_position(Exp, StimPath, 'fid', fid, 'zero_mean', true, 'normalize', true);
eyePosRaw = Exp.vpx.smo(:,2:3);
eyePosRaw = Exp.vpx.raw0(:,2:3);


            th = cmat(3);
            R = [cosd(th) -sind(th); sind(th) cosd(th)];
            S = [cmat(1) 0; 0 cmat(2)];
            A = (R*S)';
            Ainv = pinv(A);

            eyePos = (Exp.vpx.raw0(:,2:3) - cmat(4:5))*Ainv;
            

%Why doesn't this work???
fig = io.checkCalibration(Exp, eyePos);
% Now it works?
%%
C = calibGUI(Exp,'use_smo','true');
%%
C.cmat=[0.6344    0.6812   -3.0005  -42.1253   19.3845];
%  0.6645   -0.6898   -0.0005  -42.1243   19.4055
% 
% % %Put data into calibration object
% C.Exp=Exp;
%%
eyePos =C.get_eyepos();
fig = io.checkCalibration(Exp, eyePos);

%% Update the eyePos in Exp (Overwriting the raw data feels like the wrong way to do this, but subsequent analysis expect it?)
Exp.vpx.smo(:,2:3)=eyePos;

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


%% Load spikes
dataPath = '/home/huklab/Documents/SpikeSorting/Output/2023-06-08/';
Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp.st = Exp.osp.st - .5;
Exp.osp.st(Exp.osp.st<0)=0;
% move zero to the end
id = max(Exp.osp.clu)+1;

Exp.osp.clu(Exp.osp.clu==0) = id;
Exp.osp.cids = Exp.osp.cids2;
Exp.osp.cids(Exp.osp.cids==0) = id;

lbledgood=find(Exp.osp.cgs==2)+1;

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
if ~isempty(dotTrials)
    
    BIGROI = [-10 -10 10 10];
    %eyePos = C.refine_calibration();
    binSize = .2;
    Frate = 60;
    [Xstim, RobsSpace, opts] = io.preprocess_spatialmapping_data(Exp, ...
        'ROI', BIGROI*Exp.S.pixPerDeg, 'binSize', binSize*Exp.S.pixPerDeg, ...
        'eyePosExclusion', 2e3, ...
        'eyePos', eyePos, 'frate', Frate, 'fastBinning', true);
    
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
legend(h, {'All units', 'Good Units', '1 Spike/Sec'}, 'Location', 'Best')

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
figure(1); clf
for ilag = 1:nlags
    subplot(2, ceil(nlags/2), ilag)
    imagesc(opts.xax, opts.yax, reshape(stas(ilag, :), opts.dims), wm)
end

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
    imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2])
    colormap(plot.coolwarm)
    title(Exp.osp.cids(goodunits(cc)))
end
    

%% plot subset
subset = [407 408 411 419 473 477 507 520 570 571 577 580 581 585 609 611 ...
    656 657 677 695 700 708 709 714 722 725 735 748 775 799 842 844 850 889 890 891 902 908 911  ...
    933 934 944 945 946 950 955 968 975 1002 1010 1021 1022 1023 1029 1031 1041 1048 1077 1100 1152 1181 1184 1189 1192 ...
    1202 1219 1239 1252 1253 1274 1277 1284 1286 1292 1293 1307 1308 1314 1318 1350 1364 1399 1402 1414 1419 1451 1462 ...
    1467 1485 1487 1500 1505 1506 1513 1521 1655 1661 1668 1671 1672 1725 1736 1958 1959 1970 1982];
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
si=cc;
    subplot(sy,sx, si)

    I = squeeze(stas(:,:,ccid(cc)));
    I = (I - mean(I(:)) )/ std(I(:));

    [bestlag,j] = find(I==max(I(:)));
    I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
    I = mean(I);
%     I = std(I);
    xax = opts.xax/Exp.S.pixPerDeg;
    yax = opts.yax/Exp.S.pixPerDeg;
    imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2])
    colormap(plot.coolwarm)
    title(Exp.osp.cids(goodunits(ccid(cc))))
    title(Exp.osp.clusterDepths(goodunits(ccid(cc))))
    grid on
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
    
    t0=t_-.1;
    t1=t_+.8;
    idx=Exp.osp.st>t0&Exp.osp.st<t1;
    subset_st=Exp.osp.st(idx);
    subset_clu=Exp.osp.clu(idx);
    
    
    tv=linspace(t0,t1,100);
    tsp=mode(diff(tv));
    for cc = 1:NC
        spikeRate2(:,cc,ii) = histcounts(subset_st(subset_clu==spikeIds(cc)), tv)./tsp;
    end

end



msr=mean(spikeRate2,3);
figure(1);clf
%%
tv=linspace(-.1,.8,100);
plot(tv(2:end),msr)


%%
Exp.osp.st=Exp.osp.st-0.5;
Exp.osp.st(Exp.osp.st<0)=0; %don't do this, but for now 
%%
eyePos=Exp.vpx.smo(:,2:3);
