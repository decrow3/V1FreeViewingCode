%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory

%%
load('/home/huklab/Documents/SpikeSorting/Output/2022-08-09-Combined-full/rez.mat')
tt=max(rez.st3(:,1))
tbatch=tt/rez.ops.Nbatch
dshifttimes=(1:tbatch:tt);
dshift=rez.dshift;
save('/home/huklab/Documents/SpikeSorting/Output/2022-08-09-Combined-full/shifts.mat','dshift','dshifttimes')
subplot(8,1,[1:7])
plot.raster(rez.st3(:,1),rez.st3(:,2),2)
xlim([0 tt])
subplot(8,1,8)
plot(1:tbatch:tt,rez.dshift)
xlim([0 tt])


%% Try to synchronize using strobed values
HUKDATA = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');
% dataPath = fullfile(HUKDATASHARE, 'gru', 'g20211217');
% figDir = '~/Google Drive/HuklabTreadmill/imported_sessions_qa/g20211217/';
dataPath = fullfile(HUKDATA, 'BrieMTsession','20220809');
figDir = fullfile(HUKDATA, 'BrieMTsession','20220809','out');%'~

%%
% File combiner ran and sorted, files saved to:
dataPath= '/media/huklab/Data/NPX/Brie/2022-08-09/Combined/';
FPGAtimestamps = load(fullfile(dataPath, 'FPGAtimestamps.mat'));

% try real sync --> fails
strobeChannel = find(FPGAtimestamps.ts.eventVirCh==8);
assert(mod(numel(strobeChannel), 2)==0, 'strobes missing')

nStrobes = numel(strobeChannel)/2;
ix = find(diff(FPGAtimestamps.ts.fullWords)~=0)+1;

tstruct = struct('tstrobes', FPGAtimestamps.ts.timestamps_s(ix), ...
    'strobes', FPGAtimestamps.ts.fullWords(ix));

%ts = timestamps.ts.timestamps_s;
ts = double(FPGAtimestamps.ts.timestamps)/30000;
val = double(FPGAtimestamps.ts.eventVirCh-1);

ephys_info.eventId = double(sign(FPGAtimestamps.ts.eventVirChStates)==1);
[t,v] = read_ephys.convert_data_to_strobes(val, ts, ephys_info);



figure(1); clf; set(gcf, 'Color', 'w')
plot(t, v, 'o')
hold on
plot(FPGAtimestamps.ts.timestamps_s, FPGAtimestamps.ts.fullWords, 'o')
title('(Nonsense) Strobes')
xlabel('Time (ephys clock)')
ylabel('Strobe Value')        

strobes = read_ephys.fix_missing_strobes(v);

hold off
    

%Try using t and v instead
tstruct = struct('tstrobes', t, 'strobes', v);
plot(t, strobes, 'o')
title('FPGA Strobes')
xlabel('Time (ephys clock)')
ylabel('Strobe Value')

%% Load stimuli/behavior data (without sync) and align to tstruct
stimPath= '/media/huklab/Data/NPX/Brie/2022-08-09/20220809';
Exp = io.basic_marmoview_import(stimPath, 'Timestamps', tstruct);

% This creates Exp.ptb2Ephys
% Not sure if this is actually the way I want to align these


%% Check number of timestamps
ptbsw = cell2mat(cellfun(@(x) x.STARTCLOCK(:)', Exp.D, 'uni', 0));
ptbew = cell2mat(cellfun(@(x) x.ENDCLOCK(:)', Exp.D, 'uni', 0));

ptbwords = reshape([ptbsw'; ptbew'], [], 1); % all strobe words

fprintf('PTB sent %d strobed words during this session\n', numel(ptbwords))
fprintf('Ephys recovered %d strobes (flip on channel 7)\n', nStrobes)

% Force synch using xcorr
ptbst = cellfun(@(x) x.STARTCLOCKTIME, Exp.D); % ptb trial starts
ptbet = cellfun(@(x) x.ENDCLOCKTIME, Exp.D); % ptb trial ends


%%
figure(1);clf
hold on
scatter(ptbst,ptbsw(:,1))
scatter(ptbst,ptbsw(:,2))
scatter(ptbst,ptbsw(:,3))
scatter(ptbst,ptbsw(:,4))
scatter(ptbst,ptbsw(:,5))
scatter(ptbst,ptbsw(:,6))

%% Just sync for your file
% toffset=1.6600655*10^9;
% ptbst=ptbst(ptbet>toffset);
% ptbet=ptbet(ptbet>toffset);
% ptbsw=ptbsw(ptbet>toffset,:);
% 
% figure(1);clf
% hold on
% scatter(ptbst,ptbsw(:,1))
% scatter(ptbst,ptbsw(:,2))
% scatter(ptbst,ptbsw(:,3))
% scatter(ptbst,ptbsw(:,4))
% scatter(ptbst,ptbsw(:,5))
% scatter(ptbst,ptbsw(:,6))
%%
ptbt = reshape([ptbst'; ptbet'], [], 1); % all strobe times from ptb


figure(1); clf
plot(diff(ptbt), 'o')
hold on



t1 = ts(find(diff(ts) > 2e-3) + 1);
figure(2); clf
T = max([t1-t1(1);ptbt-ptbt(1)]);
%T=max(ptbt); %max([ptbt-ptbt(1)]);

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

%% get best lag to match time series
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
bins = 0:bin_size:T;

[ecnt, ~, eid] = histcounts(t-t(1)+offset, bins);
[pcnt, ~, pid] = histcounts(ptbt-ptbt(1), bins);

figure(4); clf; set(gcf, 'Color', 'w')
plot(bins(1:end-1), pcnt); hold on
plot(bins(1:end-1), ecnt, '.')
%xlim([3000, 7000])
legend({'PTB sent Strobes', 'EPHYS recieved Strobes'})
xlabel('Time from first Strobe (seconds)')
title('Binned Strobe Times')


%%
potential_matches = find(ecnt==1 & pcnt==1);
fprintf('found %d potential matches at lag %02.3fs\n', numel(potential_matches), lags(id)*bin_size)

ematches = ismember(eid,potential_matches);
pmatches = ismember(pid, potential_matches);

assert(sum(ematches)==sum(pmatches), 'numer of matching timestamps does not match')
clf

fun = synchtime.align_clocks(ptbt(pmatches), t(ematches));
plot(ptbt(pmatches), t(ematches), 'o'); hold on
plot(xlim, fun(xlim), 'k')
xlabel('PTB clock')
ylabel('Ephys clock')

%% BUT! IS THIS THE SAME EPHYS CLOCK THAT THE SPIKES ARE ON???

%% Bypass
% fun = synchtime.align_clocks(ptbt, t);
Exp.ptb2Ephys = fun;

%% We really what to use inverse timing (epyhs to ptb), not implemented
% Because the ephys timing restarts each recording
ifun = synchtime.align_clocks(t(ematches), ptbt(pmatches));
plot(t(ematches), ptbt(pmatches), 'o'); hold on
plot(xlim, ifun(xlim), 'k')
ylabel('PTB clock')
xlabel('Ephys clock')
%% inverse timing (epyhs to ptb)
func=functions(fun);
w1=func.workspace{:}.w(1);
w2=func.workspace{:}.w(2);
infun= @(t) t*w1+w2;
Exp.Ephys2ptb = infun;

%% fake strobe times (necessary to synchronize the eye tracker)
% for iTrial = 1:numel(Exp.D)
% 
%     t = Exp.ptb2Ephys(Exp.D{iTrial}.STARTCLOCKTIME);
%     if t > 0
%         Exp.D{iTrial}.START_EPHYS = t;
%     end
% 
%     t = Exp.ptb2Ephys(Exp.D{iTrial}.ENDCLOCKTIME);
%     if t > 0
%         Exp.D{iTrial}.END_EPHYS = t;
%     end
% end

    
%% Baic marmoView import. Synchronize with Ephys if it exists
disp('BASIC MARMOVIEW IMPORT')
fid = 1;

%%

% cleanup eye position files
% edfFiles = dir(fullfile(dataPath, '*edf.mat'));
% for iFile = 1:numel(edfFiles)
%     edfname = fullfile(edfFiles(iFile).folder, edfFiles(iFile).name);
%     Dtmp = load(edfname);
%     Dtmp.edf = Dtmp.e;
%     save(edfname, '-v7.3', '-struct', 'Dtmp');
% end
% 
% disp('IMPORT EYE POSITION')

%% Import eye position signals
%dataPath = 'H:\.shortcut-targets-by-id\1te-Fna8YGaocWpO9rfNSoLtLeuzSNCPq\HuklabTreadmill\gru\g20220516'; %~/Google Drive/HuklabTreadmill/gru/g20210521/';
fid = 1;
[Exp, fig] = io.import_eye_position(Exp, dataPath, 'fid', fid, 'zero_mean', true, 'normalize', true);

%% QUICKSAVE!!
% Exp0 = Exp; % backup here
% save('b202200809_MT_quicksaveln198.mat','-v7.3')


C = calibGUI(Exp);
% GUIrefine breaks?  0.0063   -0.0071    2.1853    3.3857  -63.7712
%%
C0=C; %0.0065    0.0070   -0.0020    0.0806   -0.0549
%0.0068   -0.0070   -0.0020    0.0806   -0.0549
%%
%       scalex   scaley               shiftx shifty 
%WORKING
C.cmat=[-0.260    0.1600    0.4100      .275   -0.05];
C.cmat=[0.0065    0.0070   -0.0020    0.0806   -0.0549];
%FLIP Y: Match online calibration
C.cmat=[0.250    -0.1550    0.2000      -1.2   0.17];
C.cmat=[0.0057    -0.00610   -0.0020    0.076   -0.068];


%fig = io.checkCalibration(C.Exp);
% hold on
eyePos =get_eyepos(C);
plot(x1,y1,'r.')
hold on

plot(eyePos(1:20:end,1),eyePos(1:20:end,2),'k.')
hold off
 axis equal;xlim([-25 25]);ylim([-25 25]);
 
 %% Refine by hand
 C.cmat=[0.0061    -0.00690   -0.0020    0.080   -0.0550];
 C.cmat=[0.0064   -0.0069    1.9980    0.0800   -0.0550];
%%
eyePos =get_eyepos(C);
 figure(1);hold off
 plot(eyePos(1:50:end,1),eyePos(1:50:end,2),'r.')
xx=meshgrid(-20:5:20,-20:5:20);
yy=meshgrid(20:-5:-20,20:-5:-20)';
hold on
scatter(xx(:),yy(:),'k','linewidth',2);

axis equal;xlim([-12 12]);ylim([-12 12]);

%%

for i = 1:numel(Exp.D)
    
   eyeX = Exp.D{i}.eyeData
    Exp.D{i}.dx
end
%%
%% Check on facecal set
%C.Exp=Exp1;
figure(1); clf
plot(Exp1.vpx.smo(:,2)); hold on
figure(2); clf
plot(Exp1.vpx.smo(:,2)); hold on

%eyePos = C.refine_calibration(); % C.get_eyepos();
eyePos =get_eyepos(C);%Exp0.vpx.smo(:,2:3);

%FOR SOME REASON Y is flipped here
eyePos(:,2)=eyePos(:,2);


% get smo eyetrace
figure(1);
plot(eyePos(:,1))

vtt = Exp1.vpx.raw(:,1);
vpp = Exp1.vpx.raw(:,4);
vxxd = sgolayfilt(eyePos(:,1), 2, 9);
vyyd = sgolayfilt(eyePos(:,2), 2, 9);

vxx = medfilt1(vxxd, 5);
vyy = medfilt1(vyyd, 5);

vxx = imgaussfilt(vxx, 7);
vyy = imgaussfilt(vyy, 7);

vx = [0; diff(vxx)];
vy = [0; diff(vyy)];

vx = sgolayfilt(vx, 1, 3);
vy = sgolayfilt(vy, 1, 3);

% convert to d.v.a / sec
vx = vx * 1e3;
vy = vy * 1e3;

spd = hypot(vx, vy);
Exp1.vpx.smo = [vtt vxxd vyyd vpp vx vy spd];
%%
figure(1); 
hold off

plot(Exp1.vpx.smo(1:50:end,2),Exp1.vpx.smo(1:50:end,3),'.')
xx=meshgrid(-20:5:20,-20:5:20);
yy=meshgrid(20:-5:-20,20:-5:-20)';
hold on
scatter(xx(:),yy(:),'k');

axis equal;xlim([-10 10]);ylim([-10 10]);
%%
%Y has been flipped???
io.checkCalibration(Exp1)
%% Apply to other datasets

Exp=Exp0;
C.Exp = Exp0;%C.Exp;

%%
% S.CalibMat = C.cmat;

figure(1); clf
plot(Exp.vpx.smo(:,2)); hold on
%plot(Exp1.vpx.smo(:,2)); hold on

%eyePos = C.refine_calibration(); % C.get_eyepos();
eyePos =get_eyepos(C);%Exp0.vpx.smo(:,2:3);

%FOR SOME REASON Y is flipped here
eyePos(:,2)=eyePos(:,2);

% get smo eyetrace
plot(eyePos(:,1))

vtt = Exp.vpx.raw(:,1);
vpp = Exp.vpx.raw(:,4);
vxxd = sgolayfilt(eyePos(:,1), 2, 9);
vyyd = sgolayfilt(eyePos(:,2), 2, 9);

vxx = medfilt1(vxxd, 5);
vyy = medfilt1(vyyd, 5);

vxx = imgaussfilt(vxx, 7);
vyy = imgaussfilt(vyy, 7);

vx = [0; diff(vxx)];
vy = [0; diff(vyy)];

vx = sgolayfilt(vx, 1, 3);
vy = sgolayfilt(vy, 1, 3);

% convert to d.v.a / sec
vx = vx * 1e3;
vy = vy * 1e3;

spd = hypot(vx, vy);
Exp.vpx.smo = [vtt vxxd vyyd vpp vx vy spd];

plot(Exp.vpx.smo(:,2))
%%
figure(1); 
hold off

plot(Exp.vpx.smo(1:50:end,2),Exp.vpx.smo(1:50:end,3),'.')
xx=meshgrid(-20:5:20,-20:5:20);
yy=meshgrid(20:-5:-20,20:-5:-20)';
hold on
scatter(xx(:),yy(:),'k');

axis equal;xlim([-15 15]);ylim([-15 15]);

%%
save('b20220809_quicksaveln500_re.mat','-v7.3')
%%
% %%
% [eyePos, eyePos0] = io.getEyeCalibrationFromRaw(Exp, 'cmat', S.CalibMat);
% 
% Exp.vpx.smo(:,2:3) = eyePos0;
% fig = io.checkCalibration(Exp);
% fig.Position(1) = 50;
% title('Before Auto Refine')
% 
% Exp.vpx.smo(:,2:3) = eyePos;
% fig = io.checkCalibration(Exp);
% fig.Position(1) = fig.Position(1) + 100;
% title('After Auto Refine')
% % saveas(fig, fullfile(figDir, 'eye_calibration_from_file.pdf'))
% 
% %%
% if any(isnan(S.CalibMat))
%     disp('Running Manual Calibration')
%     Exp.FileTag = S.processedFileName;
%     C = calibGUI(Exp);
%     keyboard
%     S.CalibMat = C.cmat;
% end
% % redo calibration offline using FaceCal data
% [eyePos, eyePos0] = io.getEyeCalibrationFromRaw(Exp, 'cmat', S.CalibMat);
% 
% Exp.vpx.smo(:,2:3) = eyePos0;
% fig = io.checkCalibration(Exp);
% fig.Position(1) = 50;
% title('Before Auto Refine')
% 
% Exp.vpx.smo(:,2:3) = eyePos;
% fig = io.checkCalibration(Exp);
% fig.Position(1) = fig.Position(1) + 100;
% title('After Auto Refine')
% 
% f = uifigure('Position', [500 200 400 170]);
% selection = uiconfirm(f, 'Use Auto Refine Calibration?', 'Which Calibration', ...
%     'Options', {'Before', 'After'});
% close(f)
% 
% if strcmp(selection, 'After')
%     Exp.vpx.smo(:,2:3) = eyePos;
%     fig = io.checkCalibration(Exp);
%     title('After Auto Refine')
%     saveas(fig, fullfile(figDir, 'eye_calibration_redo_autorefine.pdf'))
% else
%     Exp.vpx.smo(:,2:3) = eyePos0;
%     fig = io.checkCalibration(Exp);
%     title('Before Auto Refine')
%     saveas(fig, fullfile(figDir, 'eye_calibration_redo_before.pdf'))
% end
% 
% close all
% % 
%%

% trial = 15;
% eyepos = Exp.D{trial}.eyeData(:,2:3);
% eyepos = Exp.vpx.raw(:,2:3);
% eyepos(:,2) = 1 - eyepos(:,2);
% 
% x = (eyepos(:,1)-Exp.D{trial}.c(1)) / (Exp.D{trial}.dx*Exp.S.pixPerDeg);
% y = (eyepos(:,2)-Exp.D{trial}.c(2)) / (Exp.D{trial}.dy*Exp.S.pixPerDeg);
% 
% figure(1); clf
% plot(x, y)


%%
% 


% Saccade processing:
% Perform basic processing of eye movements and saccades
Exp = saccadeflag.run_saccade_detection_cloherty(Exp, ...
    'ShowTrials', false,...
    'accthresh', 2e4,...
    'velthresh', 10,...
    'velpeak', 10,...
    'isi', 0.02);

% track invalid sampls
Exp.vpx.Labels(isnan(Exp.vpx.raw(:,2))) = 4;

validTrials = io.getValidTrials(Exp, 'Grating');
for iTrial = validTrials(:)'
    if ~isfield(Exp.D{iTrial}.PR, 'frozenSequence')
        Exp.D{iTrial}.PR.frozenSequence = false;
    end
end
        
fprintf(fid, 'Done importing session\n');

% Discont in fits?
plot(min(Exp.slist(:,1)):1000:max(Exp.slist(:,1)),Exp.ptb2Ephys(min(Exp.slist(:,1)):1000:max(Exp.slist(:,1))))

%% load spikes
% dataPath = '/home/declan/data/BrieMTsession/20220809/ephys/2022-08-09_17-31-35/Neuropix-PXI-101.0/Kilosort/';
%dataPath = '/home/declan/data/BrieMTsession/20220809/ephys/2022-08-09_12-44-58/experiment1/continuous/Neuropix-PXI-101.0/kilosort/';
Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp.st = Exp.osp.st;% + .6;
% move zero to the end
id = max(Exp.osp.clu)+1;

Exp.osp.clu(Exp.osp.clu==0) = id;
Exp.osp.cids = Exp.osp.cids2;
Exp.osp.cids(Exp.osp.cids==0) = id;

% Exp.osp.st = st_s.spikeTimes(:,1)+.6;
% Exp.osp.clu = st_s.spikeTimes(:,2);
% Exp.osp.depths = st_s.depths;
% Exp.osp.templates = st_s.templates;
% Exp.osp.cids = unique(Exp.osp.clu);

%%
% timingloc='/home/declan/data/BrieMTsession/20220809/ephys/2022-08-09_12-44-58/experiment1/continuous/Neuropix-PXI-101.0/';
% NPXtimestamps          = readNPY([timingloc 'timestamps.npy']);
NPXtimestamps = load([dataPath 'timestamps.mat']);
d_ts = double(NPXtimestamps.timestamps);
d_ts_s = d_ts/30000;
offset=min(d_ts_s)
Exp.osp.st0=Exp.osp.st;
Exp.osp.clu0=Exp.osp.clu;
%% Check timing!
figure(1); clf
plot.raster(Exp.osp.st0+offset, Exp.osp.clu0, 4); hold on
plot(t,strobes,'o')
%%
ptbt_=ptbt;
ptbt_=(ptbt_'+[0:0.001:0.005]');
ptbt_=ptbt_(:);

subplot(111)
plot(Exp.ptb2Ephys (ptbt_),ptbwords,'.')
hold on
plot(FPGAtimestamps.ts.timestamps_s, FPGAtimestamps.ts.fullWords, 'o'); hold off

%%
% Exp.osp.st0=Exp.osp.st;
% Exp.osp.st = Exp.osp.st0+offset;%Exp.Ephys2ptb(Exp.osp.st);
% Exp.spikeIds = Exp.osp.clu;



%% If happy, pick out good labelled units
Exp.osp = load(fullfile(dataPath, 'spkilo.mat'));
Exp.osp.clu0=Exp.osp.clu;
Exp.osp.st0 =Exp.osp.st;


%phy labels 0 = noise; 1 = MUA; 2 = Good; 3 = Unsorted
labels              = Exp.osp.cgs;
labelledClusters    = Exp.osp.cids;
goodclusters        = labelledClusters(labels==2);

Exp.osp.clu = Exp.osp.clu0(ismember(Exp.osp.clu0,goodclusters));

Exp.osp.st  = Exp.osp.st0(ismember(Exp.osp.clu0,goodclusters));
Exp.osp.st = Exp.osp.st+offset;%Exp.Ephys2ptb(Exp.osp.st);
Exp.spikeIds = Exp.osp.clu;

%% Don't know why we can't just add one to all, but leaving as is
% move zero to the end
id = max([Exp.osp.clu' Exp.osp.cids])+1;
Exp.osp.clu(Exp.osp.clu==0) = id;
Exp.osp.cids(Exp.osp.cids==0) = id;
% Exp.osp.cids = Exp.osp.cids2;
% Exp.osp.cids(Exp.osp.cids==0) = id;

%The above errors out, so just remaking
Exp.osp.cids = unique(Exp.osp.clu); 

Exp.spikeIds = Exp.osp.clu;
Exp.spikeTimes = Exp.osp.st;

%% Check timing! (again)
figure(1); clf
plot.raster(Exp.osp.st, (Exp.osp.clu-900)/8, 1)
ylabel('Unit ID')
xlabel('Time (seconds)')
hold on
plot(t,strobes,'o')

%% Trials against responses
inds = io.getValidTrials(Exp, 'BackImage');
spikeRate2=[];

spikeIds=unique(Exp.osp.clu);
NC=length(spikeIds);
for trial = 1:length(inds)
    thisTrial=inds(trial);
    %if length(Exp.D{thisTrial}.PR.NoiseHistory(:,4))>0
    %tt=Exp.D{thisTrial}.PR.NoiseHistory(:,1);
    tt=Exp.D{thisTrial}.STARTCLOCKTIME;
%     plot( tt-tt(1) ,Exp.D{thisTrial}.PR.NoiseHistory(:,4)) ; hold on
    
    %Convert start time (In PTB time) to ephys time
    t_=Exp.ptb2Ephys(tt); 
    t_=t_-7;
    t0=t_-40.0;
    t1=t_+40.0;
    idx=Exp.osp.st>t0&Exp.osp.st<=t1;
    subset_st=Exp.osp.st(idx);
    subset_clu=Exp.osp.clu(idx);
    
    
    tv=linspace(t0,t1,100);
    tsp=mode(diff(tv));
    for cc = 1:NC
        spikeRate2(:,cc,trial) = histcounts(subset_st(subset_clu==spikeIds(cc)), tv)./tsp;
    end
    %end
end



msr=mean(spikeRate2,3);
figure(1);clf


tv=linspace(t0-t_,t1-t_,100);
plot(tv(2:end),mean(msr'))
%plot.raster(time, neuron_id)

%% Exp.slist has saccadetimes
% Exp slist(1) is saccadetimes, try saccade triggered av, then try
% histogram of lags

spikeIds=unique(Exp.osp.clu);
NC=length(spikeIds);
stimes=Exp.ptb2Ephys(Exp.slist(:,1));
%stimes=stimes(stimes<5000&stimes>0);
nper=10;
spikeRate2 = zeros(99, NC, round(length(stimes)/nper));
%
for ii = 1:nper:length(stimes)
    tt=stimes(ii);
    %Convert start time (In PTB time) to ephys time
    t_=(tt)-7; 
    
    t0=t_-.1;
    t1=t_+.5;
    idx=Exp.osp.st>t0&Exp.osp.st<=t1;
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
imagesc(((msr-mean(msr))./std(msr))')
%
tv=linspace(t0-t_,t1-t_,100);
plot(tv(2:end),mean(msr'))

%% FPGA/IMEC offset
offset2=7;
Exp.osp.st=Exp.osp.st+offset2;
Exp.spikeTimes = Exp.osp.st;
%% Get spatial map

dotTrials = io.getValidTrials(Exp, 'Dots');
if ~isempty(dotTrials)
    

%     
    BIGROI = [-20 -20 20 20];
    binSize = .5;


%     BIGROI = [-5 -5 5 5];
%     eyePos = C.refine_calibration();
    eyePos = Exp.vpx.smo(:,2:3);
%     eyePos(:,1) = -eyePos(:,1);
    % eyePos(:,2) = -eyePos(:,2);

    
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

figure(1); clf
plot.raster(Exp.osp.st, Exp.osp.clu, 1); hold on
ylabel('Unit ID')
xlabel('Time (seconds)')
cmap = lines;
stims = {'Dots'}%, 'DriftingGrating', 'BackImage', 'FaceCal'};

for istim = 1:numel(stims)
    validTrials = io.getValidTrials(Exp, stims{istim});
    for iTrial = validTrials(:)'
        t0 = Exp.ptb2Ephys(Exp.D{iTrial}.STARTCLOCKTIME);
        t1 = Exp.ptb2Ephys(Exp.D{iTrial}.ENDCLOCKTIME);
        yd = ylim;
        h(istim) = fill([t0 t0 t1 t1], [yd(1) yd(2) yd(2) yd(1)], 'k', 'FaceColor', cmap(istim,:), 'FaceAlpha', .5, 'EdgeColor', 'none');
    end
end

%legend(h, stims)

%
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

win = [-1 15];
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
    
%% LAGS with subset
subset = [18, 21, 22, 42, 43, 52, 76, 92, 152, 161, 174, 259, 256, 566, 559, 464, 616, 626, 633, 714, 728,756,909,994,1025,1026,1064,1071,1266];
win = [-1 15];
stas = forwardCorrelation(Xstim, mean(RobsSpace(:,subset),2), win);
stas = stas / std(stas(:)) - mean(stas(:));

stas1 = stas()
wm = [min(stas(:)) max(stas(:))];
nlags = size(stas,1);
figure(1); clf
for ilag = 1:nlags
    subplot(2, ceil(nlags/2), ilag)
    imagesc(xax, yax, reshape(stas(ilag, :), opts.dims), wm)
end

%% plot subset
% subset = [33, 37, 68, 102, 298, 631, 634, 715, 718, 775, 832];
subset = [15, 46, 68, 37, 148, 151, 601, 602, 624, 631, 777, 991];
nsub = numel(subset);
[~, ccid] = ismember(subset, Exp.osp.cids(goodunits));

fig = ceil(1); clf
figure(fig); set(gcf, 'Color', 'w')

sx = ceil(sqrt(nsub));
sy = round(sqrt(nsub));

for cc = 1:nsub

    si = mod(cc, nperfig);  
    if si==0
        si = nperfig;
    end
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
    title(Exp.osp.cids(goodunits(cc)))
    grid on
end


%%
%offset=24050.948; %reading output on line 151 "fun = synchtime.align_clocks(ptbt(pmatches), t(ematches));'
% timingloc='/home/declan/data/BrieMTsession/20220809/ephys/2022-08-09_17-31-35/Neuropix-PXI-101.0/'
% FPGAtimestamps          = readNPY([timingloc 'timestamps.npy']);
% d_ts = double(FPGAtimestamps);
% d_ts_s = d_ts/30000;
% offset=min(d_ts_s);

% Exp.spikeTimes = Exp.osp.st+offset;%Exp.Ephys2ptb(Exp.osp.st);
% Exp.spikeIds = Exp.osp.clu;
%%
% convert to D struct format
D = io.get_drifting_grating_output(Exp, 'truncate', 0);

% Note: PLDAPS sessions and MarmoV5 will load into different formats
D0=D; % save
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

% find treadmill times that correspond to this session 
treadIx = D.treadTime > gratOnsets(1) & D.treadTime < gratOffsets(end);

%Downsample (heavily)
treadTime = D.treadTime(treadIx);
%Frames per grating
nn=round(median(diff(D.GratingOnsets))/median(diff(treadTime)));
treadTime = treadTime(1:nn:end);

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

%% PLOT SPIKE TIMES
spikeIds = unique(D.spikeIds(D.sessNumSpikes==sessionId));
NC = numel(spikeIds);
spikeRate = zeros(numel(treadTime), NC);

bs = diff(treadTime);
for cc = 1:NC
    spikeRate(:,cc) = [0 histcounts(D.spikeTimes(D.spikeIds==spikeIds(cc)), treadTime)./bs'];
end
goodunits = find(mean(spikeRate)>0.005);
spikeRate = spikeRate ./ max(spikeRate); % normalize for visualization

treadSpeed = D.treadSpeed(treadIx);
treadSpeed1 = treadSpeed(1:nn:end);

% Not sure whats going on with these -ves
treadSpeed(treadSpeed<0)=0;
% Mean treadmill over nn
rows=floor(length(treadSpeed)/nn);
treadSpeed=nanmean(reshape((treadSpeed(1:rows*nn)),nn,rows))';
if rem(length(treadSpeed),nn)
    treadSpeed(rows+1)=nanmean((treadSpeed(rows*nn:end)));
end
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
imagesc(treadTime, 1:numel(goodunits), spikeRate(:,goodunits)'); hold on
for i = 1:numel(onsets)
    fill(treadTime([onsets(i) onsets(i) offsets(i) offsets(i)]), [ylim fliplr(ylim)], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
end
xlim(treadTime([1 end]))
colormap(1-gray)
title('Spikes')
ylabel('Unit #')   

% PLOT TREADMILL RUNNING SPEED
subplot(3,1,3) % tread speed
plot(treadTime, treadSpeed , 'k'); hold on
xlabel('Time (s)')
ylabel('Speed (cm/s)')

for i = 1:numel(onsets)
    fill(treadTime([onsets(i) onsets(i) offsets(i) offsets(i)]), [ylim fliplr(ylim)], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
end
xlim(treadTime([1 end]))

title('Treadmill')
trmax = prctile(D.treadSpeed(treadIx), 99);

%%
save('b20220809_quicksaveln825_re_huh.mat','-v7.3')

%% Checking for ephys/FPGA offsets via saccade triggered average
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
    
    tmin=-2.80;
    tplus=2.8;
    t0=t_+tmin;
    t1=t_+tplus;
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
imagesc(((msr-mean(msr))./std(msr))')
%
tv=linspace(tmin,tplus,100);
plot(tv(2:end),msr)
%%
 plot(tv(2:end),mean(msr'))

%remove largest units( swamping mean)


plot(tv(2:end),mean(msr(:,(mean(msr)>.025))'))
%%
figure(111); clf
subplot(8,1,[1 2 3 4 5 6 7]) % plot spike count

% imagesc(1:length(treadTime), 1:numel(goodunits), spikeRate(:,goodunits)');
imagesc(1:length(treadTime), 1:size(spikeRate,2), spikeRate');
%


c       = 256*ones(256,3);
c(:,1)  = 256:-1:1;
c(:,2)  = 256:-1:1;

c(1:128,3)  = 256*ones(128,1);
c(1:128,1)  = 256:-2:2;
c(1:128,2)  = 256:-2:2;
c(129:256,3)  = 256:-1:129;
c(129:256,1)  = 1:128;
c(129:256,2)  = 1:128;


colormap(uint8(c))

% hold on
% for i = 1:numel(onsets)
%     fill([onsets(i) onsets(i) offsets(i) offsets(i)], [ylim fliplr(ylim)], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
% end
xlim([1 length(treadTime)]);set(gca,'visible','off')

title('Spikes')
ylabel('Unit #')   
%plot.fixfigure(gcf, 12, [8 11])

% PLOT TREADMILL RUNNING SPEED
subplot(8,1,8) % tread speed
plot(1:length(treadTime), treadSpeed ,'Color', [.75 .75 .75],'LineWidth',1.5);% hold on
xlabel('Trials')
ylabel('Speed (cm/s)')

% for i = 1:numel(onsets)
%     fill(([onsets(i) onsets(i) offsets(i) offsets(i)]), [ylim fliplr(ylim)], 'r', 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
% end
 xlim([1 length(treadTime)]);set(gca,'visible','off')

title('Treadmill')
trmax = prctile(D.treadSpeed(treadIx), 99);

%plot.fixfigure(gcf, 12, [8 11])
%%
figDir='/home/declan/data/BrieMTsession/20220809/out/';
saveas(gcf, fullfile(figDir, 'grating_sess_birdseye3.pdf'))

%% play movie
vidObj = VideoWriter(fullfile(figDir, 'session.mp4'), 'MPEG-4');
vidObj.Quality = 100;
vidObj.FrameRate = 5;
open(vidObj);


% play a movie of the session
figure(gcf)
t0 = treadTime(1);
win = 25;
for t = 1:500
    xd = [t0 t0 + win];
    for i = 1:3
        subplot(3,1,i)
        xlim(xd)
    end
    ylim([0 trmax])
    drawnow
    t0 = t0 + win/10;
    if t0 > treadTime(end)
        break
    end
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
end

close(vidObj);

%% Firing rates
D.treadSpeed = max(D.treadSpeed, 0);
D.sessNumTread = double(D.treadTime > 0);

iix = ismember(D.spikeIds, Exp.osp.cids(goodunits));
D.spikeIds(~iix) = 9999;

%%
sessionId = 1;
[StimDir, spksb, runningSpeed, Dstat] = bin_population(D, sessionId);

%% PSTH
stimids = unique(StimDir);
nstim = numel(stimids);

Rtuning = zeros(nstim, Dstat.NLags, Dstat.NCells);
for istim = 1:nstim
    iix = StimDir==stimids(istim);
    Rtuning(istim,:,:) = squeeze(mean(spksb(iix,:,:)));
end

% running?
isrunning = mean(runningSpeed,2) > 3;
Rtc = zeros(nstim, 2, Dstat.NCells);
tix = Dstat.lags > .05 & Dstat.lags < .9;
R = squeeze(mean(spksb(:,tix,:), 2));
for istim = 1:nstim
    iix = StimDir==stimids(istim);

    Rtc(istim,1,:) = mean(R(iix & isrunning,:));
    Rtc(istim,2,:) = mean(R(iix & ~isrunning,:));
end

%% check out individual cell
%subset = [31,32,108,144,358, 598, 648, 892, 692,905,918,934, 949, 964, 965, 980, 1174, 1183, 1055];
%subset = [6, 31,34,37,44,107,114,124,143,147,164,227,276,384, 488, 584, 630, 761, 805,819];
% subset = [33, 37, 68, 102, 298, 631, 634, 715, 718, 775, 832];
subset = [15, 46, 68, 37, 148, 151, 601, 602, 624, 631, 777, 991];
nsub = numel(subset);
[~, ccid] = ismember(subset, Exp.osp.cids(goodunits));

%% or...
ccid=Exp.osp.cids(goodunits);
cc=0
%%
cc = cc + 1;
    
I = squeeze(stas(:,:,ccid(cc)));
I = (I - mean(I(:)) )/ std(I(:));

[bestlag,j] = find(I==max(I(:)));
I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
I = mean(I);

figure(1); clf
subplot(2,2,1)
xax = opts.xax/Exp.S.pixPerDeg;
yax = opts.yax/Exp.S.pixPerDeg;
imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2])
colormap(plot.coolwarm)
title(Exp.osp.cids(goodunits(cc)))

i = find(I == max(I(:)));
subplot(2,2,2)
tkernel = stas(:, i, ccid(cc));
plot(win(1):win(2), tkernel, 'k', 'Linewidth', 2)
%%
subplot(2,2,3)
for istim = 1:nstim
    plot(Dstat.lags, imboxfilt(squeeze(Rtuning(istim,:,ccid(cc))),11)*60, 'Color', cmap(istim,:)); hold on
end

subplot(2,2,4)
plot(stimids, Rtc(:,1,ccid(cc))*60, '-o'); hold on
plot(stimids, Rtc(:,2,ccid(cc))*60, '-o');


%%
%% plotting
sx = 15;
NC = size(stas,3);
sy = ceil(NC/sx);

figure(1); clf; set(gcf, 'Color', 'w')
ax = plot.tight_subplot(sy, sx, 0.01, 0.01);

for cc = 1:NC

    set(gcf, 'currentaxes', ax(cc))

    I = squeeze(stas(:,:,cc));
    I = (I - mean(I(:)) )/ std(I(:));

    [bestlag,j] = find(I==max(I(:)));
    I = I(min(max(bestlag(1)+[-1 0 1], 1), nlags),:);
    I = mean(I);
    xax = opts.xax/Exp.S.pixPerDeg;
    yax = opts.yax/Exp.S.pixPerDeg;
    imagesc(xax, yax, imgaussfilt(reshape(I, opts.dims), 1), [-2, 2])
    colormap(plot.coolwarm)
    title(Exp.osp.cids(goodunits(cc)))
end
    


%%
figure(2); clf

cmap = plot.viridis(nstim);
ax = plot.tight_subplot(sy, sx, 0.01, 0.01);

for cc = 1:Dstat.NCells
    set(gcf, 'currentaxes', ax(cc))
    for istim = 1:nstim
        plot(Dstat.lags, imboxfilt(squeeze(Rtuning(istim,:,cc)),11), 'Color', cmap(istim,:)); hold on
    end
    axis off
end


figure(3); clf
ax = plot.tight_subplot(sy, sx, 0.01, 0.01);
for cc = 1:Dstat.NCells
    set(gcf, 'currentaxes', ax(cc))
    imagesc(Dstat.lags, 1:nstim, squeeze(Rtuning(:,:,cc)))
    axis off
end

figure(4); clf
ax = plot.tight_subplot(sy, sx, 0.01, 0.01);
for cc = 1:Dstat.NCells
    set(gcf, 'currentaxes', ax(cc))
    plot(stimids, Rtc(:,1,cc)*60, '-o'); hold on
    plot(stimids, Rtc(:,2,cc)*60, '-o');
end
%%



%% Unit by Unit simple analysis
thresh = 1;
nboot = 100;
D.sessNumEye=1;
D.sessNumEye= ones(size(D.eyeLabels));

unitList = unique(D.spikeIds);
NC = numel(unitList);

corrRho = zeros(NC,1);
corrPval = zeros(NC,1);

frBaseR = zeros(NC,3);
frBaseS = zeros(NC,3);

frStimR = zeros(NC,3);
frStimS = zeros(NC,3);

for cc = 1:NC
    unitId = unitList(cc);

    [stimDir0, robs, runSpd, opts] = bin_ssunit(D, unitId, 'win', [-.2 .1]);
    
%     pause
    goodIx = getStableRange(sum(robs,2), 'plot', false);
%    
    %stimDir0 = stimDir0(goodIx);
    stimDir = stimDir0{1};
    stimDir = stimDir(goodIx);
    
    robs = robs(goodIx,:);
    runSpd = runSpd{1}(goodIx,:);
    
    iix = opts.lags < 0;
    frbase = sum(robs(:,iix),2) / (max(opts.lags(iix)) - min(opts.lags(iix)));
    
    spd = mean(runSpd,2);
    
    [corrRho(cc), corrPval(cc)] = corr(spd, frbase, 'type', 'Spearman');
    
    runTrials = find(spd > thresh);
    statTrials = find(abs(spd) < 1);
    mixTrials = [runTrials; statTrials];
    
    nrun = numel(runTrials);
    nstat = numel(statTrials);
    
    n = min(nrun, nstat);
    if n>0
    frBaseR(cc,:) = prctile(mean(frbase(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
    frBaseS(cc,:) = prctile(mean(frbase(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
    
    iix = opts.lags > 0.04 & opts.lags < opts.lags(end)-.15;
    frstim = sum(robs(:,iix),2) / (max(opts.lags(iix)) - min(opts.lags(iix)));
    
    frStimR(cc,:) = prctile(mean(frstim(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
    frStimS(cc,:) = prctile(mean(frstim(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
    end
end

%% plot some outcomes

incBaseIx = find(frBaseR(:,2) > frBaseS(:,3));
decBaseIx = find(frBaseR(:,2) < frBaseS(:,1));

incStimIx = find(frStimR(:,2) > frStimS(:,3));
decStimIx = find(frStimR(:,2) < frStimS(:,1));

figure(1); clf
set(gcf, 'Color', 'w')
ms = 4;
cmap = lines;
subplot(1,2,1)
plot(frBaseS(:,[2 2])', frBaseR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
plot(frBaseS(:,[1 3])', frBaseR(:,[2 2])', 'Color', .5*[1 1 1])
plot(frBaseS(:,2), frBaseR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
plot(frBaseS(incBaseIx,2), frBaseR(incBaseIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
plot(frBaseS(decBaseIx,2), frBaseR(decBaseIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))

xlabel('Stationary')
ylabel('Running')
title('Baseline Firing Rate')

set(gca, 'Xscale', 'log', 'Yscale', 'log')
plot(xlim, xlim, 'k')

subplot(1,2,2)
plot(frStimS(:,[2 2])', frStimR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
plot(frStimS(:,[1 3])', frStimR(:,[2 2])', 'Color', .5*[1 1 1])
plot(frStimS(:,2), frStimR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
plot(frStimS(incStimIx,2), frStimR(incStimIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
plot(frStimS(decStimIx,2), frStimR(decStimIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))

xlabel('Stationary Firing Rate')
ylabel('Running Firing Rate')
title('Stim-driven firing rate')

set(gca, 'Xscale', 'log', 'Yscale', 'log')
plot(xlim, xlim, 'k')

nIncBase = numel(incBaseIx);
nDecBase = numel(decBaseIx);

nIncStim = numel(incStimIx);
nDecStim = numel(decStimIx);

modUnits = unique([incBaseIx; decBaseIx; incStimIx; decStimIx]);
nModUnits = numel(modUnits);

fprintf('%d/%d (%02.2f%%) increased baseline firing rate\n', nIncBase, NC, 100*nIncBase/NC)
fprintf('%d/%d (%02.2f%%) decreased baseline firing rate\n', nDecBase, NC, 100*nDecBase/NC)

fprintf('%d/%d (%02.2f%%) increased stim firing rate\n', nIncStim, NC, 100*nIncStim/NC)
fprintf('%d/%d (%02.2f%%) decreased stim firing rate\n', nDecStim, NC, 100*nDecStim/NC)

fprintf('%d/%d (%02.2f%%) total modulated units\n', nModUnits, NC, 100*nModUnits/NC)

[pvalStim, ~, sStim] = signrank(frStimS(:,2), frStimR(:,2));
[pvalBase, ~, sBase] = signrank(frBaseS(:,2), frBaseR(:,2));

fprintf('Wilcoxon signed rank test:\n')
fprintf('Baseline rates: p = %02.10f\n', pvalBase)
fprintf('Stim-driven rates: p = %02.10f\n', pvalStim)

good = ~(frBaseR(:,2)==0 | frBaseS(:,2)==0);

m = geomean(frBaseR(good,2)./frBaseS(good,2));
ci = bootci(nboot, @geomean, frBaseR(good,2)./frBaseS(good,2));

fprintf("geometric mean baseline firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good)) 
%%
good = ~(frStimR(:,2)==0 | frStimS(:,2)==0);
m = geomean(frStimR(good,2)./frStimS(good,2));
% ci = bootci(nboot, @geomean, frStimR(:,2)./frStimS(:,2));
ci = nan(1,2);
fprintf("geometric mean stim-driven firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)




%% Do direction / orientation decoding


% example session to explore
Dstat = decode_stim(D, 1, 'figDir', figDir);

%% Bootstrapped empirical analyses and tuning curve fits

fpath = getpref('FREEVIEWING', 'HUKLAB_DATASHARE');

% D.sessNumTread = ones(size(D.treadTime))*1;

D.treadSpeed(D.treadSpeed < 0) = nan;
stat = tuning_empirical(D, 'binsize', 10e-3, ...
    'runningthresh', 3, ...
    'nboot', 500, ...
    'seed', 1234);  



%% Plot all tuning curves
fitS = stat.TCfitS;
fitR = stat.TCfitR;
NC = numel(fitS);

sx = ceil(sqrt(NC));
sy = round(sqrt(NC));


figure(1); clf
ax = plot.tight_subplot(sx, sy, 0.02);
for cc = 1:NC
    if min(fitS(cc).numTrials, fitR(cc).numTrials) < 50
        continue
    end
%     
%     if fitS(cc).llrpval > 0.05 && fitR(cc).llrpval > 0.05
%         continue
%     end
    fprintf("Unit %d/%d\n", cc, NC)

    set(gcf, 'currentaxes', ax(cc))
    
    cmap = lines;
    % STATIONARY
    h = errorbar(fitS(cc).thetas, fitS(cc).tuningCurve, fitS(cc).tuningCurveSE, '-o', 'Color', cmap(1,:));
    h.MarkerSize = 2;
    h.MarkerFaceColor = cmap(1,:);
    h.CapSize = 0;
    hold on
    if fitS(cc).llrpval > 0.05
        plot(xlim, mean(fitS(cc).tuningCurve)*[1 1], 'Color', cmap(1,:))
    else
        plot(linspace(0, 360, 100), fitS(cc).tuningFun(linspace(0, 360, 100)), 'Color', cmap(1,:))
    end
    
    % RUNNING
    h = errorbar(fitR(cc).thetas, fitR(cc).tuningCurve, fitR(cc).tuningCurveSE, '-o', 'Color', cmap(2,:));
    h.MarkerSize = 2;
    h.MarkerFaceColor = cmap(2,:);
    h.CapSize = 0;
    hold on
    if fitR(cc).llrpval > 0.05
        plot(xlim, mean(fitR(cc).tuningCurve)*[1 1], 'Color', cmap(2,:))
    else
        plot(linspace(0, 360, 100), fitR(cc).tuningFun(linspace(0, 360, 100)), 'Color', cmap(2,:))
    end

    set(gca, 'XTick', [], 'YTick', [])
%     set(gca, 'XTick', 0:180:360)
    xlim([0 360])
    text(10, .1*max(ylim), sprintf('%d', cc), 'fontsize', 5)
%     axis off
    title(sprintf('Unit: %d', Exp.osp.cids(cc)))
end

plot.fixfigure(gcf, 10, [sx sy]*2, 'offsetAxes', false);
saveas(gcf, fullfile(figDir, 'tuning_curves.pdf'))

%% If happy, save processe file
% pathparts = regexp(dataPath, '/', 'split');
% Exp.FileTag = [pathparts{end-1} '_' pathparts{end}(2:end) '.mat'];
Exp.FileTag = 'b20220809_re.mat'
save(fullfile(HUKDATA, 'HuklabTreadmill/VSS' , Exp.FileTag), '-v7.3', '-struct', 'Exp')
