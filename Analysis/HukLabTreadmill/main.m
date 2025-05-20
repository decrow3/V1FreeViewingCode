%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory

%% Basic summary of session running
thresh = 3;
nboot = 100;
subjects = {'gru', 'brie'};%, 'allen'};
nsubjs = numel(subjects);

gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/gratings/';

for i = 1:nsubjs
    subject = subjects{i};
    fprintf('Subject [%s], run threshold %02.2f cm/s\n', subject, thresh)
    D = load_subject(subject,gratingpath);
    
    fprintf("%d Unique units\n", numel(unique(D.spikeIds)))
    treadSpeed = D.treadSpeed;
    runix = treadSpeed > thresh;
    
    m = mean(runix)*100;
    n = numel(runix);
    booti = mean(runix(randi(n, [n nboot])));
    
    ci = prctile(booti, [2.5 97.5])*100;
    fprintf("%02.2f [%02.2f,%02.2f] %% of time running\n", m, ci(1), ci(2))
    
    m = mean(treadSpeed(runix));
    runinds = find(runix);
    n = sum(runix);
    booti = mean(treadSpeed(runinds(randi(n, [n nboot]))));
    ci = prctile(booti, [2.5 97.5]);
    fprintf("%02.2f [%02.2f,%02.2f] cm/s mean running speed\n", m, ci(1), ci(2))
    
end

%% Step 2: Load a super session file
gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/gratings/';
subject = 'gru';
D = load_subject(subject,gratingpath);

%% Step 2.1: plot one session
% This just demonstrates how to index into individual sessions and get
% relevant data

sessionId = 17; % pick a session id

% index into grating times from that session
gratIx = D.sessNumGratings==sessionId;
gratOnsets = D.GratingOnsets(gratIx);
gratOffsets = D.GratingOffsets(gratIx);
gratDir = D.GratingDirections(gratIx);

% find treadmill times that correspond to this session 
treadIx = D.treadTime > gratOnsets(1) & D.treadTime < gratOffsets(end);
treadTime = D.treadTime(treadIx);


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

% PLOT SPIKE TIMES
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




%% play movie
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
        t0 = treadTime(1);
    end
    pause(0.1)
end
lw=0;up=13;
bins=(up-lw)/nbins;
%% load session
gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/gratings/';
subject = 'gru';
D = load_subject(subject,gratingpath);

%% To avoid double counting, we need to remove sessions 15 and 17, which are the same as 14 and 16 but at a different part of the probe
D.eyeLabels=D.eyeLabels((D.sessNumEye~=17 &D.sessNumEye~=15));
D.eyePos=D.eyePos((D.sessNumEye~=17 &D.sessNumEye~=15),:);
D.eyeTime=D.eyeTime((D.sessNumEye~=17 &D.sessNumEye~=15),:);
D.sessNumEye=D.sessNumEye((D.sessNumEye~=17 &D.sessNumEye~=15),:);
D.spikeIds=D.spikeIds((D.sessNumSpikes~=17 & D.sessNumSpikes~=15),:);
D.spikeTimes=D.spikeTimes((D.sessNumSpikes~=17 & D.sessNumSpikes~=15),:);
D.sessNumSpikes=D.sessNumSpikes((D.sessNumSpikes~=17 & D.sessNumSpikes~=15),:);
D.treadSpeed=D.treadSpeed((D.sessNumTread~=17 & D.sessNumTread~=15),:);
D.treadTime=D.treadTime((D.sessNumTread~=17 & D.sessNumTread~=15),:);
D.sessNumTread=D.sessNumTread((D.sessNumTread~=17 & D.sessNumTread~=15),:);
D.GratingContrast=D.GratingContrast((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
D.GratingDirections=D.GratingDirections((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
D.GratingFrequency=D.GratingFrequency((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
D.GratingOffsets=D.GratingOffsets((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
D.GratingOnsets=D.GratingOnsets((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
D.GratingSpeeds=D.GratingSpeeds((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
D.sessNumGratings=D.sessNumGratings((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
%%
D.eyePos(D.eyePos(:,3)==0,:)=nan;
D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
Eyestat.(subject) = do_eye_analyses(D);
%Stat.(subject) = do_spike_count_analyses(D);
%% Regen eyelabels??
subject = 'brie';
D = load_subject(subject,gratingpath);
D.eyePos(D.eyePos(:,3)==0,:)=nan;
D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
%%
Eyestat.(subject) = do_eye_analyses(D);
%%
nrun = numel(Eyestat.(subject).runTrials);
nstat = numel(Eyestat.(subject).statTrials);
n=min([nrun nstat]);


%Histograms
figure(4);clf
    cmap = getcolormap(subject, false);
subplot(221)
lw=0;up=10;
bins=lw:(up-lw)/(up-1):up;
hold off
% histogram(Eyestat.(subject).SacHz(Eyestat.(subject).runTrials(randperm(nrun, n ))),bins,'Normalization','probability'  )
histogram(Eyestat.(subject).SacHz(Eyestat.(subject).runTrials),bins,'Normalization','probability'  )
hold on
histogram(Eyestat.(subject).SacHz(Eyestat.(subject).statTrials),bins,'Normalization','probability'  )
legend('Running','Stationary')
xlabel('Saccade Rate (Hz)')

nbins=20; %for plotting
subplot(222)
lw=0;up=20;
bins=lw:(up-lw)/(nbins-1):up;
hold off
histogram(Eyestat.(subject).MSac(Eyestat.(subject).runTrials),bins,'Normalization','probability'  )
hold on
histogram(Eyestat.(subject).MSac(Eyestat.(subject).statTrials),bins,'Normalization','probability'  )
legend('Running','Stationary')
xlabel('Saccade Mag')

subplot(223)
lw=800;up=3500;
bins=lw:(up-lw)/(nbins-1):up;
hold off
histogram(Eyestat.(subject).eyeSize(Eyestat.(subject).runTrials),bins,'Normalization','probability'  )
hold on
histogram(Eyestat.(subject).eyeSize(Eyestat.(subject).statTrials),bins,'Normalization','probability'  )
legend('Running','Stationary')
xlabel('Pupil Size')

subplot(224)
lw=0;up=50;
bins=lw:(up-lw)/(nbins-1):up;
hold off
histogram(Eyestat.(subject).varXY(Eyestat.(subject).runTrials),bins,'Normalization','probability'  )
hold on
histogram(Eyestat.(subject).varXY(Eyestat.(subject).statTrials),bins,'Normalization','probability'  )

legend('Running','Stationary')
xlabel('Eye pos variance')

% %% Anova pairs (need to correct for multiple?)
% SacHzRtrials=Eyestat.(subject).SacHz(Eyestat.gru.runTrials);
% SacHzStrials=Eyestat.(subject).SacHz(Eyestat.gru.statTrials);
% [P,ANOVATAB,STATS]=anova1([SacHzRtrials; SacHzStrials],[ones(length(SacHzRtrials),1);2*ones(length(SacHzStrials),1)]);
% 
% MSacRtrials=Eyestat.(subject).MSac(Eyestat.gru.runTrials);
% MSacStrials=Eyestat.(subject).MSac(Eyestat.gru.statTrials);
% [P,ANOVATAB,STATS]=anova1([MSacRtrials; MSacStrials],[ones(length(MSacRtrials),1);2*ones(length(MSacStrials),1)]);
% 
% eyeSizeRtrials=Eyestat.(subject).eyeSize(Eyestat.gru.runTrials);
% eyeSizeStrials=Eyestat.(subject).eyeSize(Eyestat.gru.statTrials);
% [P,ANOVATAB,STATS]=anova1([eyeSizeRtrials; eyeSizeStrials],[ones(length(eyeSizeRtrials),1);2*ones(length(eyeSizeStrials),1)]);
% 
% varXYRtrials=Eyestat.(subject).varXY(Eyestat.gru.runTrials);
% varXYStrials=Eyestat.(subject).varXY(Eyestat.gru.statTrials);
% [P,ANOVATAB,STATS]=anova1([varXYRtrials; varXYStrials],[ones(length(varXYRtrials),1);2*ones(length(varXYStrials),1)]);

%%
fprintf("Mean saccade frequency during running is %02.3f Hz [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).SacHzR(2), Eyestat.(subject).SacHzR(1), Eyestat.(subject).SacHzR(3), length(Eyestat.(subject).nSac))
fprintf("Mean saccade frequency during stationary is %02.3f Hz [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).SacHzS(2), Eyestat.(subject).SacHzS(1), Eyestat.(subject).SacHzS(3), length(Eyestat.(subject).nSac))
fprintf("Two-way anova saccade frequency in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(Eyestat.(subject).runTrials),length(Eyestat.(subject).statTrials),(Eyestat.(subject).SacHzP))

fprintf("Mean saccade magnitude during running is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).SacMagR(2), Eyestat.(subject).SacMagR(1), Eyestat.(subject).SacMagR(3), length(Eyestat.(subject).nSac))
fprintf("Mean saccade magnitude during stationary is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).SacMagS(2), Eyestat.(subject).SacMagS(1), Eyestat.(subject).SacMagS(3), length(Eyestat.(subject).nSac))
fprintf("Two-way anova saccade magnitude in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(Eyestat.(subject).runTrials),length(Eyestat.(subject).statTrials),(Eyestat.(subject).SacMagP))

fprintf("Mean pupil size during running is %02.3f units? [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).PupSzR(2), Eyestat.(subject).PupSzR(1), Eyestat.(subject).PupSzR(3), length(Eyestat.(subject).nSac))
fprintf("Mean pupil size during stationary is %02.3f units? [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).PupSzS(2), Eyestat.(subject).PupSzS(1), Eyestat.(subject).PupSzS(3), length(Eyestat.(subject).nSac))
fprintf("Two-way anova pupil size in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(Eyestat.(subject).runTrials),length(Eyestat.(subject).statTrials),(Eyestat.(subject).PupSzP))

fprintf("Mean eye position variance during running is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).EyeVarR(2), Eyestat.(subject).EyeVarR(1), Eyestat.(subject).EyeVarR(3), length(Eyestat.(subject).nSac))
fprintf("Mean eye position variance during stationary is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", Eyestat.(subject).EyeVarS(2), Eyestat.(subject).EyeVarS(1), Eyestat.(subject).EyeVarS(3), length(Eyestat.(subject).nSac))
fprintf("Two-way anova eye position variance in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(Eyestat.(subject).runTrials),length(Eyestat.(subject).statTrials),(Eyestat.(subject).EyeVarP))

%%
fprintf("Pearson correlation coefficient of saccade frequency with running is %02.3f [p = %d] (n=%d trials) \n", Eyestat.(subject).corrSacRho, Eyestat.(subject).corrSacPval,length(Eyestat.(subject).nSac))
fprintf("Pearson correlation coefficient of saccade magnitude with running is %02.3f [p = %d] (n=%d trials) \n", Eyestat.(subject).corrSacMagRho, Eyestat.(subject).corrSacMagPval,length(Eyestat.(subject).nSac))
fprintf("Pearson correlation coefficient of pupil size with running is %02.3f [p = %d] (n=%d trials) \n", Eyestat.(subject).corrSizeRho, Eyestat.(subject).corrSizePval,length(Eyestat.(subject).nSac))
fprintf("Pearson correlation coefficient of eye position variance with running is %02.3f [p = %d] (n=%d trials) \n", Eyestat.(subject).corrVarXYRho, Eyestat.(subject).corrVarXYPval,length(Eyestat.(subject).nSac))

%%
clear rstatTrials
%match nearest static to running
runTrials=Eyestat.(subject).runTrials;
statTrials=Eyestat.(subject).statTrials;
nSac=Eyestat.(subject).nSac;
SacHz=Eyestat.(subject).SacHz;
nrun=length(runTrials);
for ii=1:nrun
nrun=length(runTrials);
    [~,idx]=min(abs(statTrials-runTrials(ii)));
    rstatTrials(ii)=statTrials(idx);
end
   
%%
figure
subplot(2,2,1); hold off
plot(SacHz(runTrials), SacHz(rstatTrials),'.'); hold on

hold off
histogram((SacHz(runTrials)./ SacHz(rstatTrials)),10.^[-1.25:0.25:1.25],'Normalization','probability'  )
set(gca, 'Xscale', 'log')
title('Sacade rate ratio (running/stationary)')
%

subplot(2,2,2); hold off
plot(Eyestat.(subject).MSac(runTrials), Eyestat.(subject).MSac(rstatTrials),'.'); hold on
plot(xlim, xlim, 'k')
%
hold off
histogram((Eyestat.(subject).MSac(runTrials)./ Eyestat.(subject).MSac(rstatTrials)),10.^[-1.25:0.25:1.25],'Normalization','probability'  )
set(gca, 'Xscale', 'log')
title('Sacade magnitude ratio (running/stationary)')

subplot(2,2,3); hold off
plot(Eyestat.(subject).eyeSize(runTrials), Eyestat.(subject).eyeSize(rstatTrials),'.'); hold on

hold off
histogram((Eyestat.(subject).eyeSize(runTrials)./ Eyestat.(subject).eyeSize(rstatTrials)),10.^[-0.225:0.05:.225],'Normalization','probability'  )
set(gca, 'Xscale', 'log')
title('Eye size ratio (running/stationary)')
%

subplot(2,2,4); hold off
plot(Eyestat.(subject).varXY(runTrials), Eyestat.(subject).varXY(rstatTrials),'.'); hold on
plot(xlim, xlim, 'k')
%
hold off
histogram((Eyestat.(subject).varXY(runTrials)./ Eyestat.(subject).varXY(rstatTrials)),10.^[-2.5:0.5:2.5],'Normalization','probability'  )
set(gca, 'Xscale', 'log')
title('Eye pos variance ratio (running/stationary)')

%% Combine subjects (for eye stats vs running here)
subjects = {'gru','brie'};%, 'allen'};
nsubjs = numel(subjects);
SacHzR_all=[];
SacHzS_all=[];
MSacR_all=[];
MSacS_all=[];
eyeSizeR_all=[];
eyeSizeS_all=[];
varXYR_all=[];
varXYS_all=[];

for isubj = 1:nsubjs
    %
    subject = subjects{isubj};
    nrun(isubj) = numel(Eyestat.(subject).runTrials);
    nstat(isubj)  = numel(Eyestat.(subject).statTrials);
    n(isubj) =min([nrun nstat]);
    
    SacHzR_all=[SacHzR_all; Eyestat.(subject).SacHz(Eyestat.(subject).runTrials)];
    SacHzS_all=[SacHzS_all; Eyestat.(subject).SacHz(Eyestat.(subject).statTrials)];
    MSacR_all =[MSacR_all; Eyestat.(subject).MSac(Eyestat.(subject).runTrials)];
    MSacS_all =[MSacS_all; Eyestat.(subject).MSac(Eyestat.(subject).statTrials)];
    eyeSizeR_all=[eyeSizeR_all; Eyestat.(subject).eyeSize(Eyestat.(subject).runTrials)];
    eyeSizeS_all=[eyeSizeS_all; Eyestat.(subject).eyeSize(Eyestat.(subject).statTrials)];
    varXYR_all=[varXYR_all; Eyestat.(subject).varXY(Eyestat.(subject).runTrials)];
    varXYS_all=[varXYS_all; Eyestat.(subject).varXY(Eyestat.(subject).statTrials)];

end
nrun=sum(nrun);
nstat=sum(nstat);
%%
%Histograms
figure(4);clf
    cmap = getcolormap(subject, false);
subplot(131)
lw=0;up=8;
% bins=lw:(up-lw)/(up-1):up;
nbins=11;
bins=lw:(up-lw)/(nbins-1):up;
hold off
histogram(SacHzR_all,bins,'Normalization','probability'  )
hold on
histogram(SacHzS_all,bins,'Normalization','probability'  )
legend('Running','Stationary')
xlabel('Saccade Rate (Hz)')

nbins=20; %for plotting
subplot(132)
lw=0;up=20;
bins=lw:(up-lw)/(nbins-1):up;
hold off
histogram(MSacR_all,bins,'Normalization','probability'  )
hold on
histogram(MSacS_all,bins,'Normalization','probability'  )
legend('Running','Stationary')
xlabel('Saccade Mag (deg)')

subplot(133)
lw=0.5*100;up=175;
bins=lw:(up-lw)/(nbins-1):up;
hold off
histogram(100*eyeSizeR_all(eyeSizeR_all>0)./mean([eyeSizeR_all(eyeSizeR_all>0);eyeSizeS_all(eyeSizeS_all>0)]),bins,'Normalization','probability'  )
hold on
histogram(100*eyeSizeS_all(eyeSizeS_all>0)./mean([eyeSizeR_all(eyeSizeR_all>0);eyeSizeS_all(eyeSizeS_all>0)]),bins,'Normalization','probability'  )
legend('Running','Stationary')
xlabel('Pupil Size (%)')

% subplot(224)
% lw=0;up=50;
% bins=lw:(up-lw)/(nbins-1):up;
% hold off
% histogram(varXYR_all,bins,'Normalization','probability'  )
% hold on
% histogram(varXYS_all,bins,'Normalization','probability'  )
% 
% legend('Running','Stationary')
% xlabel('Eye pos variance')

%% Stats for combined data
run_thresh=1;
nboot=100;


SacHz=[];
MSac=[];
eyeSize=[];
varXY=[];
spd=[];

for isubj = 1:nsubjs
    %
    nrun(isubj) = numel(Eyestat.(subject).runTrials);
    nstat(isubj)  = numel(Eyestat.(subject).statTrials);
    n(isubj) =min([nrun nstat]);
    subject = subjects{isubj};

    SacHz=[SacHz; Eyestat.(subject).SacHz];    
    MSac =[MSac; Eyestat.(subject).MSac];
    eyeSize=[eyeSize; Eyestat.(subject).eyeSize];
    varXY=[varXY; Eyestat.(subject).varXY];

    spd=[spd; Eyestat.(subject).spd];
end


%Behaviour correlates with speed
[corrSacRho, corrSacPval] = corr(spd(~isnan(spd)), SacHz(~isnan(spd)), 'type', 'Spearman');
[corrSizeRho, corrSizePval] = corr(spd(~isnan(spd)&~isnan(eyeSize)), eyeSize(~isnan(spd)&~isnan(eyeSize)), 'type', 'Spearman');
[corrVarXYRho, corrVarXYPval] = corr(spd(~isnan(spd)&~isnan(varXY)), varXY(~isnan(spd)&~isnan(varXY)), 'type', 'Spearman');
[corrMSacRho, corrMSacPval] = corr(spd(~isnan(spd)&~isnan(MSac)), MSac(~isnan(spd)&~isnan(MSac)), 'type', 'Spearman');

runTrials = find(spd > run_thresh);
statTrials = find(spd < run_thresh);
mixTrials = [runTrials; statTrials];


nrun = numel(runTrials);
nstat = numel(statTrials);
n = min(nrun, nstat);

SacHzR = prctile(mean(SacHz(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
SacHzS = prctile(mean(SacHz(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
SacHzAll = prctile(mean(SacHz((randi(nstat, [n nboot])))), [2.5 50 97.5]);

PupSzR = prctile(nanmean(eyeSize(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
PupSzS = prctile(nanmean(eyeSize(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
PupSzAll = prctile(nanmean(eyeSize((randi(nstat, [n nboot])))), [2.5 50 97.5]);

EyeVarR = prctile(nanmean(varXY(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
EyeVarS = prctile(nanmean(varXY(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
EyeVarAll = prctile(nanmean(varXY((randi(nstat, [n nboot])))), [2.5 50 97.5]);


% 
% runsactri=runTrials(ismember(runTrials,sacmagtri));
% statsactri=statTrials(ismember(statTrials,sacmagtri));
% mixscatri=mixTrials(ismember(mixTrials,sacmagtri));
% 
% nrunsac = numel(runsactri);
% nstatsac = numel(statsactri);
% nmixsac = numel(mixscatri);
% 
% SacMagR = prctile(nanmean(MSac(runsactri(randi(nrunsac, [n opts.nboot])))), [2.5 50 97.5]);
% SacMagS = prctile(nanmean(MSac(statsactri(randi(nstatsac, [n opts.nboot])))), [2.5 50 97.5]);
% SacMagAll = prctile(nanmean(MSac(mixscatri(randi(nmixsac, [n opts.nboot])))), [2.5 50 97.5]);

SacMagR = prctile(nanmean(MSac(runTrials(randi(nrun, [n nboot])))), [2.5 50 97.5]);
SacMagS = prctile(nanmean(MSac(statTrials(randi(nstat, [n nboot])))), [2.5 50 97.5]);
SacMagAll = prctile(nanmean(MSac((randi(n, [n nboot])))), [2.5 50 97.5]);


%%

group=[repmat({'Running'},nrun,1); repmat({'Stationary'},nstat,1)];
[SacHzP,SacHzANOVATAB,SacHzSTATS]=anova1([SacHzR_all; SacHzS_all],group);
[PupSzP,PupSzANOVATAB,PupSzSTATS]=anova1([eyeSizeR_all; eyeSizeS_all],group);
[EyeVarP,EyeVarANOVATAB,EyeVarSTATS]=anova1([varXYR_all; varXYS_all],group);
[SacMagP,SacMagANOVATAB,SacMagSTATS]=anova1([MSacR_all; MSacS_all],group);

%%

fprintf("Mean saccade frequency during running is %02.3f Hz [%02.3f, %02.3f] (n=%d trials) \n", SacHzR(2), SacHzR(1), SacHzR(3), length(SacHz))
fprintf("Mean saccade frequency during stationary is %02.3f Hz [%02.3f, %02.3f] (n=%d trials) \n", SacHzS(2), SacHzS(1), SacHzS(3), length(SacHz))
fprintf("Two-way anova saccade frequency in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(SacHzP))

fprintf("Mean saccade magnitude during running is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", SacMagR(2), SacMagR(1), SacMagR(3), length(SacHz))
fprintf("Mean saccade magnitude during stationary is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", SacMagS(2), SacMagS(1), SacMagS(3), length(SacHz))
fprintf("Two-way anova saccade magnitude in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(SacMagP))

fprintf("Mean pupil size during running is %02.3f units? [%02.3f, %02.3f] (n=%d trials) \n", PupSzR(2), PupSzR(1), PupSzR(3), length(SacHz))
fprintf("Mean pupil size during stationary is %02.3f units? [%02.3f, %02.3f] (n=%d trials) \n", PupSzS(2), PupSzS(1), PupSzS(3), length(SacHz))
fprintf("Two-way anova pupil size in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(PupSzP))

fprintf("Mean eye position variance during running is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", EyeVarR(2), EyeVarR(1), EyeVarR(3), length(SacHz))
fprintf("Mean eye position variance during stationary is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", EyeVarS(2), EyeVarS(1), EyeVarS(3), length(SacHz))
fprintf("Two-way anova eye position variance in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(EyeVarP))

%%
fprintf("Pearson correlation coefficient of saccade frequency with running is %02.3f [p = %d] (n=%d trials) \n", corrSacRho, corrSacPval,length(SacHz))
fprintf("Pearson correlation coefficient of saccade magnitude with running is %02.3f [p = %d] (n=%d trials) \n", corrMSacRho, corrMSacPval,length(SacHz))
fprintf("Pearson correlation coefficient of pupil size with running is %02.3f [p = %d] (n=%d trials) \n", corrSizeRho, corrSizePval,length(SacHz))
fprintf("Pearson correlation coefficient of eye position variance with running is %02.3f [p = %d] (n=%d trials) \n", corrVarXYRho, corrVarXYPval,length(SacHz))

%% Reload Gru's data with sessions 15 and 17 back in for cell-by-cell analysis
subject = 'gru';
D = load_subject(subject,gratingpath);
D.eyePos(D.eyePos(:,3)==0,:)=nan;
D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
%%

Stat.(subject) = do_spike_count_eye_analyses(D);
%%
Stat = struct();

subjects = {'gru', 'brie'};
nsubjs = numel(subjects);
for isubj = 1:nsubjs
    subject = subjects{isubj};
    D = load_subject(subject,gratingpath);

    Stat.(subject) = do_spike_count_eye_analyses(D);
end

%%
figdir = '/home/huklab/Documents/NPX_pilot/V1FreeViewingCode/Figures';
nboot = 100;



% [207, 179, 144; 195, 97, 66]
% cmap = lines;

% colors = {[(cmap(1,:) + [1 1 1]/2)/2; cmap(1,:)]*255, ...
%     [(cmap(5,:) + [1 1 1]/2)/2; cmap(5,:)]*255, ...
%     [(cmap(2,:) + [1 1 1]/2)/2; cmap(2,:)]*255};

subjects = {'gru','brie'};%, 'allen'};
nsubjs = numel(subjects);


% This is a messy quick way to combine over the loop
    frBaseRall = [];
    frBaseSall = [];
    frStimRall = [];
    frStimSall = [];
    frPrefRall = [];
    frPrefSall = [];

    incBaseIxall = []; 
    decBaseIxall = []; 
    notSigBaseIxall = [];

    incStimIxall = []; 
    decStimIxall = []; 
    notSigStimIxall = []; 

    incPrefIxall = []; 
    decPrefIxall = []; 
    notSigPrefIxall = []; 

        % Eye data
        SacHzRall = []; 
        SacHzSall = []; 
        
        incSacHzIxall = []; 
        decSacHzIxall = []; 
        
        %Saccade Magnitude (SacMg)
        SacMgRall = []; 
        SacMgSall = []; 

        incSacMgIxall = [];
        decSacMgIxall = []; 

        %Pupil size (PupSz)
        PupSzRall=[];
        PupSzSall=[];
        incPupSzIxall=[];
        decPupSzIxall=[];
        
        % Eye Variance (EyeVar)
        EyeVarRall=[];
        EyeVarSall=[];
        incEyeVarIxall=[];
        decEyeVarIxall=[];

        NCtillnow=0;


for isubj = 1:nsubjs
    %
    subject = subjects{isubj};
    cmap = getcolormap(subject, false);

    fprintf('***************************\n\n')
    fprintf('***************************\n\n')
    fprintf('%s\n', subject)

    frBaseS = Stat.(subject).frBaseS;
    frBaseR = Stat.(subject).frBaseR;
    frStimR = Stat.(subject).frStimR;
    frStimS = Stat.(subject).frStimS;
    frPrefR = Stat.(subject).frPrefR;
    frPrefS = Stat.(subject).frPrefS;

    good = ~isnan(Stat.(subject).meanrate); %~(frStimS(:,2) < 1 | frStimR(:,2) < 1);
    frBaseR = frBaseR(good,:);
    frBaseS = frBaseS(good,:);
    frStimR = frStimR(good,:);
    frStimS = frStimS(good,:);
    frPrefR = frPrefR(good,:);
    frPrefS = frPrefS(good,:);
    NC = size(frBaseR,1);

    incBaseIx = find(frBaseR(:,2) > frBaseS(:,3));
    decBaseIx = find(frBaseR(:,2) < frBaseS(:,1));
    notSigBaseIx = setdiff( (1:NC)', [incBaseIx; decBaseIx]);

    incStimIx = find(frStimR(:,2) > frStimS(:,3));
    decStimIx = find(frStimR(:,2) < frStimS(:,1));
    notSigStimIx = setdiff( (1:NC)', [incStimIx; decStimIx]);

    incPrefIx = find(frPrefR(:,2) > frPrefS(:,3));
    decPrefIx = find(frPrefR(:,2) < frPrefS(:,1));
    notSigPrefIx = setdiff( (1:NC)', [incPrefIx; decPrefIx]);

    figure(isubj); clf
    set(gcf, 'Color', 'w')
    ms = 2;
    
    subplot(1,3,1)
    plot(frBaseS(:,[2 2])', frBaseR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
    plot(frBaseS(:,[1 3])', frBaseR(:,[2 2])', 'Color', .5*[1 1 1])

    plot(frBaseS(incBaseIx,2), frBaseR(incBaseIx,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frBaseS(decBaseIx,2), frBaseR(decBaseIx,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frBaseS(notSigBaseIx,2), frBaseR(notSigBaseIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
    xlabel('Stationary')
    ylabel('Running')
    title('Baseline Firing Rate')

     set(gca, 'Xscale', 'log', 'Yscale', 'log')
    plot([.1 100], [.1 100], 'k')
    xlim([0.1 100])
    ylim([0.1 100])

    subplot(1,3,2)
    plot(frStimS(:,[2 2])', frStimR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
    plot(frStimS(:,[1 3])', frStimR(:,[2 2])', 'Color', .5*[1 1 1])

    plot(frStimS(incStimIx,2), frStimR(incStimIx,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frStimS(decStimIx,2), frStimR(decStimIx,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frStimS(notSigStimIx,2), frStimR(notSigStimIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))

    xlabel('Stationary Firing Rate')
    ylabel('Running Firing Rate')
    title('Stim-driven firing rate')

     set(gca, 'Xscale', 'log', 'Yscale', 'log')
    plot([.1 100], [.1 100], 'k')
    xlim([0.1 100])
    ylim([0.1 100])

    subplot(1,3,3)
    plot(frPrefS(:,[2 2])', frPrefR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
    plot(frPrefS(:,[1 3])', frPrefR(:,[2 2])', 'Color', .5*[1 1 1])

    plot(frPrefS(incPrefIx,2), frPrefR(incPrefIx,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frPrefS(decPrefIx,2), frPrefR(decPrefIx,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frPrefS(notSigPrefIx,2), frPrefR(notSigPrefIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))

    xlabel('Stationary Firing Rate')
    ylabel('Running Firing Rate')
    title('Pref-driven firing rate')

     set(gca, 'Xscale', 'log', 'Yscale', 'log')
    plot([.1 100], [.1 100], 'k')
    xlim([0.1 100])
    ylim([0.1 100])


    nIncBase = numel(incBaseIx);
    nDecBase = numel(decBaseIx);

    nIncStim = numel(incStimIx);
    nDecStim = numel(decStimIx);

    nIncPref = numel(incPrefIx);
    nDecPref = numel(decPrefIx);

    modUnits = unique([incBaseIx; decBaseIx; incStimIx; decStimIx; incPrefIx; decPrefIx]);
    nModUnits = numel(modUnits);

    fprintf('%d/%d (%02.2f%%) increased baseline firing rate\n', nIncBase, NC, 100*nIncBase/NC)
    fprintf('%d/%d (%02.2f%%) decreased baseline firing rate\n', nDecBase, NC, 100*nDecBase/NC)

    fprintf('%d/%d (%02.2f%%) increased stim firing rate\n', nIncStim, NC, 100*nIncStim/NC)
    fprintf('%d/%d (%02.2f%%) decreased stim firing rate\n', nDecStim, NC, 100*nDecStim/NC)

    fprintf('%d/%d (%02.2f%%) increased Pref firing rate\n', nIncPref, NC, 100*nIncPref/NC)
    fprintf('%d/%d (%02.2f%%) decreased Pref firing rate\n', nDecPref, NC, 100*nDecPref/NC)

    fprintf('%d/%d (%02.2f%%) total modulated units\n', nModUnits, NC, 100*nModUnits/NC)

    [pvalStim, ~, sStim] = signrank(frStimS(:,2), frStimR(:,2));
    [pvalBase, ~, sBase] = signrank(frBaseS(:,2), frBaseR(:,2));
    [pvalPref, ~, sPref] = signrank(frPrefS(:,2), frPrefR(:,2));

    fprintf('Wilcoxon signed rank test:\n')
    fprintf('Baseline rates: p = %02.10f\n', pvalBase)
    fprintf('Stim-driven rates: p = %02.10f\n', pvalStim)
    fprintf('Pref-driven rates: p = %02.10f\n', pvalPref)

    good = ~(frBaseR(:,2)==0 | frBaseS(:,2)==0);
    m = geomean(frBaseR(good,2)./frBaseS(good,2));
    ci = bootci(nboot, @geomean, frBaseR(good,2)./frBaseS(good,2));

    fprintf("geometric mean baseline firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good))

    good = find(~(frStimS(:,2)==0 | frStimS(:,2)==0));
    good = good(~isinf(log(frStimR(good,2)./frStimS(good,2))));
    m = geomean(frStimR(good,2)./frStimS(good,2));
    ci = bootci(nboot, @geomean, frStimR(good,2)./frStimS(good,2));

    fprintf("geometric mean stim-driven firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)

    good = find(~(frPrefS(:,2)==0 | frPrefR(:,2)==0 | isnan(frPrefS(:,2)) | isnan(frPrefR(:,2))));
    good = good(~isinf(log(frPrefR(good,2)./frPrefS(good,2))));
    m = geomean(frPrefR(good,2)./frPrefS(good,2));
    ci = bootci(nboot, @geomean, frPrefR(good,2)./frPrefS(good,2));

    fprintf("geometric mean Pref-driven firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)



    %plot.formatFig(gcf, [4 3], 'nature')
    %saveas(gcf, fullfile(figdir, sprintf('rate_compare_%s.pdf', subject)))

    good = ~isnan(Stat.(subject).meanrate); 
    
    % EYE DATA
    % Saccade Rate (SacHz)

    SacHzR = Stat.(subject).SacHzR(good,:);
    SacHzS = Stat.(subject).SacHzS(good,:);
        
        incSacHzIx = find(SacHzR(:,2) > SacHzS(:,3));
        decSacHzIx = find(SacHzR(:,2) < SacHzS(:,1));
        
        %Saccade Magnitude (SacMg)
    SacMgR = Stat.(subject).SacMgR(good,:);
    SacMgS = Stat.(subject).SacMgS(good,:);

        incSacMgIx = find(SacMgR(:,2) > SacMgS(:,3));
        decSacMgIx = find(SacMgR(:,2) < SacMgS(:,1));
        
        figure(isubj*2); clf
        set(gcf, 'Color', 'w')
        ms = 4;
        cmap = lines;
        subplot(2,2,1)
        plot(SacHzS(:,[2 2])', SacHzR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
        plot(SacHzS(:,[1 3])', SacHzR(:,[2 2])', 'Color', .5*[1 1 1])
        plot(SacHzS(:,2), SacHzR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
        plot(SacHzS(incSacHzIx,2), SacHzR(incSacHzIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
        plot(SacHzS(decSacHzIx,2), SacHzR(decSacHzIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))
        
        xlabel('Stationary')
        ylabel('Running')
        title('Saccade rate (Hz)')
        
        set(gca, 'Xscale', 'log', 'Yscale', 'log')
        plot(xlim, xlim, 'k')
        
        subplot(2,2,2)
        plot(SacMgS(:,[2 2])', SacMgR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
        plot(SacMgS(:,[1 3])', SacMgR(:,[2 2])', 'Color', .5*[1 1 1])
        plot(SacMgS(:,2), SacMgR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
        plot(SacMgS(incSacMgIx,2), SacMgR(incSacMgIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
        plot(SacMgS(decSacMgIx,2), SacMgR(decSacMgIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))
        
        xlabel('Stationary')
        ylabel('Running')
        title('Saccade magnitude')
        
        set(gca, 'Xscale', 'log', 'Yscale', 'log')
        plot(xlim, xlim, 'k')
        
        nIncBase = numel(incSacHzIx);
        nDecBase = numel(decSacHzIx);
        
        nIncStim = numel(incSacMgIx);
        nDecStim = numel(decSacMgIx);
        
        modUnits = unique([incSacHzIx; decSacHzIx; incSacMgIx; decSacMgIx]);
        nModUnits = numel(modUnits);
        
        fprintf('%d/%d (%02.2f%%) increased saccade rate\n', nIncBase, NC, 100*nIncBase/NC)
        fprintf('%d/%d (%02.2f%%) decreased saccade rate\n', nDecBase, NC, 100*nDecBase/NC)
        
        
        %fprintf('%d/%d (%02.2f%%) total modulated units\n', nModUnits, NC, 100*nModUnits/NC)
        
        [pvalStim, ~, sStim] = signrank(SacMgS(:,2), SacMgR(:,2));
        [pvalBase, ~, sBase] = signrank(SacHzS(:,2), SacHzR(:,2));
        
        %fprintf('Wilcoxon signed rank test:\n')
        fprintf('Saccade rates: p = %02.10f\n', pvalBase)
        
        
        good = ~(SacHzR(:,2)==0 | SacHzS(:,2)==0);
        
        m = geomean(SacHzR(good,2)./SacHzS(good,2));
        ci = bootci(nboot, @geomean, SacHzR(good,2)./SacHzS(good,2));
        
        fprintf("geometric mean saccade rate (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good)) 
        fprintf('%d/%d (%02.2f%%) increased saccade mag\n', nIncStim, NC, 100*nIncStim/NC)
        fprintf('%d/%d (%02.2f%%) decreased saccade mag\n', nDecStim, NC, 100*nDecStim/NC)
        fprintf('Saccade mag: p = %02.10f\n', pvalStim)

        m = geomean(SacMgR(:,2)./SacMgS(:,2));
%        ci = bootci(nboot, @geomean, SacMgR(:,2)./SacMgS(:,2));
        
        fprintf("geometric mean saccade mag (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)
        
        good = ~isnan(Stat.(subject).meanrate); 
        
        %Pupil size (PupSz)
        PupSzR = Stat.(subject).PupSzR(good,:);
        PupSzS = Stat.(subject).PupSzS(good,:);
        incPupSzIx = find(PupSzR(:,2) > PupSzS(:,3));
        decPupSzIx = find(PupSzR(:,2) < PupSzS(:,1));
        
        % Eye Variance (EyeVar)
        EyeVarR = Stat.(subject).EyeVarR(good,:);
        EyeVarS = Stat.(subject).EyeVarS(good,:);        
        incEyeVarIx = find(EyeVarR(:,2) > EyeVarS(:,3));
        decEyeVarIx = find(EyeVarR(:,2) < EyeVarS(:,1));
        
        
        ms = 4;
        cmap = lines;
        subplot(2,2,3)
        plot(PupSzS(:,[2 2])', PupSzR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
        plot(PupSzS(:,[1 3])', PupSzR(:,[2 2])', 'Color', .5*[1 1 1])
        plot(PupSzS(:,2), PupSzR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
        plot(PupSzS(incPupSzIx,2), PupSzR(incPupSzIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
        plot(PupSzS(decPupSzIx,2), PupSzR(decPupSzIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))
        
        xlabel('Stationary')
        ylabel('Running')
        title('Pupil Size')
        
        set(gca, 'Xscale', 'log', 'Yscale', 'log')
        plot(xlim, xlim, 'k')
        
        subplot(2,2,4)
        plot(EyeVarS(:,[2 2])', EyeVarR(:,[1 3])', 'Color', .5*[1 1 1]); hold on
        plot(EyeVarS(:,[1 3])', EyeVarR(:,[2 2])', 'Color', .5*[1 1 1])
        plot(EyeVarS(:,2), EyeVarR(:,2), 'o', 'Color', cmap(1,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(1,:))
        plot(EyeVarS(incEyeVarIx,2), EyeVarR(incEyeVarIx,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
        plot(EyeVarS(decEyeVarIx,2), EyeVarR(decEyeVarIx,2), 'o', 'Color', cmap(3,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(3,:))
        
        xlabel('Stationary ')
        ylabel('Running ')
        title('Eye position variance')
        
        set(gca, 'Xscale', 'log', 'Yscale', 'log')
        plot(xlim, xlim, 'k')
        
        nIncBase = numel(incPupSzIx);
        nDecBase = numel(decPupSzIx);
        
        nIncStim = numel(incEyeVarIx);
        nDecStim = numel(decEyeVarIx);
        
        modUnits = unique([incPupSzIx; decPupSzIx; incEyeVarIx; decEyeVarIx]);
        nModUnits = numel(modUnits);
        
        fprintf('%d/%d (%02.2f%%) increased pupil size\n', nIncBase, NC, 100*nIncBase/NC)
        fprintf('%d/%d (%02.2f%%) decreased pupil size\n', nDecBase, NC, 100*nDecBase/NC)
        

        
        %fprintf('%d/%d (%02.2f%%) total modulated units\n', nModUnits, NC, 100*nModUnits/NC)
        
        [pvalStim, ~, sStim] = signrank(EyeVarS(:,2), EyeVarR(:,2));
        [pvalBase, ~, sBase] = signrank(PupSzS(:,2), PupSzR(:,2));
        
        %fprintf('Wilcoxon signed rank test:\n')
        fprintf('Pupil size: p = %02.10f\n', pvalBase)
        
        
        good = ~(PupSzR(:,2)==0 | PupSzS(:,2)==0);
        
        m = geomean(PupSzR(good,2)./PupSzS(good,2));
        ci = bootci(nboot, @geomean, PupSzR(good,2)./PupSzS(good,2));
        
        fprintf("geometric mean pupil size ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good)) 
 
        fprintf('%d/%d (%02.2f%%) increased eye pos variance\n', nIncStim, NC, 100*nIncStim/NC)
        fprintf('%d/%d (%02.2f%%) decreased eye pos variance\n', nDecStim, NC, 100*nDecStim/NC)
        fprintf('Eye position variance: p = %02.10f\n', pvalStim)


        m = geomean(EyeVarR(:,2)./EyeVarS(:,2));
        ci = bootci(nboot, @geomean, EyeVarR(:,2)./EyeVarS(:,2));
        
        fprintf("geometric mean eye pos variance ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)


% This is a messy quick way to combine over the loop
    frBaseRall = [frBaseRall; frBaseR];
    frBaseSall = [frBaseSall; frBaseS];
    frStimRall = [frStimRall;frStimR];
    frStimSall = [frStimSall;frStimS];
    frPrefRall = [frPrefRall;frPrefR];
    frPrefSall = [frPrefSall;frPrefS];

    
    incBaseIxall = [incBaseIxall;incBaseIx+NCtillnow]; 
    decBaseIxall = [decBaseIxall;decBaseIx+NCtillnow]; 
    notSigBaseIxall = [notSigBaseIxall;notSigBaseIx+NCtillnow]; 

    incStimIxall = [incStimIxall;incStimIx+NCtillnow];  
    decStimIxall = [decStimIxall;decStimIx+NCtillnow];  
    notSigStimIxall = [notSigStimIxall;notSigStimIx+NCtillnow];  

    incPrefIxall = [incPrefIxall;incPrefIx+NCtillnow];  
    decPrefIxall = [decPrefIxall;decPrefIx+NCtillnow];  
    notSigPrefIxall = [notSigPrefIx;notSigPrefIxall]; 

        % Eye data
        SacHzRall = [SacHzRall;SacHzR]; 
        SacHzSall = [SacHzSall;SacHzS]; 
        
        incSacHzIxall = [incSacHzIxall;incSacHzIx+NCtillnow];  
        decSacHzIxall = [decSacHzIxall;decSacHzIx+NCtillnow];  
        
        %Saccade Magnitude (SacMg)
        SacMgRall = [SacMgRall;SacMgR]; 
        SacMgSall = [SacMgSall;SacMgS]; 

        incSacMgIxall = [incSacMgIxall;incSacMgIx+NCtillnow]; 
        decSacMgIxall = [decSacMgIxall;decSacMgIx+NCtillnow];  

        %Pupil size (PupSz)
        PupSzRall=[PupSzRall;PupSzR];
        PupSzSall=[PupSzSall;PupSzS];
        incPupSzIxall=[incPupSzIxall;incPupSzIx+NCtillnow]; 
        decPupSzIxall=[decPupSzIxall;decPupSzIx+NCtillnow]; 
        
        % Eye Variance (EyeVar)
        EyeVarRall=[EyeVarRall;EyeVarR];
        EyeVarSall=[EyeVarSall;EyeVarS];
        incEyeVarIxall=[incEyeVarIxall;incEyeVarIx+NCtillnow]; 
        decEyeVarIxall=[decEyeVarIxall;decEyeVarIx+NCtillnow]; 

        NCtillnow=NCtillnow+NC;


end
%% combined running effect
figure(1);clf
NC=NCtillnow;
        set(gcf, 'Color', 'w')
        ms = 2;
        cmap = getcolormap('gru', false);
   subplot(1,3,1)
    plot(frBaseSall(:,[2 2])', frBaseRall(:,[1 3])', 'Color', .5*[1 1 1]); hold on
    plot(frBaseSall(:,[1 3])', frBaseRall(:,[2 2])', 'Color', .5*[1 1 1])

    plot(frBaseSall(incBaseIxall,2), frBaseRall(incBaseIxall,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frBaseSall(decBaseIxall,2), frBaseRall(decBaseIxall,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frBaseSall(notSigBaseIxall,2), frBaseRall(notSigBaseIxall,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))
    xlabel('Stationary')
    ylabel('Running')
    title('Baseline Firing Rate')

     set(gca, 'Xscale', 'log', 'Yscale', 'log')
    plot([.1 100], [.1 100], 'k')
    xlim([0.1 100])
    ylim([0.1 100])

    subplot(1,3,2)
    plot(frStimSall(:,[2 2])', frStimRall(:,[1 3])', 'Color', .5*[1 1 1]); hold on
    plot(frStimSall(:,[1 3])', frStimRall(:,[2 2])', 'Color', .5*[1 1 1])

    plot(frStimSall(incStimIxall,2), frStimRall(incStimIxall,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frStimSall(decStimIxall,2), frStimRall(decStimIxall,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frStimSall(notSigStimIxall,2), frStimRall(notSigStimIxall,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))

    xlabel('Stationary Firing Rate')
    ylabel('Running Firing Rate')
    title('Stim-driven firing rate')

     set(gca, 'Xscale', 'log', 'Yscale', 'log')
    plot([.1 100], [.1 100], 'k')
    xlim([0.1 100])
    ylim([0.1 100])

    subplot(1,3,3)
    plot(frPrefSall(:,[2 2])', frPrefRall(:,[1 3])', 'Color', .5*[1 1 1]); hold on
    plot(frPrefSall(:,[1 3])', frPrefRall(:,[2 2])', 'Color', .5*[1 1 1])

    plot(frPrefSall(incPrefIxall,2), frPrefRall(incPrefIxall,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frPrefSall(decPrefIxall,2), frPrefRall(decPrefIxall,2), 'o', 'Color', cmap(6,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(6,:))
    plot(frPrefSall(notSigPrefIxall,2), frPrefRall(notSigPrefIxall,2), 'o', 'Color', cmap(2,:), 'MarkerSize', ms, 'MarkerFaceColor', cmap(2,:))

    xlabel('Stationary Firing Rate')
    ylabel('Running Firing Rate')
    title('Pref-driven firing rate')

     set(gca, 'Xscale', 'log', 'Yscale', 'log')
    plot([.1 100], [.1 100], 'k')
    xlim([0.1 100])
    ylim([0.1 100])


    nIncBase = numel(incBaseIxall);
    nDecBase = numel(decBaseIxall);

    nIncStim = numel(incStimIxall);
    nDecStim = numel(decStimIxall);

    nIncPref = numel(incPrefIxall);
    nDecPref = numel(decPrefIxall);

    modUnits = unique([incBaseIxall; decBaseIxall; incStimIxall; decStimIxall; incPrefIxall; decPrefIxall]);
    nModUnits = numel(modUnits);

    fprintf('\n-----COMBINED DATA-----\n')
    fprintf('%d/%d (%02.2f%%) increased baseline firing rate\n', nIncBase, NC, 100*nIncBase/NC)
    fprintf('%d/%d (%02.2f%%) decreased baseline firing rate\n', nDecBase, NC, 100*nDecBase/NC)

    fprintf('%d/%d (%02.2f%%) increased stim firing rate\n', nIncStim, NC, 100*nIncStim/NC)
    fprintf('%d/%d (%02.2f%%) decreased stim firing rate\n', nDecStim, NC, 100*nDecStim/NC)

    fprintf('%d/%d (%02.2f%%) increased Pref firing rate\n', nIncPref, NC, 100*nIncPref/NC)
    fprintf('%d/%d (%02.2f%%) decreased Pref firing rate\n', nDecPref, NC, 100*nDecPref/NC)

    fprintf('%d/%d (%02.2f%%) total modulated units\n', nModUnits, NC, 100*nModUnits/NC)

    [pvalStim, ~, sStim] = signrank(frStimSall(:,2), frStimRall(:,2));
    [pvalBase, ~, sBase] = signrank(frBaseSall(:,2), frBaseRall(:,2));
    [pvalPref, ~, sPref] = signrank(frPrefSall(:,2), frPrefRall(:,2));

    fprintf('Wilcoxon signed rank test:\n')
    fprintf('Baseline rates: p = %02.10f\n', pvalBase)
    fprintf('Stim-driven rates: p = %02.10f\n', pvalStim)
    fprintf('Pref-driven rates: p = %02.10f\n', pvalPref)

    good = ~(frBaseRall(:,2)==0 | frBaseSall(:,2)==0);
    m = geomean(frBaseRall(good,2)./frBaseSall(good,2));
    ci = bootci(nboot, @geomean, frBaseRall(good,2)./frBaseSall(good,2));

    fprintf("geometric mean baseline firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good))

    good = find(~(frStimSall(:,2)==0 | frStimSall(:,2)==0));
    good = good(~isinf(log(frStimRall(good,2)./frStimSall(good,2))));
    m = geomean(frStimRall(good,2)./frStimSall(good,2));
    ci = bootci(nboot, @geomean, frStimRall(good,2)./frStimSall(good,2));

    fprintf("geometric mean stim-driven firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)

    good = find(~(frPrefSall(:,2)==0 | frPrefRall(:,2)==0 | isnan(frPrefSall(:,2)) | isnan(frPrefRall(:,2))));
    good = good(~isinf(log(frPrefRall(good,2)./frPrefSall(good,2))));
    m = geomean(frPrefRall(good,2)./frPrefSall(good,2));
    ci = bootci(nboot, @geomean, frPrefRall(good,2)./frPrefSall(good,2));

    fprintf("geometric mean Pref-driven firing rate ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), NC)


%%
%pulling the following out into a new script 'asidecombine.m' for clarity
asidecombine_good


% ---v--- Skip this

%% Baseline firing
ntrials=numel(Stat.(subject).corrSacFRPvalB);
figure(95);clf
subplot(2,2,1);
hold off
histogram(Stat.(subject).corrSacFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade rate corr with baseline firing')
pc=nnz(Stat.(subject).corrSacFRPvalB>0.05)/numel(Stat.(subject).corrSacFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with saccade rate \n", pc, nnz(Stat.(subject).corrSacFRPvalB>0.05), ntrials)

%

subplot(2,2,2);
hold off
histogram(Stat.(subject).corrSacMagFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade mag corr with baseline firing')
pc=nnz(Stat.(subject).corrSacmagFRPvalB>0.05)/numel(Stat.(subject).corrSacmagFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with saccade magnitude \n", pc, nnz(Stat.(subject).corrSacmagFRPvalB>0.05), ntrials)


subplot(2,2,3); hold off
histogram(Stat.(subject).corrSizeFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Pupil size corr with baseline firing')
pc=nnz(Stat.(subject).corrSizeFRPvalB>0.05)/numel(Stat.(subject).corrSizeFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with pupil size \n", pc, nnz(Stat.(subject).corrSizeFRPvalB>0.05), ntrials)

%

subplot(2,2,4); hold off
histogram(Stat.(subject).corrVarXYFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Eyevar corr with baseline firing')
pc=nnz(Stat.(subject).corrVarXYFRPvalB>0.05)/numel(Stat.(subject).corrVarXYFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with eye position variance \n", pc, nnz(Stat.(subject).corrVarXYFRPvalB>0.05), ntrials)

%% Stimuli driven
figure(96);clf
subplot(2,2,1);
hold off
histogram(Stat.(subject).corrSacFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade rate corr with stim driven firing')
pc=nnz(Stat.(subject).corrSacFRPvalSt>0.05)/numel(Stat.(subject).corrSacFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with saccade rate \n", pc, nnz(Stat.(subject).corrSacFRPvalSt>0.05), ntrials)


subplot(2,2,2);
hold off
histogram(Stat.(subject).corrSacMagFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade mag corr with stim driven firing')
pc=nnz(Stat.(subject).corrSacmagFRPvalSt>0.05)/numel(Stat.(subject).corrSacmagFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with saccade magnitude \n", pc, nnz(Stat.(subject).corrSacmagFRPvalSt>0.05), ntrials)

subplot(2,2,3); hold off
histogram(Stat.(subject).corrSizeFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Pupil size corr with stim driven firing')
pc=nnz(Stat.(subject).corrSizeFRPvalSt>0.05)/numel(Stat.(subject).corrSizeFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with pupil size \n", pc, nnz(Stat.(subject).corrSizeFRPvalSt>0.05), ntrials)

%
subplot(2,2,4); hold off
histogram(Stat.(subject).corrVarXYFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Eyevar corr with stim driven firing')
pc=nnz(Stat.(subject).corrVarXYFRPvalSt>0.05)/numel(Stat.(subject).corrVarXYFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with eye position variance \n", pc, nnz(Stat.(subject).corrVarXYFRPvalSt>0.05), ntrials)


%%
fprintf("Running speed and firing rate: %.4f, r^2=%.4f \n\n",nanmean(Stat.(subject).runrho),nanmean(Stat.(subject).runrho.^2))

fprintf("Saccade rate and firing rate: %.4f, r^2=%.4f \n",(nanmean(Stat.(subject).corrSacFRRho)),nanmean(Stat.(subject).corrSacFRRho.^2))
fprintf("Saccade rate and base firing rate: %.4f, r^2=%.4f \n",(nanmean(Stat.(subject).corrSacFRRhoB)),nanmean(Stat.(subject).corrSacFRRhoB.^2))
fprintf("Saccade rate and stim-driven firing rate: %.4f, r^2=%.4f \n\n",(nanmean(Stat.(subject).corrSacFRRhoSt)),nanmean(Stat.(subject).corrSacFRRhoSt.^2))

fprintf("Pupil size and firing rate: %.4f, r^2=%.4f \n",(nanmean(Stat.(subject).corrSizeFRRho)),(nanmean(Stat.(subject).corrSizeFRRho))^2)
fprintf("Pupil size and base firing rate: %.4f, r^2=%.4f \n",(nanmean(Stat.(subject).corrSizeFRRhoB)),(nanmean(Stat.(subject).corrSizeFRRhoB))^2)
fprintf("Pupil size and stim-driven firing rate: %.4f, r^2=%.4f \n\n",(nanmean(Stat.(subject).corrSizeFRRhoSt)),nanmean(Stat.(subject).corrSizeFRRhoSt.^2))

fprintf("Eyevar and firing rate: %.4f, r^2=%.4f \n",(nanmean(Stat.(subject).corrSacMagFRRho)),(nanmean(Stat.(subject).corrSacMagFRRho))^2)
fprintf("Eyevar and base firing rate: %.4f, r^2=%.4f \n",(nanmean(Stat.(subject).corrSacMagFRRhoB)),(nanmean(Stat.(subject).corrSacMagFRRhoB))^2)
fprintf("Eyevar and stim-driven firing rate: %.4f, r^2=%.4f \n\n",(nanmean(Stat.(subject).corrSacMagFRRhoSt)),nanmean(Stat.(subject).corrSacMagFRRhoSt.^2))



% fprintf("",(nanmean(Stat.(subject).corrSacMagFRRhoSt))
% fprintf("",(nanmean(Stat.(subject).corrSacMagFRRhoSt))^2

% Careful about outliers
baseR=Stat.(subject).frBaseR(:,2);
baseR((baseR)>nanmean(baseR)+3*nanstd(baseR))=nan;
baseR((baseR)<nanmean(baseR)-3*nanstd(baseR))=nan;
baseS=Stat.(subject).frBaseS(:,2);
baseS((baseS)>nanmean(baseS)+3*nanstd(baseS))=nan;
baseS((baseS)<nanmean(baseS)-3*nanstd(baseS))=nan;
baseratio=(baseR./baseS);
baseratio=baseratio(~isnan(baseratio));
baseratio=baseratio(~isnan(baseratio));
baseratio=baseratio(baseratio~=0);

Basegm=geomean(baseratio,'omitnan' );
fprintf("Geomean (running/stationary) during baseline: %.4f\n" ,Basegm)
stimratio=Stat.(subject).frStimR(:,2)./Stat.(subject).frStimS(:,2);
stimratio=stimratio(stimratio~=0);
Stimgm=geomean(stimratio,'omitnan' );
fprintf("Geomean (running/stationary) during stimuli: %.4f\n" ,Stimgm)
prefratio=Stat.(subject).murun(:,1)./(Stat.(subject).mustat(:,1));
prefratio=prefratio(prefratio~=0);
prefdirgm=geomean(prefratio,'omitnan' );
fprintf("Geomean at preferred direction (running/stationary): %.4f\n" ,prefdirgm)

%%

SacHzDelta=(Stat.(subject).SacHzR-Stat.(subject).SacHzS); %Per cell diff in saccades
ExpectedFRDelta=SacHzDelta.*Stat.(subject).regRnSac_B(:,1); %Per cell expected change in firing rate

FRDeltaStim=Stat.(subject).frStimR(:,2)-Stat.(subject).frStimS(:,2);
fprintf("Mean delta (running-stationary) observed during stimuli: %.4f Hz \n" ,nanmean(FRDeltaStim))
fprintf("Expected mean delta(running-stationary) from reg. on nSaccades: %.4f Hz \n" ,nanmean(ExpectedFRDelta(:,2)))
%fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-ExpectedFRDelta(:,2)))

clear meanspd
for ii=1:size(Stat.(subject).runningspeed)
if ~isempty(Stat.(subject).runningspeed{ii})
meanspd(ii)=nanmean(Stat.(subject).runningspeed{ii});
else
meanspd(ii)=nan;
end
end

%Mean amount of running on trial * B for running effect
ExpectedFRRunDelta=meanspd'.*Stat.(subject).regRnSpd_B(:,1);

fprintf("Expected mean delta(running-stationary) from reg. on run speed: %.4f Hz \n" ,nanmean(ExpectedFRRunDelta))
%fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-ExpectedFRRunDelta))


%fprintf("Expected running rate from geomean: mu*statrate %.4f Hz \n" ,nanmean(Stimgm*Stat.(subject).frStimS(:,2)))
meandelta=(Stat.(subject).frStimS(:,2).*nanmean(stimratio) - Stat.(subject).frStimS(:,2));
fprintf("Expected mean delta(running-stationary) from mean running effect: %.4f  Hz\n", nanmean(meandelta))
%fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-meandelta))
fprintf("(Mean only uses one number for all cells) \n \n")
% FRDeltaStim./Stat.(subject).frStimS;
%  nanmean(FRDeltaStim./Stat.(subject).frStimS(:,2))

%%

% FRDeltaStim=Stat.(subject).frStimR(:,2)-Stat.(subject).frStimS(:,2);
% meangain=mean(FRDeltaStim./Stat.(subject).frStimS(:,2),'omitnan');
% fprintf("Mean gain (running-stationary)/stationary observed during stimuli: %.4f  \n" ,meangain)
% 
% SacHzDelta=(Stat.(subject).SacHzR-Stat.(subject).SacHzS); %Per cell diff in saccades
% ExpectedFRDelta=SacHzDelta.*Stat.(subject).regRnSac_B(:,1); %Per cell expected change in firing rate
% meangainsac=mean(ExpectedFRDelta(:,2)./Stat.(subject).frStimS(:,2),'omitnan');
% fprintf("Expected gain (running-stationary)/stationary from reg. on nSaccades: %.4f  \n" ,meangainsac)
% 
% meangainrun=mean((ExpectedFRRunDelta)./Stat.(subject).frStimS(:,2),'omitnan');
% fprintf("Expected gain (running-stationary)/stationary from reg. on run speed: %.4f  \n" ,meangainrun)
% 
% 
% meandelta=(Stat.(subject).frStimS(:,2).*nanmean(stimratio) - Stat.(subject).frStimS(:,2));
% meanmeangain=mean((meandelta)./Stat.(subject).frStimS(:,2),'omitnan');
% fprintf("Expected gain (running-stationary)/stationary from mean running effect: %.4f  \n", meanmeangain)
% %fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-meandelta))
% fprintf("(Mean only uses one number for all cells) \n \n")
% % FRDeltaStim./Stat.(subject).frStimS;
%% Delta redon
FRDeltaStim=Stat.(subject).frStimR(:,2)-Stat.(subject).frStimS(:,2);
meangain=mean(Stat.(subject).frStimR(:,2)./Stat.(subject).frStimS(:,2),'omitnan');
fprintf("Mean delta (running-stationary) observed during stimuli: %.4f Hz\n" ,nanmean(FRDeltaStim))

SacHzDelta=(Stat.(subject).SacHzR-Stat.(subject).SacHzS); %Per cell diff in saccades
ExpectedFRDelta=SacHzDelta.*Stat.(subject).regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRSacHzR=Stat.(subject).SacHzR.*Stat.(subject).regRnSac_B(:,1)+Stat.(subject).regRnSac_B(:,2);
ExpectedFRSacHzS=Stat.(subject).SacHzS.*Stat.(subject).regRnSac_B(:,1)+Stat.(subject).regRnSac_B(:,2);
ExpectedFRDelta2=ExpectedFRSacHzR-ExpectedFRSacHzS;
fprintf("Expected delta (running-stationary) from reg. on nSaccades: %.4f Hz \n" ,nanmean(ExpectedFRDelta2(:,2)))


speedDelta=(Stat.(subject).speedR-Stat.(subject).speedS); %Per cell diff in saccades
ExpectedFRspeedDelta=speedDelta.*Stat.(subject).regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRspeedR=Stat.(subject).speedR.*Stat.(subject).regRnSpd_B(:,1)+Stat.(subject).regRnSpd_B(:,2);
ExpectedFRspeedS=Stat.(subject).speedS.*Stat.(subject).regRnSpd_B(:,1)+Stat.(subject).regRnSpd_B(:,2);
meangainsac=mean(ExpectedFRspeedR(:,2)./ExpectedFRspeedS(:,2),'omitnan');
fprintf("Expected delta (running-stationary) from reg. on run speed: %.4f Hz \n" ,nanmean(ExpectedFRspeedDelta(:,2)))

%
FRDeltaStim=Stat.(subject).frStimR(:,2)-Stat.(subject).frStimS(:,2);
meangain=mean(Stat.(subject).frStimR(:,2)./Stat.(subject).frStimS(:,2),'omitnan');
fprintf("Mean gain (running/stationary) observed during stimuli: %.4f  \n" ,meangain)

SacHzDelta=(Stat.(subject).SacHzR-Stat.(subject).SacHzS); %Per cell diff in saccades
ExpectedFRDelta=SacHzDelta.*Stat.(subject).regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRSacHzR=Stat.(subject).SacHzR.*Stat.(subject).regRnSac_B(:,1)+Stat.(subject).regRnSac_B(:,2);
ExpectedFRSacHzS=Stat.(subject).SacHzS.*Stat.(subject).regRnSac_B(:,1)+Stat.(subject).regRnSac_B(:,2);
meangainsac=mean(ExpectedFRSacHzR(:,2)./ExpectedFRSacHzS(:,2),'omitnan');
fprintf("Expected gain (running/stationary) from reg. on nSaccades: %.4f  \n" ,meangainsac)


speedDelta=(Stat.(subject).speedR-Stat.(subject).speedS); %Per cell diff in saccades
ExpectedFRDelta=speedDelta.*Stat.(subject).regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRspeedR=Stat.(subject).speedR.*Stat.(subject).regRnSpd_B(:,1)+Stat.(subject).regRnSpd_B(:,2);
ExpectedFRspeedS=Stat.(subject).speedS.*Stat.(subject).regRnSpd_B(:,1)+Stat.(subject).regRnSpd_B(:,2);
meangainrun=mean(ExpectedFRspeedR(:,2)./ExpectedFRspeedS(:,2),'omitnan');
fprintf("Expected gain (running/stationary) from reg. on run speed: %.4f  \n" ,meangainrun)



%%

subjects = {'gru', 'brie'};
nsubjs = numel(subjects);
field = 'rateweight';
for isubj = 1:nsubjs
    subject = subjects{isubj};
    cmap = getcolormap(subject, false);
    fprintf('***************************\n\n')
    fprintf('***************************\n\n')
    fprintf('%s\n', subject)

    NC = numel(Stat.(subject).meanrate);
    good = find(~isnan(Stat.(subject).meanrate));
    dsi = nan(NC,1);
    osi = nan(NC,2);
    maxrate = nan(NC,1);

    for cc = good(:)'
        mu = squeeze(Stat.(subject).(['d' field])(cc,:,2));
        maxrate(cc) = max(mu);

        baseline = min(Stat.(subject).baselinerate(cc), 1);
%         baseline = 0; %min(min(mu(:)), baseline);
        dsi(cc) = direction_selectivity_index(Stat.(subject).directions, mu(:)-baseline, 2);
        mu = squeeze(Stat.(subject).(['o' field])(cc,:,2));
%         baseline = min(mu(:));
        osi(cc,1) = orientation_selectivity_index(Stat.(subject).orientations, mu(:)-baseline, 2);
        mu = squeeze(Stat.(subject).(['o' field])(cc,:,1));
        osi(cc,2) = orientation_selectivity_index(Stat.(subject).orientations, mu(:)-baseline, 2);
    end

    osi(osi(:,1) < 0 | osi(:,1) > 1,:) = nan;

    figure(isubj); clf

%     histogram(dsi, 'binEdges', linspace(0,1,100))
%     hold on
%     histogram(osi, 'binEdges', linspace(0,1,100), 'FaceColor', cmap(6,:), 'EdgeColor', 'none', 'FaceAlpha', .5)
%     hold on
    iix = Stat.(subject).pvis(:,3)<0.05;
    histogram(osi(iix), 'binEdges', linspace(0,1,100), 'FaceColor', cmap(6,:), 'EdgeColor', 'none', 'FaceAlpha', 1)
    xlabel('OSI')
    set(gca, 'XTick', 0:.25:1)
    plot.formatFig(gcf, [1.2 1], 'nature')
    saveas(gcf, fullfile(figdir, sprintf('osi_dist_%s.pdf', subject)))


    mean(Stat.(subject).pvis(good,3)<0.05)

    osi(maxrate < 5 | maxrate > 50,:) = nan;

    figure(isubj + 10); clf
    excludefromplotting = ~((osi(:,2)./osi(:,1)) > .8 & (osi(:,2)./osi(:,1)) < 1.5);
    osi(excludefromplotting,:) = nan;
    metric = osi(:,2);

    [~, inds] = nanmin((metric(:) - [.1 .5 .9]).^2);
    

    for i = 1:3
        cc = inds(i);
        subplot(1,3,i)
        mu = squeeze(Stat.(subject).(['d' field])(cc,:,:));
        plot(Stat.(subject).directions, mu, 'Color', cmap(6,:))
        hold on
        plot(Stat.(subject).directions, mu(:,2), '.', 'Color', cmap(6,:))

        plot(xlim, Stat.(subject).baselinerate(cc)*[1 1], '--', 'Color', cmap(1,:))
        title(sprintf('OSI: %02.2f', metric(cc)))
        
        plot.offsetAxes(gca, false, 10)
        ylim([0 max(ylim)])
        set(gca, 'XTick', 0:45:330)
        xlim(Stat.(subject).directions([1 end]))
        xlabel('Orienation (deg.)')
        if i==1
            ylabel('Firing Rate (sp s^{-1})')
        end
    end
    plot.formatFig(gcf, [3 1], 'nature')
    saveas(gcf, fullfile(figdir, sprintf('example_tuning_%s.pdf', subject)))

end

%%
subjects = {'gru', 'brie'};
nsubjs = numel(subjects);
field = 'rateweight';
for isubj = 1:nsubjs
    subject = subjects{isubj};
    cmap = getcolormap(subject, false);
    fprintf('***************************\n\n')
    fprintf('***************************\n\n')
    fprintf('%s\n', subject)


end
%%
mean(corrPval < 0.05)

figure(1); clf
histogram(corrRho, -.5:.025:.5); hold on
histogram(corrRho(corrPval < 0.05), -.5:.025:.5)
xlabel("Spearman's Rho")
ylabel('Count')
legend({'All', 'p<0.05'})

%%
plot(sum(robs,2))


%%

goodix = getStableRange(R, 'plot', false); % key function: finds stable region of firing

%% Do direction / orientation decoding

% example session to explore
Dstat = decode_stim(D, 12);


%%
Nsess = max(D.sessNumGratings);
clear Dstat
for sessionId = 1:Nsess
    Dstat(sessionId) = decode_stim(D, sessionId, 'runThreshold', 3);
end

%%

%% decoding error
circdiff = @(x,y) angle(exp(1i*(x - y)/180*pi))/pi*180;

mR = zeros(Nsess,1);
mRCi = zeros(Nsess,2);
mS = zeros(Nsess,1);
mSCi = zeros(Nsess,2);
nS = zeros(Nsess,1);
nR = zeros(Nsess,1);

figure(1); clf
for iSess = 1:Nsess
    aerr = abs(circdiff(Dstat(iSess).Stim, Dstat(iSess).decoderStimTot));
    
    inds = Dstat(iSess).runTrials;
    nR(iSess) = numel(inds);
    mR(iSess) = median(aerr(inds));
    mRCi(iSess,:) = bootci(1000, @median, aerr(inds));
    
    inds = setdiff(1:Dstat(iSess).NTrials, Dstat(iSess).runTrials);
    nS(iSess) = numel(inds);
    mS(iSess) = median(aerr(inds));
    mSCi(iSess,:) = bootci(1000, @median, aerr(inds));
    
    plot(mS(iSess)*[1 1], mRCi(iSess,:), 'Color', .5*[1 1 1]); hold on
    plot(mSCi(iSess,:), mR(iSess)*[1 1], 'Color', .5*[1 1 1]);
    h = plot(mS(iSess), mR(iSess), 'o');
    h.MarkerFaceColor = h.Color;
    
end

plot(xlim, xlim, 'k')
xlim([0 20])
ylim([0 20])
xlabel('Stationary')
ylabel('Running')
title('Median Decoding Error (degrees)')

%%
sessionId=12;
clear Dstat
Dstat(sessionId) = decode_running(D, sessionId, 'Decoder', 'svm', 'runThreshold', 3);
%% Decode running
Nsess = max(D.sessNumGratings);
clear Dstat
for sessionId = 1:Nsess
    Dstat(sessionId) = decode_running(D, sessionId, 'Decoder', 'svm', 'runThreshold', 3);
end
%%

%%


figure(1); clf
for s = 1:Nsess
    h = plot(Dstat(s).chance*[1 1], Dstat(s).accCi);
    hold on
    plot(Dstat(s).chance, Dstat(s).acc, 'o', 'Color', h.Color, 'MarkerFaceColor', h.Color);
end
plot(xlim, xlim, 'k')
xlabel('Chance level (Based on % running)')
ylabel('Decoding Accuracy')

%% Bootstrapped empirical analyses and tuning curve fits

fpath = getpref('FREEVIEWING', 'HUKLAB_DATASHARE');
fname = fullfile(fpath, [subject 'TCstats.mat']);

if exist(fname, 'file')
    stat = load(fname);
else
    stat = tuning_empirical(D, 'binsize', 10e-3, ...
        'runningthresh', 3, ...
        'nboot', 500, ...
        'seed', 1234);
    save(fname, '-v7.3', '-struct', 'stat')
end

%% combine two subjects?
subjects = {'gru', 'brie'};
stat = [];
for i = 1:2
    subject = subjects{i};
    fname = fullfile(fpath, [subject 'TCstats.mat']);
    s = load(fname);
    
    if isempty(stat)
        stat = s;
        continue
    end
    
    nUnits = numel(stat.TCempirical.TCdiff);
    nNew = numel(s.TCempirical.TCdiff);
    
    snew = struct();

    % ---- TC EMPIRICAL
    fields = fieldnames(stat.TCempirical);
    thetas = unique([stat.TCempirical.thetas; s.TCempirical.thetas]);
    nthetas = numel(thetas);
    
    for ifield = 1:numel(fields)
        sz = size(stat.TCempirical.(fields{ifield}));
        sznew = size(s.TCempirical.(fields{ifield}));
        unitDim = find(sz==nUnits);
        if ~isempty(unitDim)
            nonUnitDims = setdiff(1:numel(sz), unitDim);
            if sz(1)==numel(stat.TCempirical.thetas)
                snew.TCempirical.(fields{ifield}) = zeros(nthetas, nUnits+nNew, 3);
                snew.TCempirical.(fields{ifield})(ismember(stat.TCempirical.thetas, thetas), 1:nUnits,:) = stat.TCempirical.(fields{ifield});
                snew.TCempirical.(fields{ifield})(ismember(s.TCempirical.thetas, thetas), nUnits+(1:nNew),:) = s.TCempirical.(fields{ifield});
            else
                snew.TCempirical.(fields{ifield}) = zeros(nUnits+nNew, sz(nonUnitDims));
                snew.TCempirical.(fields{ifield})(1:nUnits,:) = stat.TCempirical.(fields{ifield});
                snew.TCempirical.(fields{ifield})(nUnits + (1:nNew),:) = s.TCempirical.(fields{ifield});
            end
        end
    end

    snew.TCempirical.thetas = thetas;
    
    % --- struct arrays
    snew.running = [stat.running; s.running];
    snew.TCfitR = [stat.TCfitR s.TCfitR];
    snew.TCfitS = [stat.TCfitS s.TCfitS];
    
    % --- SPEED TUNING
    fields = fieldnames(stat.speedTuning);
    
    for ifield = 1:numel(fields)
        sz = size(stat.speedTuning.(fields{ifield}));
        sznew = size(s.speedTuning.(fields{ifield}));
        unitDim = find(sz==nUnits);
        if ~isempty(unitDim)
            nonUnitDims = setdiff(1:numel(sz), unitDim);
            
            snew.speedTuning.(fields{ifield}) = zeros(sz(nonUnitDims), nUnits+nNew);
            snew.speedTuning.(fields{ifield})(:,1:nUnits) = stat.speedTuning.(fields{ifield});
            snew.speedTuning.(fields{ifield})(:, nUnits + (1:nNew)) = s.speedTuning.(fields{ifield});
        else
            snew.speedTuning.(fields{ifield}) = s.speedTuning.(fields{ifield});
        end
    end
    
    % --- PSTHS
    fields = fieldnames(stat.psths);
    
    for ifield = 1:numel(fields)
        sz = size(stat.psths.(fields{ifield}));
        sznew = size(s.psths.(fields{ifield}));
        unitDim = find(sz==nUnits);
        if ~isempty(unitDim)
            nonUnitDims = setdiff(1:numel(sz), unitDim);
            newd = arrayfun(@(x) x, sz(nonUnitDims), 'uni', 0);
            snew.psths.(fields{ifield}) = zeros(newd{:}, nUnits+nNew);
            snew.psths.(fields{ifield})(:,ismember(stat.TCempirical.thetas, thetas),:,1:nUnits) = stat.psths.(fields{ifield});
            snew.psths.(fields{ifield})(:,ismember(s.TCempirical.thetas, thetas),:, nUnits + (1:nNew)) = s.psths.(fields{ifield});
        else
            snew.psths.(fields{ifield}) = s.psths.(fields{ifield});
        end
    end
    
    
end
    
stat = snew;
%% plot running tuning?

NC = numel(stat.TCfitR);
runrat = stat.speedTuning.rateSpdMu./stat.speedTuning.rateSpdMu(1,:);
[~, ind] = sort(runrat(end,:));


sx = ceil(sqrt(NC));
sy = round(sqrt(NC));


figure(1); clf
ax = plot.tight_subplot(sx, sy, 0.001, 0.001);
cmap = lines;
cmap(1,:) = .2*[1 1 1];
for i = 1:(sx*sy)
    
    
    set(gcf, 'currentaxes', ax(i))
    
    if i > NC
        axis off
        continue
    end
    
    fprintf('Unit %d/%d\n', i, NC)
    
    cc = ind(i);
    
    plot.errorbarFill(stat.speedTuning.bins, stat.speedTuning.rateSpdMu(:,cc), stat.speedTuning.rateSpdSe(:,cc), 'b', 'FaceColor', cmap(1,:), 'EdgeColor', 'none', 'FaceAlpha', .5); hold on
    plot(stat.speedTuning.bins, stat.speedTuning.rateSpdMu(:,cc), 'o', 'Color', cmap(1,:))
    
%     plot.errorbarFill(stat.speedTuning.bins, stat.speedTuning.rateSpdMuStim(:,cc), stat.speedTuning.rateSpdSeStim(:,cc), 'r')
    plot(xlim, (stat.speedTuning.rateSpdMu(1,cc) + stat.speedTuning.rateSpdSe(1,cc))*[1 1], 'k--')
    plot(xlim, (stat.speedTuning.rateSpdMu(1,cc) - stat.speedTuning.rateSpdSe(1,cc))*[1 1], 'k--')
%     plot(xlim, stat.speedTuning.rateSpdMu(1,cc)*[1 1], 'k')
    xlabel('Speed (cm / s)')
    ylabel('Firing Rate')
    axis off
    
    drawnow
    
end




%% plot all running speed tuning curves on top of each other

figure(10); clf
runrat = stat.speedTuning.rateSpdMu - mean(stat.speedTuning.rateSpdMu(1:3,:));

plot(stat.speedTuning.bins, runrat, '-', 'MarkerSize', 2, 'Color', [cmap(1,:) .25] ); hold on
xlabel('Running Speed (cm/s)')
ylabel('\Delta Firing Rate')
plot(stat.speedTuning.bins, nanmean(runrat, 2), 'r', 'Linewidth', 2)
plot(xlim, [0 0], 'b--', 'Linewidth', 2)
ylim([-5 5])


%% Raw Running Modulation (ignore stimulus entirely)
% Just look at the entire session. Using labeled epochs of running and
% stationary, count the mean firing rate (while accounting for issues with
% stationarity by resampling from the epochs to hopefully match)

nExamples = 5;
figure(1); clf
set(gcf, 'Color', 'w')
NC = numel(stat.running);
rateS = nan(NC, 3);
rateR = nan(NC, 1);
cmap = lines;

nbins = numel(stat.running(1).psthMu);
psthRunOnset = nan(nbins, NC);
psthBins = stat.running(1).psthBins;
numEpochs = nan(NC, 1);
isvalid = find(arrayfun(@(x) ~isempty(x.spikerate), stat.running));

subplot(2,2,1)
for cc = 1:NC
    if isempty(stat.running(cc).spikerate)
        continue
    end
    
    rateS(cc,:) = prctile(stat.running(cc).rateStatNull, [2.5 50 97.5]);
    rateR(cc) = stat.running(cc).rateRun;
    
    plot(rateS(cc,[1 3]),rateR(cc)*[1 1], 'Color', .5*[1 1 1]); hold on
    plot(rateS(cc,2), rateR(cc), 'ow', 'MarkerFaceColor', cmap(1,:))
   
    psthRunOnset(:,cc) = stat.running(cc).psthMu;
    
    numEpochs(cc) = numel(stat.running(cc).goodix);
end

plot(xlim, xlim, 'k')
title('Mean Firing Rate', 'Fontweight', 'normal')
xlabel('Stationary')
ylabel('Running')


suppressed = find(rateR < rateS(:,1));
enhanced = find(rateR > rateS(:,3));
nEnc = numel(enhanced);
nSup = numel(suppressed);
nTot = sum(~isnan(rateR));

fprintf('Found %d/%d enhanced (%02.2f%%)\n', nEnc, nTot, 100*nEnc/nTot)
fprintf('Found %d/%d suppressed (%02.2f%%)\n', nSup, nTot, 100*nSup/nTot)

for cc = suppressed(:)'
    plot(rateS(cc,[1 3]),rateR(cc)*[1 1], 'Color', .5*[1 1 1]); hold on
    plot(rateS(cc,2), rateR(cc), 'ow', 'MarkerFaceColor', cmap(2,:))
end

for cc = enhanced(:)'
    plot(rateS(cc,[1 3]),rateR(cc)*[1 1], 'Color', .5*[1 1 1]); hold on
    plot(rateS(cc,2), rateR(cc), 'ow', 'MarkerFaceColor', cmap(4,:))
end

subplot(2,2,2)
nfun = @(x) x./mean(x(psthBins<0,:));
plot(psthBins, nfun(psthRunOnset(:,enhanced)), 'Color', (1+cmap(4,:))/2); hold on
plot(psthBins, nfun(psthRunOnset(:,suppressed)), 'Color', (1+cmap(2,:))/2); hold on

plot(psthBins, mean(nfun(psthRunOnset(:,enhanced)),2), 'Color', cmap(4,:), 'Linewidth', 2)
plot(psthBins, mean(nfun(psthRunOnset(:,suppressed)),2), 'Color', cmap(2,:), 'Linewidth', 2)
xlabel('Time from running onset (s)')
ylabel('Relative Rate (mean normalized)')
ylim([.5 2])
xlim(psthBins([1 end]))

subplot(2,2,3)
nfun = @(x) (x - mean(x(psthBins<0,:)))./std(x(psthBins<0,:));
plot(psthBins, nfun(psthRunOnset(:,enhanced)), 'Color', (1+cmap(4,:))/2); hold on
plot(psthBins, nfun(psthRunOnset(:,suppressed)), 'Color', (1+cmap(2,:))/2); hold on

plot(psthBins, mean(nfun(psthRunOnset(:,enhanced)),2), 'Color', cmap(4,:), 'Linewidth', 2)
plot(psthBins, mean(nfun(psthRunOnset(:,suppressed)),2), 'Color', cmap(2,:), 'Linewidth', 2)
xlabel('Time from running onset (s)')
ylabel('Normalized Rate (z score)')
ylim([-5 5])
xlim(psthBins([1 end]))

subplot(2,2,4)
nfun = @(x) (x - mean(x(psthBins<0,:)));
plot(psthBins, nfun(psthRunOnset(:,enhanced)), 'Color', (1+cmap(4,:))/2); hold on
plot(psthBins, nfun(psthRunOnset(:,suppressed)), 'Color', (1+cmap(2,:))/2); hold on

plot(psthBins, mean(nfun(psthRunOnset(:,enhanced)),2), 'Color', cmap(4,:), 'Linewidth', 2)
plot(psthBins, mean(nfun(psthRunOnset(:,suppressed)),2), 'Color', cmap(2,:), 'Linewidth', 2)
xlabel('Time from running onset (s)')
ylabel('\Delta Rate (spikes/sec)')
ylim([-5 5])
xlim(psthBins([1 end]))

% Sanity Check: Check that this effect isn't a function of the number of running epochs
figure(2); clf
set(gcf, 'Color', 'w')
plot([1; 1]*numEpochs', (rateR-rateS(:,[1 3]))', '-k', 'Linewidth', 2); hold on
plot([1; 1]*numEpochs(suppressed)', (rateR(suppressed)-rateS(suppressed,[1 3]))', '-', 'Color', cmap(2,:), 'Linewidth', 2);
plot([1; 1]*numEpochs(enhanced)', (rateR(enhanced)-rateS(enhanced,[1 3]))', '-', 'Color', cmap(4,:), 'Linewidth', 2);
plot(xlim, [0 0], 'k--')
xlabel('Num Running Epochs')
ylabel('\Delta Firing Rate')

% display top 5 examples of suppression and enhancement
deltaFR = rateR(isvalid) - rateS(isvalid,2);
[~, ind] = sort(deltaFR);
ind = isvalid(ind);

figure(3); clf
set(gcf, 'Color', 'w')

spacing = 0.05;
ax = plot.tight_subplot(nExamples, 2, spacing, 0.05, 0.05);
sm = 20;


for i = 1:nExamples
    cc = ind(i);
    
    set(gcf, 'currentaxes', ax((i-1)*2+1))
    yyaxis left
    plot(imgaussfilt(stat.running(cc).spikerate, sm))
    ylabel('Firing Rate')
    axis tight
    
    yyaxis right
    plot(imgaussfilt(stat.running(cc).runningspeed, sm))
    ylabel('Running Speed')
    axis tight
    if i==1
        title('Most Suppressed')
    end
end
xlabel('Time')

for i = 1:nExamples
    cc = ind(end-(i-1));

    set(gcf, 'currentaxes', ax((i-1)*2+2))
    yyaxis left
    plot(imgaussfilt(stat.running(cc).spikerate, sm))
    ylabel('Firing Rate')
    axis tight
    
    yyaxis right
    plot(imgaussfilt(stat.running(cc).runningspeed, sm))
    ylabel('Running Speed')
    axis tight
    
    if i==1
        title('Most Enhanced')
    end
end

xlabel('Time')

% PLOT running onset-aligned PSTHs 
figure(4); clf
set(gcf, 'Color', 'w')
ax = plot.tight_subplot(nExamples, 2, spacing, 0.05, 0.05);

for i = 1:nExamples
    cc = ind(i);
    
    set(gcf, 'currentaxes', ax((i-1)*2+1))
    plot.errorbarFill(stat.running(cc).psthBins, stat.running(cc).psthMu, stat.running(cc).psthSe); hold on
    plot(stat.running(cc).psthBins, stat.running(cc).psthMu, 'k', 'Linewidth', 2)
    plot(stat.running(cc).psthBins, stat.running(cc).psthNullCi, 'r--')
    ylabel('Firing Rate')
    axis tight
    
    if i==1
        title('Most Suppressed')
    end
end
xlabel('Time from Running Onset (s)')

for i = 1:nExamples
    cc = ind(end-(i-1));

    set(gcf, 'currentaxes', ax((i-1)*2+2))
    plot.errorbarFill(stat.running(cc).psthBins, stat.running(cc).psthMu, stat.running(cc).psthSe); hold on
    plot(stat.running(cc).psthBins, stat.running(cc).psthMu, 'k', 'Linewidth', 2)
    plot(stat.running(cc).psthBins, stat.running(cc).psthNullCi, 'r--')
    axis tight
    
    if i==1
        title('Most Enhanced')
    end
end
xlabel('Time from Running Onset (s)')

% TUning curve analysis
figure(5); clf
set(gcf, 'Color', 'w')
ax = plot.tight_subplot(nExamples, 2, spacing, 0.05, 0.05);

for i = 1:nExamples
    cc = ind(i);
    
    set(gcf, 'currentaxes', ax((i-1)*2+1))
    plot.errorbarFill(stat.TCfitS(cc).thetas, stat.TCfitS(cc).tuningCurve, stat.TCfitS(cc).tuningCurveSE*2, 'k', 'FaceColor', cmap(1,:), 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
    plot(stat.TCfitS(cc).thetas, stat.TCfitS(cc).tuningCurve, 'o', 'Color', cmap(1,:))
    plot.errorbarFill(stat.TCfitR(cc).thetas, stat.TCfitR(cc).tuningCurve, stat.TCfitR(cc).tuningCurveSE*2, 'k', 'FaceColor', cmap(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
    plot(stat.TCfitR(cc).thetas, stat.TCfitR(cc).tuningCurve, 'o', 'Color', cmap(2,:))
    
    ylabel('Firing Rate')
    axis tight
    
    if i==1
        title('Most Suppressed')
    end
end
xlabel('Direction')

for i = 1:nExamples
    cc = ind(end-(i-1));

    set(gcf, 'currentaxes', ax((i-1)*2+2))
    plot.errorbarFill(stat.TCfitS(cc).thetas, stat.TCfitS(cc).tuningCurve, stat.TCfitS(cc).tuningCurveSE, 'k', 'FaceColor', cmap(1,:), 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
    plot(stat.TCfitS(cc).thetas, stat.TCfitS(cc).tuningCurve, 'o', 'Color', cmap(1,:))
    plot.errorbarFill(stat.TCfitR(cc).thetas, stat.TCfitR(cc).tuningCurve, stat.TCfitR(cc).tuningCurveSE, 'k', 'FaceColor', cmap(2,:), 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
    plot(stat.TCfitR(cc).thetas, stat.TCfitR(cc).tuningCurve, 'o', 'Color', cmap(2,:))
    ylabel('Firing Rate')
    axis tight
    
    if i==1
        title('Most Enhanced')
    end
end
xlabel('Time from Running Onset (s)')



%% TC diff

nullci = prctile(stat.TCempirical.maxFRdiffNull, [2.5 97.5], 2);

nInc = sum(stat.TCempirical.maxFRdiff > nullci(:,2));
nDec = sum(stat.TCempirical.maxFRdiff < nullci(:,1));
fprintf('%d/%d (%02.2f%%) units had significantly increased firing (outside null) rate while running\n', nInc, NC, 100*nInc/NC)
fprintf('%d/%d (%02.2f%%) units had significantly decreased firing (outside null) rate while running\n', nDec, NC, 100*nDec/NC)
fprintf('%d/%d (%02.2f%%) units had any significant modulation of max firing rate\n', nDec+nInc, NC, 100*(nDec+nInc)/NC)



maxR = [];
maxS = [];
for i = 1:3
    maxR = [maxR max(stat.TCempirical.TCrun(:,:,i))'];
    maxS = [maxS max(stat.TCempirical.TCstat(:,:,i))'];
end

figure(1); clf
cmap = lines;

plot(maxS(:,[1 3])', maxR(:,[2 2])', 'Color', .5*[1 1 1]); hold on
plot(maxS(:,[2 2])', maxR(:,[1 3])', 'Color', .5*[1 1 1])
plot(maxS(:,2), maxR(:,2), 'o', 'Color', cmap(1,:), 'MarkerFaceColor', cmap(1,:), 'MarkerSize', 2);
plot(xlim, xlim, 'k')
xlabel('Stationary')
ylabel('Running')
title('Max Firing Rate')


nInc = sum(maxR(:,2) > maxS(:,3));
nDec = sum(maxR(:,2) < maxS(:,1));
fprintf('%d/%d (%02.2f%%) units had significantly increased firing rate while running\n', nInc, NC, 100*nInc/NC)
fprintf('%d/%d (%02.2f%%) units had significantly decreased firing rate while running\n', nDec, NC, 100*nDec/NC)


% --- MEDIAN --- %
% DIFFERENCE
maxDiff = maxR(:,2)-maxS(:,2);
m = nanmedian(maxDiff);
ci = bootci(1000, @nanmedian, maxDiff);
fprintf('Median difference (Running-Stationary) = %02.2f [%02.2f, %02.2f] (n=%d)\n', m, ci(1), ci(2), sum(~isnan(maxDiff)))

% RATIO
maxRat = maxR(:,2)./maxS(:,2);
good = ~(isnan(maxRat) | isinf(maxRat));
m = nanmedian(maxRat(good));
ci = bootci(1000, @nanmedian, maxRat(good));
fprintf('Median ratio (Running:Stationary) = %02.3f [%02.3f, %02.3f] (n=%d)\n', m, ci(1), ci(2), sum(good))

% --- MEAN --- %
% DIFFERENCE
m = nanmean(maxDiff);
ci = bootci(1000, @nanmean, maxDiff);
fprintf('Mean difference (Running-Stationary) = %02.2f [%02.2f, %02.2f] (n=%d)\n', m, ci(1), ci(2), sum(~isnan(maxDiff)))

% RATIO
good = ~(isnan(maxRat) | isinf(maxRat));
m = nanmean(maxRat(good));
ci = bootci(1000, @nanmean, maxRat(good));
fprintf('Mean ratio (Running:Stationary) = %02.3f [%02.3f, %02.3f] (n=%d)\n', m, ci(1), ci(2), sum(good))

%% examples
[~, ind] = sort(maxRat);

cc = cc + 1;


tcr = squeeze(stat.TCempirical.TCrun(:,ind(cc),:));
tcs = squeeze(stat.TCempirical.TCstat(:,ind(cc),:));

figure(1); clf
subplot(1,2,1)

plot(tcr, 'r'); hold on
plot(tcs, 'b')
plot(tcr(:,2), 'r', 'Linewidth', 2); hold on
plot(tcs(:,2), 'b', 'Linewidth', 2)

title(maxRat(ind(cc)))

stat.TCempirical.thetas

% osi = 

%%

TCdiffNull = stat.TCempirical.TCdiffNull;
TCdiff = stat.TCempirical.TCdiff;
maxFRdiffNull = stat.TCempirical.maxFRdiffNull;
maxFRdiff = stat.TCempirical.maxFRdiff;

nullLevel = prctile(TCdiffNull, 95, 2);

figure(2); clf
plot(nullLevel, TCdiff, '.'); hold on
plot(xlim, xlim, 'k')
xlabel("95th percentile for null running modulation")
ylabel("Empirical running modulation")
text(min(xlim) + .8*range(xlim), min(ylim) + .2*range(ylim), 'Favors Null')
text(min(xlim) + .2*range(xlim), min(ylim) + .8*range(ylim), 'Reject Null')

sigTCdiff = TCdiff > nullLevel;
fprintf('%d/%d units have TC modulation (%02.2f)%%\n', sum(sigTCdiff), numel(sigTCdiff), mean(sigTCdiff))
title('TC diff')
% max FR
nullLevel = prctile(maxFRdiffNull, [2.5 97.5], 2);

figure(3); clf
set(gcf, 'Color', 'w')
subplot(1,2,1)
plot(nullLevel(:,1), maxFRdiff, '.'); hold on
plot(xlim, xlim, 'k')
xlabel("2.5th percentile for null max FR modulation")
ylabel("Empirical FR modulation")
text(min(xlim) + .8*range(xlim), min(ylim) + .2*range(ylim), 'Reject Null')
text(min(xlim) + .2*range(xlim), min(ylim) + .8*range(ylim), 'Favors Null')

subplot(1,2,2)
plot(nullLevel(:,2), maxFRdiff, '.'); hold on
plot(xlim, xlim, 'k')

xlabel("97.5th percentile for null max FR modulation")
ylabel("Empirical FR modulation")
text(min(xlim) + .8*range(xlim), min(ylim) + .2*range(ylim), 'Favors Null')
text(min(xlim) + .2*range(xlim), min(ylim) + .8*range(ylim), 'Reject Null')

sigFRmod = (maxFRdiff < nullLevel(:,1) | maxFRdiff > nullLevel(:,2));
fprintf('%d/%d units have Max FR modulation (%02.2f)%%\n', sum(sigFRmod), numel(sigFRmod), mean(sigFRmod))

modUnits = union(find(sigFRmod), find(sigTCdiff));
nMod = numel(modUnits);

fprintf('%d units have potential modulation\n', nMod)

iUnit = 1;
%%


% needs: spls, dfilt, lags, ths
% iUnit = iUnit + 1;
% if iUnit > nMod
%     iUnit = 1;
% end
% modUnits = 1:size(dfilt,2);
% nMod = numel(modUnits);
for iUnit = 1:nMod
    cc = modUnits(iUnit);
    
    fprintf('Unit: %d\n', cc)
    
    % find stable region of firing rate
    unitix = dfilt(:,cc);
    dur = median(D.GratingOffsets(unitix) - D.GratingOnsets(unitix));
    dur = max(dur, .1);
    win = [0.04 dur];
    
    nStim = numel(ths);
    FrateR = nan(numel(lags), nStim);
    FrateS = nan(numel(lags), nStim);
    TCR = nan(nStim, 3);
    TCS = nan(nStim, 3);
    
    figure(1); clf
    subplot(4,2,[1 3]) % no running
    
    spkS = [];
    spkR = [];
    thctr = 0;
    for th = 1:nStim
        iix = find(GratingDirections==ths(th) & ~runningTrial & unitix);
        
        nt = numel(iix);
        spk = squeeze(spks(iix,cc,:));
        if binsize == 1e-3
            [ii,jj] = find(spk);
            plot.raster(lags(jj), ii+thctr, 2, 'Color', 'k' ); hold on
        else
            spk = imboxfilt(spk, [1 3]);
            imagesc(lags, (1:nt)+thctr, spk); hold on
        end
        if size(spk,2) == 1
            spk = spk';
        end
        
        spkS = [spkS; spk];
        thctr = thctr + nt;
        
        R = sum(spk(:,tix),2);
        if isempty(R)
            continue
        end
        TCS(th,1) = mean(R);
        TCS(th,2:3) = bootci(100, @mean, R)';
        FrateS(:,th) = mean(spk);
    end
    title('Stationary')
    ylabel('Trials (sorted by direction)')
    axis tight
    
    subplot(4,2,[5 7]) % no running
    
    thctr = 0;
    for th = 1:nStim
        iix = find(GratingDirections==ths(th) & runningTrial & unitix);
        
        nt = numel(iix);
        spk = squeeze(spks(iix,cc,:));
        
        if binsize == 1e-3
            [ii,jj] = find(spk);
            plot.raster(lags(jj), ii+thctr, 2, 'Color', 'k' ); hold on
        else
            spk = imboxfilt(spk, [1 3]);
            imagesc(lags, (1:nt)+thctr, spk); hold on
        end
        
        if size(spk,2) == 1
            spk = spk';
        end
        
        spkR = [spkR; spk];
        
        thctr = thctr + nt;
        
        R = sum(spk(:,tix),2);
        if isempty(R) || numel(R) < 5
            continue
        end
        TCR(th,1) = mean(R);
        TCR(th,2:3) = bootci(100, @mean, R)';
        FrateR(:,th) = mean(spk);
    end
    title('Running')
    axis tight
    colormap(1-gray)
    ylabel('Trials (sorted by direction)')
    xlabel('Time from Grating Onset')
    
    
    
    vals = [reshape(psthsRunning(:,:,cc), [], 1); reshape(psthsNoRunning(:,:,cc), [], 1)];
    clim = [min(vals(:)) max(vals(:))];
    
    subplot(4,2,2)
    m = FrateS';
    m = imboxfilt(m, [1 3]);
    imagesc(lags, ths, m, clim)
    axis tight
    ylabel('Direction')
    xlabel('Time')
    title('PSTH Stationary')
    
    subplot(4,2,4)
    m = FrateR';
    m = imboxfilt(m, [1 3]);
    imagesc(lags, ths, m, clim)
    axis tight
    ylabel('Direction')
    xlabel('Time')
    title('PSTH Running')
    
    colormap(1-gray)
    
    subplot(4,2,6)
    cmap = lines;
    plot(ths, TCS(:,1), 'k', 'Color', cmap(1,:)); hold on
    fill([ths' fliplr(ths')], [TCS(:,2)' fliplr(TCS(:,3)')], 'k', 'EdgeColor', cmap(1,:))
    
    plot(ths, TCR(:,1), 'k', 'Color', cmap(2,:)); hold on
    fill([ths' fliplr(ths')], [TCR(:,2)' fliplr(TCR(:,3)')], 'k', 'EdgeColor', cmap(2,:))
    title('Tuning Curve')
    xlabel('Direction')
    ylabel('Spike Count')
    xlim([0 360])
    set(gca, 'box', 'off')
    
    subplot(4,2,8)
    plot(lags, nanmean(FrateS, 2)/binsize, 'Color', cmap(1,:)); hold on
    plot(lags, nanmean(FrateR, 2)/binsize, 'Color', cmap(2,:))
    axis tight
    
    title('Mean across directions')
    xlabel('Time from Grat Onset')
    ylabel('Firing Rate')
    
    plot.suplabel(sprintf('Unit %d', cc), 't');
    plot.fixfigure(gcf, 10, [6 8]);
    saveas(gcf, fullfile('Figures', 'HuklabTreadmill', sprintf('examplemod%02.0f.png', cc)))
end


%% Some summaries

figure(1); clf
histogram(arrayfun(@(x) x.llrpval, fitS), 100); hold on
histogram(arrayfun(@(x) x.llrpval, fitR), 100);
legend({'Stationary', 'Running'})
xlabel('LL ratio pval')
title('How many cells are "tuned"?')




%% Plot all tuning curves
fitS = stat.TCfitS;
fitR = stat.TCfitR;
NC = numel(fitS);

sx = ceil(sqrt(NC));
sy = round(sqrt(NC));


figure(1); clf
ax = plot.tight_subplot(sx, sy, 0.001, 0.001);
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
%     title(cc)
end

%%
ntrials = arrayfun(@(x,y) min(x.numTrials, y.numTrials), fitR, fitS);
figure(2); clf
istuned = arrayfun(@(x,y) (x.llrpval < 0.05) & (y.llrpval < 0.05), fitR, fitS);
istuned = istuned & ntrials > 50;
fprintf('%d units that are tuned\n', sum(istuned))

wrappi = @(x) mod(x/pi, 1)*pi;
wrap2pi = @(x) mod(x/2/pi, 1)*2*pi;

mfr = arrayfun(@(x,y) max([x.tuningCurve; y.tuningCurve]), fitS(istuned), fitR(istuned));

bS = arrayfun(@(x) x.paramsML(4), fitS(istuned));
bR = arrayfun(@(x) x.paramsML(4), fitR(istuned));
bSsd = arrayfun(@(x) x.paramsSD(4), fitS(istuned));
bRsd = arrayfun(@(x) x.paramsSD(4), fitR(istuned));

AS = arrayfun(@(x) x.paramsML(3), fitS(istuned));
AR = arrayfun(@(x) x.paramsML(3), fitR(istuned));
ASsd = arrayfun(@(x) x.paramsSD(3), fitS(istuned));
ARsd = arrayfun(@(x) x.paramsSD(3), fitR(istuned));

thS = arrayfun(@(x) x.paramsML(1), fitS(istuned));
thR = arrayfun(@(x) x.paramsML(1), fitR(istuned));
thSsd = arrayfun(@(x) x.paramsSD(1), fitS(istuned));
thRsd = arrayfun(@(x) x.paramsSD(1), fitR(istuned));

thS = wrap2pi(thS);
thR = wrap2pi(thR);

vS = arrayfun(@(x) x.paramsML(2), fitS(istuned));
vR = arrayfun(@(x) x.paramsML(2), fitR(istuned));

lS = arrayfun(@(x) x.paramsML(end), fitS(istuned));
lR = arrayfun(@(x) x.paramsML(end), fitR(istuned));
lSsd = arrayfun(@(x) x.paramsSD(end), fitS(istuned));
lRsd = arrayfun(@(x) x.paramsSD(end), fitR(istuned));

thS(lS > .5) = wrappi(thS(lS > .5));
thR(lR > .5) = wrappi(thR(lR > .5));


subplot(2,2,1)
errorbar(bS, bR, bSsd, bSsd, bRsd, bRsd, 'o', 'Color', .5*[1 1 1], 'MarkerSize', 5, 'MarkerFaceColor', cmap(1,:), 'CapSize', 0); hold on
xlim([0 20])
ylim([0 20])
plot(xlim, xlim, 'k')
title('Baseline')
xlabel('Stationary')
ylabel('Running')

subplot(2,2,2)
% plot
mfr = max(mfr, 10);
errorbar(AS./mfr, AR./mfr, ASsd./mfr, ASsd./mfr, ARsd./mfr, ARsd./mfr, 'o', 'Color', .5*[1 1 1], 'MarkerSize', 5, 'MarkerFaceColor', cmap(1,:), 'CapSize', 0); hold on
% errorbar(AS, AR, ASsd, ASsd, ARsd, ARsd, 'o', 'Color', .5*[1 1 1], 'MarkerSize', 5, 'MarkerFaceColor', cmap(1,:), 'CapSize', 0); hold on
plot(xlim, xlim, 'k')
title('Amplitude (normalized by max FR)')
xlabel('Stationary')
ylabel('Running')
xlim([0 1])
ylim([0 1])


subplot(2,2,3)
errorbar(thS, thR, thSsd, thSsd, thRsd, thRsd, 'o', 'Color', .5*[1 1 1], 'MarkerSize', 5, 'MarkerFaceColor', cmap(1,:), 'CapSize', 0); hold on
plot(xlim, xlim, 'k')
title('Ori Pref')
xlabel('Stationary')
ylabel('Running')
xlim([0 1]*pi)
ylim([0 1]*pi)

subplot(2,2,4)
errorbar(lS, lR, lSsd, lSsd, lRsd, lRsd, 'o', 'Color', .5*[1 1 1], 'MarkerSize', 5, 'MarkerFaceColor', cmap(1,:), 'CapSize', 0); hold on
title('Lambda')
xlabel('Stationary')
ylabel('Running')
xlim([0 1])
ylim([0 1])
plot(xlim, xlim, 'k')


%% Units that became more direction tuned
figure(10); clf

lrat = max(lS, .1)./max(lR, .1);

tunedList = find(istuned);
idx = tunedList(lrat > 2);
n = numel(idx);
sx = ceil(sqrt(n));
sy = round(sqrt(n));
for i = 1:n
    subplot(sx, sy, i)
    cc = idx(i);
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

    set(gca, 'XTick', 0:180:360)
    xlim([0 360])
    text(10, .8*range(ylim)+min(ylim), sprintf('%d', cc), 'fontsize', 8)
end

plot.suplabel('Direction', 'x');
plot.suplabel('Firing Rate', 'y');
plot.suplabel('Units that became more direction tuned', 't')

%% Units that became less direction tuned
figure(10); clf

lrat = max(lS, .1)./max(lR, .1);

tunedList = find(istuned);
idx = tunedList(lrat < .5);
n = numel(idx);
sx = ceil(sqrt(n));
sy = round(sqrt(n));
for i = 1:n
    subplot(sx, sy, i)
    cc = idx(i);
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

    set(gca, 'XTick', 0:180:360)
    xlim([0 360])
    text(10, .8*range(ylim)+min(ylim), sprintf('%d', cc), 'fontsize', 8)
end

plot.suplabel('Direction', 'x');
plot.suplabel('Firing Rate', 'y');
plot.suplabel('Units that became less direction tuned', 't')


%% Units that increased amplitude
figure(10); clf

amprat = max(AS, .1)./max(AR, .1);

tunedList = find(istuned);
idx = tunedList(amprat < .8);
n = numel(idx);
sx = ceil(sqrt(n));
sy = round(sqrt(n));
for i = 1:n
    subplot(sx, sy, i)
    cc = idx(i);
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

    set(gca, 'XTick', 0:180:360)
    xlim([0 360])
    text(10, .8*range(ylim)+min(ylim), sprintf('%d', cc), 'fontsize', 8)
end

plot.suplabel('Direction', 'x');
plot.suplabel('Firing Rate', 'y');
plot.suplabel('Units that increased amplitude', 't')


%% Units that decreased amplitude
figure(10); clf

amprat = max(AS, .1)./max(AR, .1);

tunedList = find(istuned);
idx = tunedList(amprat < 1.2);
n = numel(idx);
sx = ceil(sqrt(n));
sy = round(sqrt(n));
for i = 1:n
    subplot(sx, sy, i)
    cc = idx(i);
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

    set(gca, 'XTick', 0:180:360)
    xlim([0 360])
    text(10, .8*range(ylim)+min(ylim), sprintf('%d', cc), 'fontsize', 8)
end

plot.suplabel('Direction', 'x');
plot.suplabel('Firing Rate', 'y');
plot.suplabel('Units that decreased amplitude', 't')
%% plot tuning curves sorted

figure(33); clf

thetas = linspace(0, 360, 100);

S = cell2mat(arrayfun(@(x,y) x.tuningFun(thetas)./y, fitS(istuned)', mfr', 'uni', 0));
R = cell2mat(arrayfun(@(x,y) x.tuningFun(thetas)./y, fitR(istuned)', mfr', 'uni', 0));
n = sum(istuned);

[~, ind] = sort(vS);

figure(10); clf
subplot(1,3,1)
imagesc(thetas, 1:n, S(ind,:))
title('Stationary')
subplot(1,3,2)
imagesc(thetas, 1:n, R(ind,:))
title('Running')
subplot(1,3,3)
imagesc(thetas, 1:n, R(ind,:)-S(ind,:))
title('Difference')

[~, ind] = sort(thR);

figure(11); clf
subplot(1,3,1)
imagesc(thetas, 1:n, S(ind,:))
title('Stationary')
subplot(1,3,2)
imagesc(thetas, 1:n, R(ind,:))
title('Running')
subplot(1,3,3)
imagesc(thetas, 1:n, R(ind,:)-S(ind,:))
title('Difference')

%% Firing rate analysis

istuned = true(numel(fitS), 1);
figure(10); clf

x = arrayfun(@(x) mean(x.tuningCurve), fitS(istuned));
y = arrayfun(@(x) mean(x.tuningCurve), fitR(istuned));

subplot(1,2,1)
plot(x, y, 'o'); hold on
xlabel('Stationary')
ylabel('Running')
title('Mean Firing Rate')
plot(xlim, xlim, 'k')

mxci = bootci(100, @median, x);
myci = bootci(100, @median, y);

fprintf('MEAN FIRING RATE\n')
fprintf('Stationary FR median = %02.2f [%02.2f, %02.2f]\n', median(x), mxci(1), mxci(2))
fprintf('Running FR median = %02.2f [%02.2f, %02.2f]\n', median(y), myci(1), myci(2))

[pval, h, stats] = ranksum(x, y);

fprintf('wilcoxon pval = %02.5f\n', pval)

figure(2); clf
set(gcf, 'Color', 'w')

subplot(1,2,1)

m = geomean(y./x);
mci = bootci(1000, @geomean, y./x);

fprintf('Ratio of FR = %02.2f [%02.2f, %02.2f]\n', m, mci(1), mci(2))

[cnt, bins] = histcounts(y./x, 100);
bins = (bins(1:end-1) + bins(2:end))/2;

bar(bins, cnt, 'FaceColor', .5*[1 1 1]);
hold on

fill(mci([1 1 2 2]), [ylim, fliplr(ylim)], 'r', 'FaceAlpha', .5)

xlabel('mean FR Ratio Running : Stationary')
ylabel('Count')


% MAX FIRING RATE
fprintf('MAX FIRING RATE\n')

figure(10);

subplot(1,2,2)

x = arrayfun(@(x) max(x.tuningCurve), fitS(istuned));
y = arrayfun(@(x) max(x.tuningCurve), fitR(istuned));

plot(x, y, 'o'); hold on
xlabel('Stationary')
ylabel('Running')
title('Max Firing Rate')
plot(xlim, xlim, 'k')

mxci = bootci(100, @median, x);
myci = bootci(100, @median, y);
fprintf('Stationary FR median = %02.2f [%02.2f, %02.2f]\n', median(x), mxci(1), mxci(2))
fprintf('Running FR median = %02.2f [%02.2f, %02.2f]\n', median(y), myci(1), myci(2))

[pval, h, stats] = ranksum(x, y);
fprintf('wilcoxon pval = %02.5f\n', pval)

figure(2);

subplot(1,2,2)
m = geomean(y./x);
mci = bootci(1000, @geomean, y./x);

fprintf('Ratio of FR = %02.2f [%02.2f, %02.2f]\n', m, mci(1), mci(2))

[cnt, bins] = histcounts(y./x, 100);
bins = (bins(1:end-1) + bins(2:end))/2;

bar(bins, cnt, 'FaceColor', .5*[1 1 1]);
hold on

fill(mci([1 1 2 2]), [ylim, fliplr(ylim)], 'r', 'FaceAlpha', .5)

xlabel('Max FR Ratio Running : Stationary')
ylabel('Count')



%%

figure(33); clf

thetas = fitS(1).thetas;

S = cell2mat(arrayfun(@(x,y) x.tuningCurve(:)'./y, fitS(istuned)', mfr', 'uni', 0));
R = cell2mat(arrayfun(@(x,y) x.tuningCurve(:)'./y, fitR(istuned)', mfr', 'uni', 0));
n = sum(istuned);

[~, ind] = sort(vS);

figure(10); clf
subplot(1,3,1)
imagesc(thetas, 1:n, S(ind,:))
title('Stationary')
subplot(1,3,2)
imagesc(thetas, 1:n, R(ind,:))
title('Running')
subplot(1,3,3)
imagesc(thetas, 1:n, R(ind,:)-S(ind,:))
title('Difference')

[~, ind] = sort(thR);

figure(11); clf
subplot(1,3,1)
imagesc(thetas, 1:n, S(ind,:))
title('Stationary')
subplot(1,3,2)
imagesc(thetas, 1:n, R(ind,:))
title('Running')
subplot(1,3,3)
imagesc(thetas, 1:n, R(ind,:)-S(ind,:))
title('Difference')




%%



thetas = linspace(0, 360, 100);
for cc = find(istuned)
    plot3(cc*ones(100,1), thetas, fitS(cc).tuningFun(thetas)./max(fitS(cc).tuningFun(thetas))); hold on
end

%%
figure(10); clf;
plot(AS, mfr, '.'); hold on
plot(AR, mfr, '.')
xlabel('Amplitude')
ylabel('Max Firing Rate')

figure(11); clf
plot(bS, AS, '.'); hold on
plot(bR, AR, '.')

figure(12); clf
% plot(

%% Step over cells, plot PSTH as image
cc = cc + 1;
if cc > NC
    cc = 1;
end
vals = [reshape(psthsRunning(:,:,cc), [], 1); reshape(psthsNoRunning(:,:,cc), [], 1)];
clim = [min(vals(:)) max(vals(:))];
figure(1); clf
subplot(1,2,1)

imagesc(lags, ths, psthsRunning(:,:,cc)', clim)
ylabel('Direction')
xlabel('Time')
title('Running')

subplot(1,2,2)
imagesc(lags, ths, psthsNoRunning(:,:,cc)', clim)
ylabel('Direction')
xlabel('Time')
title('No Running')
colormap(plot.viridis)


%% plot Tuning Curves
win = [0.04 .4];
iix = lags > win(1) & lags < win(2);
tdur = lags(find(iix,1,'last'))-lags(find(iix,1,'first'));
tcRun = squeeze(nansum(psthsRunning(iix,:,:)))/tdur;
tcNoRun = squeeze(nansum(psthsNoRunning(iix,:,:)))/tdur;

sx = ceil(sqrt(NC));
sy = round(sqrt(NC));


figure(1); clf
ax = plot.tight_subplot(sx, sy, 0.01, 0.01);
for cc = 1:NC
    fprintf("Unit %d/%d\n", cc, NC)
%     subplot(sx, sy, cc)
%     inds = tcRun(:,cc)>0;
    if tcRun(1,cc) == 0
        inds = 2:numel(ths);
    else
        inds = 1:(numel(ths)-1);
    end
    set(gcf, 'currentaxes', ax(cc))
    plot(ths(inds), tcRun(inds,cc), 'k', 'Linewidth', 2); hold on
    plot(ths(inds), tcNoRun(inds,cc), 'r', 'Linewidth', 2); hold on
    set(gca, 'XTick', 0:180:360)
    xlim([0 360])
    axis off
%     title(cc)
end

% plot.fixfigure(gcf, 12, [14 14])
% a = plot.suplabel('Spike Rate', 'y'); 
% a.FontSize = 20;
% a = plot.suplabel('Direction', 'x');
% a.FontSize = 20;

%%
figure(2); clf
set(gcf, 'Color', 'w')
plot(max(tcNoRun), max(tcRun), 'ow', 'MarkerFaceColor', .5*[1 1 1])
hold on
plot(xlim, xlim, 'k')
xlabel('Max Rate (Stationary)')
ylabel('Max Rate (Running)')
%%
cc = 1;

%%

cc = cc + 1;
if cc > NC
    cc = 1;
end

figure(10); clf
nth = numel(unique(D.GratingDirections));
cmap = parula(nth);
ax = plot.tight_subplot(2, nth, 0.01, 0.01);
    
vals = [reshape(psthsRunning(:,:,cc), [], 1); reshape(psthsNoRunning(:,:,cc), [], 1)]/binsize;
clim = [min(vals(:)) max(vals(:))];


for ith = 1:nth
    
    set(gcf, 'currentAxes', ax(ith));
    plot(lags, imgaussfilt(psthsNoRunning(:,ith,cc)/binsize, 2), 'Color', cmap(ith,:), 'Linewidth', 2); hold on
    clr = (cmap(ith,:) + [1 1 1])/2;
    plot(lags, imgaussfilt(psthsRunning(:,ith,cc)/binsize, 2), '-', 'Color', clr, 'Linewidth', 2); hold on
    ylim(clim)
    axis off
    if ith==1
        text(lags(1), .9*clim(2), sprintf('Unit: %d', cc))
        text(lags(1), .8*clim(2), 'Running', 'Color', clr)
        text(lags(1), .7*clim(2), 'No Running', 'Color', cmap(ith,:))
    end
    set(gcf, 'currentAxes', ax(ith+nth));
    [dx, dy] = pol2cart(ths(ith)/180*pi, 1);
    q = quiver(0,0,dx,dy,'Color', cmap(ith,:), 'Linewidth', 5, 'MaxHeadSize', 2); hold on
%     plot([0 dx], [0 dy], 'Color', cmap(ith,:), 'Linewidth', 5); hold on
%     plot(dx, dy, 'o', 'Color', cmap(ith,:), 'MarkerFaceColor', cmap(ith,:), 'MarkerSize', 20)
    
%     R = [cos(pi/2) sin(pi/2); -sin(pi/2) -cos(pi/2)];
    
    
%     for i = [90 270]
%         [dx, dy] = pol2cart((ths(ith) + i)/180*pi, .1);
%         plot(-dx, -dy, 'o', 'Color', cmap(ith,:), 'MarkerFaceColor', cmap(ith,:), 'MarkerSize', 10)
% %     S = [1 0; 0 1];
% %     dxdy = [dx dy] * R*S;
% %     plot(dxdy(1), dxdy(2), 
% %     dxdy = [dx dy] * -R*S;
% %     plot(dxdy(1), dxdy(2), 'o', 'Color', cmap(ith,:), 'MarkerFaceColor', cmap(ith,:), 'MarkerSize', 20)
%     
%     end
    xlim([-1 1]*2)
    ylim([-1 1]*2)
    axis off
end

set(gcf, 'Color', 'w')


%%

[Exp,S] = io.dataFactoryTreadmill(6);
% add unit quality (since some analyses require this field)
Exp.osp.cgs = ones(size(Exp.osp.cids))*2;
io.checkCalibration(Exp);

D = io.get_drifting_grating_output(Exp);

exname = Exp.FileTag;
outdir = fullfile(getpref('FREEVIEWING', 'HUKLAB_DATASHARE'), 'processed');
fname = fullfile(outdir, exname);

save(fname, '-v7', '-struct', 'D')

%% copy to server (for python analyses)
old_dir = pwd;

cd(outdir)
server_string = 'jake@bancanus'; %'jcbyts@sigurros';
output_dir = '/home/jake/Data/Datasets/HuklabTreadmill/processed/';

data_dir = getpref('FREEVIEWING', 'PROCESSED_DATA_DIR');
command = 'scp ';
command = [command exname ' '];
command = [command server_string ':' output_dir];

system(command)

fprintf('%s\n', fname)

cd(old_dir)

