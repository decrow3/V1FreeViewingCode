%% Additional eye analysis for paper

%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory


%% load session
% gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/gratings/';
% subject = 'gru';
% D = load_subject(subject,gratingpath);
% 
% %% To avoid double counting, we need to remove sessions 15 and 17, which are the same as 14 and 16 but at a different part of the probe
% D.eyeLabels=D.eyeLabels((D.sessNumEye~=17 &D.sessNumEye~=15));
% D.eyePos=D.eyePos((D.sessNumEye~=17 &D.sessNumEye~=15),:);
% D.eyeTime=D.eyeTime((D.sessNumEye~=17 &D.sessNumEye~=15),:);
% D.sessNumEye=D.sessNumEye((D.sessNumEye~=17 &D.sessNumEye~=15),:);
% D.spikeIds=D.spikeIds((D.sessNumSpikes~=17 & D.sessNumSpikes~=15),:);
% D.spikeTimes=D.spikeTimes((D.sessNumSpikes~=17 & D.sessNumSpikes~=15),:);
% D.sessNumSpikes=D.sessNumSpikes((D.sessNumSpikes~=17 & D.sessNumSpikes~=15),:);
% D.treadSpeed=D.treadSpeed((D.sessNumTread~=17 & D.sessNumTread~=15),:);
% D.treadTime=D.treadTime((D.sessNumTread~=17 & D.sessNumTread~=15),:);
% D.sessNumTread=D.sessNumTread((D.sessNumTread~=17 & D.sessNumTread~=15),:);
% D.GratingContrast=D.GratingContrast((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
% D.GratingDirections=D.GratingDirections((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
% D.GratingFrequency=D.GratingFrequency((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
% D.GratingOffsets=D.GratingOffsets((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
% D.GratingOnsets=D.GratingOnsets((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
% D.GratingSpeeds=D.GratingSpeeds((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
% D.sessNumGratings=D.sessNumGratings((D.sessNumGratings~=17 & D.sessNumGratings~=15),:);
%%
% D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
% D.eyePos(D.eyePos(:,3)==0,:)=nan;
% 
% Eyestat.(subject) = do_eye_analyses(D);
% %Stat.(subject) = do_spike_count_analyses(D);

%% Regen eyelabels??
% subject = 'brie';
% D = load_subject(subject,gratingpath);
%%
subject='rocky'
D.sessNumTread=(D.treadSpeed<Inf);
D.sessNumEye=(D.eyeTime<Inf);
D.sessNumTread=abs(double(D.treadSpeed<Inf));
D.sessNumEye=double(D.eyeTime<Inf);

%%
D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
D.eyePos(D.eyePos(:,3)==0,:)=nan;

Eyestat.(subject) = do_eye_analyses(D);

%% Single animal plots
% subject = 'brie';
nrun = numel(Eyestat.(subject).runTrials);
nstat = numel(Eyestat.(subject).statTrials);
n=min([nrun nstat]);


%Histograms
figure(4);clf
    %cmap = getcolormap(subject, false);
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
lw=0;800;up=3500;
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
% [P,ANOVATAB*,STATS]=anova1([SacHzRtrials; SacHzStrials],[ones(length(SacHzRtrials),1);2*ones(length(SacHzStrials),1)]);
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
subjects = {'rocky'};%, 'allen'};
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
   % cmap = getcolormap(subject, false);
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
figure(4);
subplot(131)
yy=ylim;
hold on
plot(SacHzR(:,2),yy(2),'bx')
plot(SacHzS(:,2),yy(2),'rx')
hold off

subplot(132)
yy=ylim;
hold on
plot(SacMagR(:,2),yy(2),'bx')
plot(SacMagS(:,2),yy(2),'rx')
hold off

subplot(133)
yy=ylim;
hold on
plot(100*PupSzR(:,2)./mean([eyeSizeR_all(eyeSizeR_all>0);eyeSizeS_all(eyeSizeS_all>0)]),yy(2),'bx')
plot(100*PupSzS(:,2)./mean([eyeSizeR_all(eyeSizeR_all>0);eyeSizeS_all(eyeSizeS_all>0)]),yy(2),'rx')
hold off

%%

group=[repmat({'Running'},nrun,1); repmat({'Stationary'},nstat,1)];
[SacHzP,SacHzANOVATAB,SacHzSTATS]=anova1([SacHzR_all; SacHzS_all],group);
[PupSzP,PupSzANOVATAB,PupSzSTATS]=anova1([eyeSizeR_all; eyeSizeS_all],group);
[EyeVarP,EyeVarANOVATAB,EyeVarSTATS]=anova1([varXYR_all; varXYS_all],group);
[SacMagP,SacMagANOVATAB,SacMagSTATS]=anova1([MSacR_all; MSacS_all],group);

%%

fprintf("Mean saccade frequency during running is %02.3f Hz [%02.3f, %02.3f] (n=%d trials) \n", SacHzR(2), SacHzR(1), SacHzR(3), nrun)
fprintf("Mean saccade frequency during stationary is %02.3f Hz [%02.3f, %02.3f] (n=%d trials) \n", SacHzS(2), SacHzS(1), SacHzS(3), nstat)
fprintf("Two-way anova saccade frequency in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(SacHzP))

fprintf("Mean saccade magnitude during running is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", SacMagR(2), SacMagR(1), SacMagR(3), nrun)
fprintf("Mean saccade magnitude during stationary is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", SacMagS(2), SacMagS(1), SacMagS(3), nstat)
fprintf("Two-way anova saccade magnitude in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(SacMagP))
meaneyeSize=mean([eyeSizeR_all(eyeSizeR_all>0);eyeSizeS_all(eyeSizeS_all>0)]);
fprintf("Mean pupil size during running is %02.3f%% of mean [%02.3f, %02.3f] (n=%d trials) \n", 100*PupSzR(2)./meaneyeSize, 100*PupSzR(1)./meaneyeSize, 100*PupSzR(3)./meaneyeSize, nrun)
fprintf("Mean pupil size during stationary is %02.3f%% of mean [%02.3f, %02.3f] (n=%d trials) \n", 100*PupSzS(2)./meaneyeSize, 100*PupSzS(1)./meaneyeSize, 100*PupSzS(3)./meaneyeSize, nstat)
fprintf("Two-way anova pupil size in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(PupSzP))

fprintf("Mean eye position variance during running is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", EyeVarR(2), EyeVarR(1), EyeVarR(3), nrun)
fprintf("Mean eye position variance during stationary is %02.3f degrees [%02.3f, %02.3f] (n=%d trials) \n", EyeVarS(2), EyeVarS(1), EyeVarS(3), nstat)
fprintf("Two-way anova eye position variance in running (n=%d) vs stationary (n=%d) trials, p=%d \n", length(runTrials),length(statTrials),(EyeVarP))

%%
fprintf("Pearson correlation coefficient of saccade frequency with running is %02.3f [p = %d] (n=%d trials) \n", corrSacRho, corrSacPval,length(SacHz))
fprintf("Pearson correlation coefficient of saccade magnitude with running is %02.3f [p = %d] (n=%d trials) \n", corrMSacRho, corrMSacPval,length(SacHz))
fprintf("Pearson correlation coefficient of pupil size with running is %02.3f [p = %d] (n=%d trials) \n", corrSizeRho, corrSizePval,length(SacHz))
fprintf("Pearson correlation coefficient of eye position variance with running is %02.3f [p = %d] (n=%d trials) \n", corrVarXYRho, corrVarXYPval,length(SacHz))
% 
% %% Reload Gru's data with sessions 15 and 17 back in for cell-by-cell analysis
% subject = 'brie';
% D = load_subject(subject,gratingpath);
% D.eyePos(D.eyePos(:,3)==0,:)=nan;
% D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
%%

Stat.(subject) = do_spike_count_eye_analyses(D);
%%
Stat = struct();

%subjects = {'gru', 'brie'};
nsubjs = numel(subjects);
for isubj = 1:nsubjs
    subject = subjects{isubj};
    D = load_subject(subject,gratingpath);
    D.eyePos(D.eyePos(:,3)==0,:)=nan;
    D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
    Stat.(subject) = do_spike_count_eye_analyses(D);
end
%% Just low sf
 lowsfset=(D.GratingFrequency<2);
Dsave=D;
%%
D=Dsave;
%%
lowtfset=(D.GratingSpeeds<2);
lowsfset=(D.GratingFrequency<2);
%%
D.GratingOnsets=D.GratingOnsets(lowsfset&lowtfset);
D.GratingOffsets=D.GratingOffsets(lowsfset&lowtfset);
D.GratingDirections=D.GratingDirections(lowsfset&lowtfset);
D.GratingFrequency=D.GratingFrequency(lowsfset&lowtfset);
D.GratingSpeeds=D.GratingSpeeds(lowsfset&lowtfset);
D.GratingContrast=D.GratingContrast(lowsfset&lowtfset);
D.sessNumGratings=D.sessNumGratings(lowsfset&lowtfset);
%%
 subject='rocky_lowtf'
%%
figdir = '/home/huklab/Documents/NPX_pilot/V1FreeViewingCode/Figures';
nboot = 100;



% [207, 179, 144; 195, 97, 66]
% cmap = lines;

% colors = {[(cmap(1,:) + [1 1 1]/2)/2; cmap(1,:)]*255, ...
%     [(cmap(5,:) + [1 1 1]/2)/2; cmap(5,:)]*255, ...
%     [(cmap(2,:) + [1 1 1]/2)/2; cmap(2,:)]*255};

subjects = {'rocky_lowtf'};
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
%    cmap = getcolormap(subject, false);

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
        
        
        good = ~(PupSzR(:,2)==0 | PupSzS(:,2)==0 | isnan(PupSzR(:,2)) | isnan(PupSzS(:,2)));
        
        m = geomean(PupSzR(good,2)./PupSzS(good,2));
        ci = bootci(nboot, @geomean, PupSzR(good,2)./PupSzS(good,2));
        
        fprintf("geometric mean pupil size ratio (Running:Stationary) is %02.3f [%02.3f, %02.3f] (n=%d)\n", m, ci(1), ci(2), sum(good)) 
 
        fprintf('%d/%d (%02.2f%%) increased eye pos variance\n', nIncStim, NC, 100*nIncStim/NC)
        fprintf('%d/%d (%02.2f%%) decreased eye pos variance\n', nDecStim, NC, 100*nDecStim/NC)
        fprintf('Eye position variance: p = %02.10f\n', pvalStim)

        good= ~(isnan(EyeVarR(:,2)) | isnan(EyeVarS(:,2)));
        m = geomean(EyeVarR(good,2)./EyeVarS(good,2));
        ci = bootci(nboot, @geomean, EyeVarR(good,2)./EyeVarS(good,2));
        
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
% NC=NCtillnow;
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

