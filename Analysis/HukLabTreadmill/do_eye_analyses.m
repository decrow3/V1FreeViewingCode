function EyeStat=do_eye_analyses(D,opts)

if nargin < 2
    opts = struct();
end

defopts = struct();
defopts.winstart = 0.035;
defopts.prewin = .2;
defopts.postwin = .1;
defopts.binsize = .01;
defopts.run_thresh = 3;
defopts.debug = false;
defopts.nboot = 100;
defopts.spike_rate_thresh = 1;
defopts.baseline_subtract = true; % subtract baseline before computing OSI

opts = mergeStruct(defopts, opts);

cids = unique(D.spikeIds);
NC = numel(cids);


EyeStat = struct();

binsize=opts.binsize;
passwin=[-opts.prewin opts.postwin];

% %% Don't need to do this for each cell, but want to use same time points
% [~, ~, behav, opts] = bin_ssunit(D, unitList(1), 'win', [-.2 .1]);

    
    sessNums = unique(D.sessNumSpikes); %All sessions
    
    ix = ismember(D.sessNumGratings, sessNums);
    
    StimOnset  = D.GratingOnsets(ix);
    StimOffset = D.GratingOffsets(ix);
    StimDir    = D.GratingDirections(ix);
    StimSpeed = D.GratingSpeeds(ix);
    StimFreq = D.GratingFrequency(ix);
    
    treadSessIx = ismember(D.sessNumTread, sessNums) & ~isnan(D.treadTime);
    treadTime = D.treadTime(treadSessIx);
    treadSpeed = D.treadSpeed(treadSessIx);
    
    frameIx = D.frameTimes >= (min(StimOnset) - 1) & D.frameTimes < (max(StimOffset) + 1);
    frameTime = D.frameTimes(frameIx);
    framePhase = D.framePhase(frameIx);
    
    % resample time with new binsize
    newTreadTime = treadTime(1):binsize:treadTime(end);
    newTreadSpeed = interp1(treadTime, treadSpeed, newTreadTime);
    
    newFrameTime = newTreadTime;
    newFramePhase = interp1(frameTime, framePhase, newFrameTime);
    
    
    eyeIx = ismember(D.sessNumEye, sessNums) & ~isnan(D.eyePos(:,3));
    eyeSess = D.sessNumEye(eyeIx);
    eyeTime = D.eyeTime(eyeIx);
    eyePupil = D.eyePos(eyeIx,3);
    eyeX = D.eyePos(eyeIx,1);
    eyeY = D.eyePos(eyeIx,2);
    eyeLabels = D.eyeLabels(eyeIx);
    
    pupil = nan(size(newFrameTime));
    eyeX_ = nan(size(newFrameTime));
    eyeY_ = nan(size(newFrameTime));
    eyeLabels_ = nan(size(newFrameTime));
    eyeSess_ = nan(size(newFrameTime));

    %DPR edit - need to put all eye data onto same indices as frame time for
    %iix in loop later
    

    if ~isempty(eyePupil)
        pupil = interp1(eyeTime,eyePupil,newFrameTime);
        eyeX_ = interp1(eyeTime,eyeX,newFrameTime);
        eyeY_ = interp1(eyeTime,eyeY,newFrameTime);
        eyeLabels_ = interp1(eyeTime,eyeLabels,newFrameTime,'nearest');
        eyeSess_ = interp1(eyeTime,eyeSess,newFrameTime);
    end

    %If newFrameTime creates time points between sessions interpolant will 
    % have a non integer session number
    unsupported=eyeSess_~=round(eyeSess_);
    pupil(unsupported)=nan;
    eyeX_(unsupported)=nan;
    eyeY_(unsupported)=nan;
    eyeLabels_(unsupported)=nan;
    eyeSess_(unsupported)=nan;
    newTreadSpeed(unsupported)=nan;
    newFramePhase(unsupported)=nan;

    
    treadTime = newTreadTime;
    treadSpeed = newTreadSpeed;
    framePhase = newFramePhase;
    
    % find index into onsets and offsets
    [~, ~, idOn] = histcounts(StimOnset, treadTime);
    [~, ~, idOff] = histcounts(StimOffset, treadTime);
    
    % check stim duration and only include gratings that were fully shown
    stimDuration = idOff - idOn;
    durBins = mode(stimDuration);
    
    validStim = find( abs(stimDuration-durBins)<2 & ~isnan(StimDir));
    idOff = idOff(validStim);
    idOn = idOn(validStim);
    StimDur = mode(StimOffset(validStim) - StimOnset(validStim));
    StimDir = StimDir(validStim);
    StimSpeed = StimSpeed(validStim);
    StimFreq = StimFreq(validStim);
    

    win = [passwin(1) StimDur+passwin(2)];
    bins = win(1):binsize:win(2);
    
    nbins = numel(bins)-1;
    
    blags = floor(bins/binsize); blags = blags(1:nbins);
    
    NT = numel(validStim);
    
%     opts = struct();
%     opts.NTrials = NT;
%     opts.NLags = nbins;
    
    % get running speed aligned to stim onset
    runSpeed = nan(NT, nbins);
    GratPhase = nan(NT, nbins);
    PupilArea = nan(NT, nbins);
    eyeFlag = nan(NT, nbins);
    eyePosX = nan(NT, nbins);
    eyePosY = nan(NT, nbins);
    eyeSessSet = nan(NT, nbins); 
    nt = numel(treadSpeed);
    
    for i = 1:NT
        iix = blags + idOn(i);
        valid = iix > 0 & iix < nt;
        runSpeed(i,valid) = treadSpeed(iix(valid));
        GratPhase(i,valid) = framePhase(iix(valid));
        PupilArea(i,valid) = pupil(iix(valid));
        eyeFlag(i,valid) = eyeLabels_(iix(valid));
        eyePosX(i,valid) = eyeX_(iix(valid));
        eyePosY(i,valid) = eyeY_(iix(valid));
        eyeSessSet(i,valid) = eyeSess_(iix(valid));
    end

%Discount trials were track was lost for more than half the time
invalid=((sum(eyeFlag'==4))>(nbins/2));
runSpeed(invalid,:) = nan;
GratPhase(invalid,:) = nan;
PupilArea(invalid,:) = nan;
eyeFlag(invalid,:) = nan;
eyePosX(invalid,:) = nan;
eyePosY(invalid,:) = nan;
eyeSessSet(invalid,:) = nan;

%%

% runSpd = behav{1}(goodIx,:);
spd = mean(runSpeed,2);

% eyeFlag = behav{4}(goodIx,:);
nSac=sum((diff(eyeFlag'))==1)';

 eyeSize = mean(PupilArea,2);

% eyePosX = behav{5}(goodIx,:);
% eyePosY = behav{6}(goodIx,:);
varX= var(eyePosX,[],2);
varY= var(eyePosY,[],2);

varXY= hypot(varX,varY);


%% For each saccade find magnitude
sacflags=diff(eyeFlag');
nsc=nnz(sacflags==1);
%transpose to loop over time (rows now cols) in order
saconset=find(sacflags==1);
sacoffset=find(sacflags==-1);
[onrw,oncol]=ind2sub(size(sacflags),saconset);
[offrw,offcol]=ind2sub(size(sacflags),sacoffset);
%
jj=1;
for ii=1:nsc
    trialmatch=offrw(offcol==oncol(ii)); %same column (trials), find rows(times) when saccade stopped
    lastsind=find(trialmatch>onrw(ii),1);
    if ~isempty(lastsind)
    %     presac onrw(ii) oncl(ii)
    %     offsac onrw(ii) trialmatch(lastsind)
        %now need to switch rows and columns back for eyepos
        xpre=eyePosX(oncol(ii),max([onrw(ii)-1 1]));%position before saccade
        ypre=eyePosY(oncol(ii),max([onrw(ii)-1 1]));
        xpost=eyePosX(oncol(ii),trialmatch(lastsind)+1);
        ypost=eyePosX(oncol(ii),trialmatch(lastsind)+1);

        sacmag(jj,1)=hypot(xpost-xpre,ypost-ypre);
        sacmagspd(jj,1)=spd(oncol(ii));
        sacmagtri(jj,1)=(oncol(ii));
        jj=jj+1;
    end
end


nt=size(spd,1);
for ii=1:nt
 MSac(ii,1)=mean(sacmag(sacmagtri==ii));
end

%% Bad trials, no eyetrack
% badT=sum(isnan(eyeFlag)|eyeFlag==4,2)>=(nbins/2) ... %more than half the track was lost
%     |(varXY<=.1); % no change at all, suspicious 

% nSac(badT)=nan;
% MSac(badT)=nan;
% eyeSize(badT)=nan;
% varX(badT)=nan;
% varY(badT)=nan;
% varXY(badT)=nan;
%%

[EyeStat.corrSacMagRho, EyeStat.corrSacMagPval] = corr(sacmagspd(~isnan(sacmagspd)), sacmag(~isnan(sacmagspd)), 'type', 'Spearman');
%Behaviour correlates with speed
[EyeStat.corrSacRho, EyeStat.corrSacPval] = corr(spd(~isnan(spd)), nSac(~isnan(spd)), 'type', 'Spearman');
[EyeStat.corrSizeRho, EyeStat.corrSizePval] = corr(spd(~isnan(spd)), eyeSize(~isnan(spd)), 'type', 'Spearman');
[EyeStat.corrVarXRho, EyeStat.corrVarXPval] = corr(spd(~isnan(spd)), varX(~isnan(spd)), 'type', 'Spearman');
[EyeStat.corrVarYRho, EyeStat.corrVarYPval] = corr(spd(~isnan(spd)), varY(~isnan(spd)), 'type', 'Spearman');
[EyeStat.corrVarXYRho, EyeStat.corrVarXYPval] = corr(spd(~isnan(spd)), varXY(~isnan(spd)), 'type', 'Spearman');

runTrials = find(spd > opts.run_thresh);
statTrials = find(spd < opts.run_thresh);
mixTrials = [runTrials; statTrials];


nrun = numel(runTrials);
nstat = numel(statTrials);
n = min(nrun, nstat);
SacHz=nSac./abs(diff(win));

EyeStat.SacHzR = prctile(mean(SacHz(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.SacHzS = prctile(mean(SacHz(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.SacHzAll = prctile(mean(SacHz((randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);

EyeStat.PupSzR = prctile(mean(eyeSize(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.PupSzS = prctile(mean(eyeSize(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.PupSzAll = prctile(mean(eyeSize((randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);

EyeStat.EyeVarR = prctile(mean(varXY(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.EyeVarS = prctile(mean(varXY(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.EyeVarAll = prctile(mean(varXY((randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);


% MSac is already back in terms of trials 
% runsactri=runTrials(ismember(runTrials,sacmagtri));
% statsactri=statTrials(ismember(statTrials,sacmagtri));
% mixscatri=mixTrials(ismember(mixTrials,sacmagtri));
% 
% nrunsac = numel(runsactri);
% nstatsac = numel(statsactri);
% nmixsac = numel(mixscatri);

EyeStat.SacMagR = prctile(nanmean(MSac(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.SacMagS = prctile(nanmean(MSac(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);
EyeStat.SacMagAll = prctile(nanmean(MSac((randi(n, [n opts.nboot])))), [2.5 50 97.5]);



%% Anova pairs (need to correct for multiple?)
SacHzRtrials=SacHz(runTrials);
SacHzStrials=SacHz(statTrials);
group=[repmat({'Running'},nrun,1); repmat({'Stationary'},nstat,1)];
[EyeStat.SacHzP,EyeStat.SacHzANOVATAB,EyeStat.SacHzSTATS]=anova1([SacHzRtrials; SacHzStrials],group);



eyeSizeRtrials=eyeSize(runTrials);
eyeSizeStrials=eyeSize(statTrials);
[EyeStat.PupSzP,EyeStat.PupSzANOVATAB,EyeStat.PupSzSTATS]=anova1([eyeSizeRtrials; eyeSizeStrials],group);

varXYRtrials=varXY(runTrials);
varXYStrials=varXY(statTrials);
[EyeStat.EyeVarP,EyeStat.EyeVarANOVATAB,EyeStat.EyeVarSTATS]=anova1([varXYRtrials; varXYStrials],group);

% MSacRtrials=MSac(runsactri);
% MSacStrials=MSac(statsactri);
% group=[repmat({'Running'},length(MSacRtrials),1); repmat({'Stationary'},length(MSacStrials),1)];
% [EyeStat.SacMagP,EyeStat.SacMagANOVATAB,EyeStat.SacMagSTATS]=anova1([MSacRtrials; MSacStrials],group);
MSacRtrials=MSac(runTrials);
MSacStrials=MSac(statTrials);
[EyeStat.SacMagP,EyeStat.SacMagANOVATAB,EyeStat.SacMagSTATS]=anova1([MSacRtrials; MSacStrials],group);




%For histograms, eg
% hold off
% histogram(nSac(runTrials(randi(nrun, [n 1]))),0:15)
% hold on
% histogram(nSac(statTrials(randi(nstat, [n 1]))),0:15)
EyeStat.SacHz=SacHz;
EyeStat.nSac=nSac;
EyeStat.MSac=MSac;
EyeStat.eyeSize=eyeSize;
EyeStat.varXY=varXY;
EyeStat.spd=spd;

EyeStat.runTrials=runTrials;
EyeStat.statTrials=statTrials;






