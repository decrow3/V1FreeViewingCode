function [stim, robs, behavior, opts] = bin_ssunit(D, unitId, varargin)
% BIN SUPER SESSION UNIT SPIKE TRAINS ALIGNED TO STIMULUS ONSET
% [stim, robs, behavior, opts] = bin_ssunit(D, unitID, varargin)
% INPUTS:
%   D:      supersession
%   unitID: the id of the unit to analyze
% OUTPUTS:
%   stim:      {StimDir, StimSpeed, StimFreq};
%   robs:      Binned spikes;
%   behavior:  {runSpeed, GratPhase, PupilArea};
%   opts:      the parameters of the analysis
%
% Optional Arguments:
% plot
% win [seconds before (pass in negative), seconds after]
% binsize (in seconds)

ip = inputParser();
ip.addParameter('plot', true);
ip.addParameter('win', [-.1 .1])
ip.addParameter('binsize', 10e-3)
ip.parse(varargin{:})


%% bin spikes
binsize = ip.Results.binsize; % 10 ms bins
if isnan(unitId)
    unitIx = true(size(D.spikeIds));
else
    unitIx = D.spikeIds == unitId;
end

sessNums = unique(D.sessNumSpikes(unitIx));

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
    if sum(isnan(eyeTime))>0
        badeyetime=isnan(eyeTime);
        eyeTime=eyeTime(~badeyetime);
        eyePupil=eyePupil(~badeyetime);
        eyeX=eyeX(~badeyetime);
        eyeY=eyeY(~badeyetime);
        eyeLabels=eyeLabels(~badeyetime);
        eyeSess=eyeSess(~badeyetime);
    end
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

win = [ip.Results.win(1) StimDur+ip.Results.win(2)];
bins = win(1):binsize:win(2);

nbins = numel(bins)-1;

blags = floor(bins/binsize); blags = blags(1:nbins);

NT = numel(validStim);

opts = struct();
opts.NTrials = NT;
opts.NLags = nbins;

% get running speed aligned to stim onset
runSpeed = nan(NT, nbins);
GratPhase = nan(NT, nbins);
PupilArea = nan(NT, nbins);
eyeFlag = nan(NT, nbins);
eyePosX = nan(NT, nbins);
eyePosY = nan(NT, nbins);
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
end

SpikeTimes = D.spikeTimes(unitIx);
SpikeIds = D.spikeIds(unitIx);
UnitList = unique(SpikeIds);

NC = numel(UnitList); % number of neurons
assert(NC == 1 | isnan(unitId), "You requested one unit and we have more than that")

opts.unitId = unitId;
opts.stimDuration = StimDur;

plotDuringImport = ip.Results.plot;

% count spike times aligned to STIM onset
% fprintf('Counting spikes aligned to stim onset \n')    
% count spikes
[scnt, bins] = decoding.binSpTimes(SpikeTimes, StimOnset(validStim), win, binsize);

    
if plotDuringImport
    figure(10); clf
    % visualize tuning
    [~, ind] = sort(StimDir);
        
    smcnt = imgaussfilt(scnt(ind,:), [5 3], 'FilterSize', [11 3]); % smooth along trials
    
    if binsize > 1e-3
        imagesc(bins, StimDir(ind), smcnt); colormap(1-gray)
    else % if 1ms bins, plot raster
        [iTrial,j] = find(scnt(ind,:));
        plot.raster(bins(j), iTrial, 1);
    end
    
    
    drawnow
end

if plotDuringImport
    title(sprintf('Unit %d', unitId))
    xlabel('Time from stim onset', 'Color', 'k')
    ylabel('Trials (sorted by Direction)', 'Color', 'k')
end

% output
stim = {StimDir, StimSpeed, StimFreq};
robs = scnt;
behavior = {runSpeed, GratPhase, PupilArea, eyeFlag, eyePosX, eyePosY};
opts.lags = bins;
opts.binsize = binsize;
