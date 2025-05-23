function D = get_drifting_grating_output(Exp, varargin)
% Output for GLM analyses from the drifting grating protocol

ip = inputParser();
ip.addParameter('plot', true)
ip.addParameter('wheelsmoothing', 15)
ip.addParameter('remove_breaks', false)
ip.addParameter('truncate', true)
ip.parse(varargin{:})
%% Grating tuning

if isfield(Exp, 'GratingOnsets')
    D = Exp;
else
    validTrials = io.getValidTrials(Exp, 'DriftingGrating');
    
    tonsets = [];
    toffsets = [];
    directions = [];
    contrasts = [];
    frameContrast = [];
    phases = [];
    frameTimes = [];
    speeds = [];
    freq = [];
    treadTime = [];
    treadSpeed = [];
    wheelpos = [];
    % figure(1); clf
    
    for iTrial = 1:numel(Exp.D)
        
        % treadmill position and time
        pos = Exp.D{iTrial}.treadmill.locationSpace(:,4); % convert to meters
        ttime = Exp.D{iTrial}.treadmill.locationSpace(:,1);
        
        bad=isnan(pos)|isnan(ttime);
        pos=pos(~bad);
        ttime=ttime(~bad);

        % find the initial messed up period
        init = find(diff(diff(pos)==0)==-1, 1)+1;
        if isempty(init)
            init = 3; % third sample
        end

        %Bypass - dpr 6/5/2023
        %init=min([3, length(pos)-1]);
        trialStart = ttime(init);
        
        
        treadIx = find(ttime > trialStart);
        
        pos = pos - pos(treadIx(1));
        
        dpos = diff(imgaussfilt(pos(treadIx), ip.Results.wheelsmoothing)) ./ diff(ttime(treadIx));
        dpos = [0; dpos];
        spd = dpos;
        
        if ip.Results.plot
            figure(77); clf
            subplot(2,1,1)
            plot(ttime - ttime(1), pos, 'k', 'Linewidth', 2); hold on
            fill(ttime([1 1 init init]) - ttime(1), [ylim fliplr(ylim)], 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
            ylabel('Wheel Position')
            xlabel('Trial Time')
            legend({'Wheel', 'Bad'})
            title(sprintf('Trial %d', iTrial), 'Fontweight', 'normal')
            
            subplot(2,1,2)
            plot(ttime(treadIx)-ttime(1), dpos, 'k', 'Linewidth', 2)
            xlim(ttime([1 end])-ttime(1))
            fprintf('Trial %d) removing %02.2f seconds from start\n', iTrial, trialStart - ttime(1))
            drawnow
            
%             if max(spd) > 60
%                 pause
%             end
        end
        
        
        
        bad = isnan(Exp.D{iTrial}.treadmill.locationSpace(treadIx,2));
        %tmp removing for VSS 19/5/2024
        %bad = imboxfilt(double(bad), 101)>0;
        spd(bad) = nan;
        
        ttime(bad) = nan;
        
        ylim([0 5])
        
        assert(numel(ttime(treadIx)) == numel(spd))
        treadTime = [treadTime; ttime(treadIx)];
        treadSpeed = [treadSpeed; spd];
        wheelpos = [wheelpos; pos(treadIx)];
        
        if ~ismember(iTrial,validTrials)
            continue
        end
        
        % Stimulus Fields:
        % time, orientation, cpd, phase, direction, speed, contrast
        contrast = Exp.D{iTrial}.PR.NoiseHistory(:,7);
        
        onsets = find(diff(contrast)>0)+1;
        offsets = find(diff(contrast)<0)+1;
        
        if contrast(1) > 0
            onsets = [1; onsets];
        end
        
        if contrast(end) > 0
            offsets = [offsets; numel(contrast)];
        end
        
        TrialOnsets = Exp.D{iTrial}.PR.NoiseHistory(onsets,1);
        TrialOffsets = Exp.D{iTrial}.PR.NoiseHistory(offsets,1);
        TrialDirections = Exp.D{iTrial}.PR.NoiseHistory(onsets,5);
        TrialSpeeds = Exp.D{iTrial}.PR.NoiseHistory(onsets,6);
        TrialFreq = Exp.D{iTrial}.PR.NoiseHistory(onsets,3);
        
        TrialPhase = Exp.D{iTrial}.PR.NoiseHistory(:,4);
        TrialFrames = Exp.D{iTrial}.PR.NoiseHistory(:,1);
        goodFrames = TrialFrames > trialStart;
        
        goodStim = TrialOnsets > trialStart; % ignore stimuli in the bad zone at start of trial
        
        tonsets = [tonsets; TrialOnsets(goodStim)];
        toffsets = [toffsets; TrialOffsets(goodStim)];
        directions = [directions;  TrialDirections(goodStim)];
        speeds = [speeds; TrialSpeeds(goodStim)];
        freq = [freq; TrialFreq(goodStim)];
        contrasts = [contrasts; contrast(goodStim)];
        
        frameTimes = [frameTimes; TrialFrames(goodFrames)];
        frameContrast = [frameContrast; contrast(goodFrames)];
        phases = [phases; TrialPhase(goodFrames)];
        
    end
    
    frameTimes = Exp.ptb2Ephys(frameTimes);
    tonsets = Exp.ptb2Ephys(tonsets);
    toffsets = Exp.ptb2Ephys(toffsets);
    treadTime = Exp.ptb2Ephys(treadTime);
    
    % --- spike times / saccade times / eye pos
    st = Exp.spikeTimes;
    clu = Exp.spikeIds;
    eyeTime = Exp.vpx2ephys(Exp.vpx.smo(:,1));
    eyePos = Exp.vpx.smo(:,2:4);
    
    % build struct
    D = struct();
    D.spikeTimes = st;
    D.spikeIds = clu;
    D.treadTime = treadTime;
    D.treadPos = wheelpos;
    D.treadSpeed = treadSpeed;
    D.GratingOnsets = tonsets;
    D.GratingOffsets = toffsets;
    D.GratingDirections = directions;
    D.GratingFrequency = freq;
    D.GratingSpeeds = speeds;
    D.GratingContrast = contrasts;
    D.frameTimes = frameTimes;
    D.framePhase = phases;
    D.frameContrast = frameContrast;
    D.eyeTime = eyeTime;
    D.eyePos = eyePos;
    D.eyeLabels = Exp.vpx.Labels;
    ix = hypot(D.eyePos(:,1), D.eyePos(:,2)) > 20;
    D.eyePos(ix,:) = nan;
    
end

plotGratingData(D)

if ip.Results.remove_breaks
    % --- cut out breaks in data (from other trial types)
    breaks = diff(D.GratingOnsets);
    breakStart = D.GratingOnsets(breaks>10) + 1;
    breakStop = D.GratingOnsets(find(breaks>10)+1) - 1;
    
    
    t = breakStart;
    plot((t*[1 1])', (ones(numel(t),1)*ylim)', 'k', 'Linewidth', 2)
    t = breakStop;
    plot((t*[1 1])', (ones(numel(t),1)*ylim)', 'k--', 'Linewidth', 2)
    drawnow
    
    % remove spike times during brea
    lastOffset = D.GratingOffsets(end) + 2;
    firstOnset = D.GratingOnsets(1) - 2;
    ix = getTimeIdx(D.spikeTimes, breakStart, breakStop) | D.spikeTimes > lastOffset | D.spikeTimes < firstOnset;
    fprintf("Removing %d spikes that occured outside the stimulus trials\n", sum(ix));
    D.spikeTimes(ix) = [];
    D.spikeIds(ix) = [];
    
    ix = getTimeIdx(D.eyeTime, breakStart, breakStop) | D.eyeTime > lastOffset | D.eyeTime < firstOnset;
    fprintf("Removing %d eye pos samples that occured outside the stimulus trials\n", sum(ix));
    D.eyeTime(ix) = [];
    D.eyePos(ix,:) = [];
    D.eyeLabels(ix) = [];
    
    ix = getTimeIdx(D.treadTime, breakStart, breakStop) | D.treadTime > lastOffset | D.treadTime < firstOnset;
    fprintf("Removing %d eye pos samples that occured outside the stimulus trials\n", sum(ix));
    D.treadTime(ix) = [];
    D.treadSpeed(ix) = [];
    
    ix = getTimeIdx(D.frameTimes, breakStart, breakStop) | D.frameTimes > lastOffset | D.frameTimes < firstOnset;
    fprintf("Removing %d frames that occured outside the stimulus trials\n", sum(ix));
    D.frameTimes(ix) = [];
    D.framePhase(ix) = [];
    D.frameContrast(ix) = [];
    
    plotGratingData(D)
    t = breakStart;
    plot((t*[1 1])', (ones(numel(t),1)*ylim)', 'k', 'Linewidth', 2)
    t = breakStop;
    plot((t*[1 1])', (ones(numel(t),1)*ylim)', 'k--', 'Linewidth', 2)
    drawnow
    
    nBreaks = numel(breakStop);
    fprintf('Adjusting times to remove gaps\n')
    for iBreak = 1:nBreaks
        
        breakDur = breakStop(iBreak) - breakStart(iBreak);
        
        plotGratingData(D)
        plot(breakStart(iBreak)*[1 1], ylim, 'b', 'Linewidth', 2)
        plot(breakStop(iBreak)*[1 1], ylim, 'b', 'Linewidth', 2)
        drawnow
        
        % adjust all subsequent times
        D.spikeTimes(D.spikeTimes > breakStop(iBreak)) = D.spikeTimes(D.spikeTimes > breakStop(iBreak)) - breakDur;
        D.frameTimes(D.frameTimes > breakStop(iBreak)) = D.frameTimes(D.frameTimes > breakStop(iBreak)) - breakDur;
        D.GratingOnsets(D.GratingOnsets > breakStop(iBreak)) = D.GratingOnsets(D.GratingOnsets > breakStop(iBreak)) - breakDur;
        D.GratingOffsets(D.GratingOffsets > breakStop(iBreak)) = D.GratingOffsets(D.GratingOffsets > breakStop(iBreak)) - breakDur;
        D.eyeTime(D.eyeTime > breakStop(iBreak)) = D.eyeTime(D.eyeTime > breakStop(iBreak)) - breakDur;
        D.treadTime(D.treadTime > breakStop(iBreak)) = D.treadTime(D.treadTime > breakStop(iBreak)) - breakDur;
        
        if iBreak < nBreaks
            breakStop((iBreak + 1):end) = breakStop((iBreak + 1):end) - breakDur;
            breakStart((iBreak + 1):end) = breakStart((iBreak + 1):end) - breakDur;
        end
    end
    
    plotGratingData(D)
    drawnow
end


if ip.Results.truncate
    % remove all events that happen before the first spike or after the
    % last spike
    
    % first spike
    firstOnset = min(D.spikeTimes);
    lastOffset = max(D.spikeTimes); 
    
    ix = D.eyeTime > lastOffset | D.eyeTime < firstOnset;
    fprintf("Removing %d eye pos samples that occured before the first spike or after the last\n", sum(ix));
    D.eyeTime(ix) = [];
    D.eyePos(ix,:) = [];
    D.eyeLabels(ix) = [];
    
    ix = D.treadTime > lastOffset | D.treadTime < firstOnset;
    fprintf("Removing %d eye pos samples that occured before the first spike or after the last\n", sum(ix));
    D.treadTime(ix) = [];
    D.treadSpeed(ix) = [];
    D.treadPos(ix) = []; %DPR added 10/2024
    
    ix = D.frameTimes > lastOffset | D.frameTimes < firstOnset;
    fprintf("Removing %d frames that occured outside before the first spike or after the last\n", sum(ix));
    D.frameTimes(ix) = [];
    D.framePhase(ix) = [];
    D.frameContrast(ix) = [];
    
    
    ix = D.GratingOnsets < firstOnset | D.GratingOffsets > lastOffset;
    
    gratingFields = {'GratingOnsets', 'GratingOffsets', 'GratingDirections', ...
            'GratingFrequency', 'GratingSpeeds', 'GratingContrast'};
        
    for ifield = 1:numel(gratingFields)
        D.(gratingFields{ifield})(ix) = [];
    end
    
    plotGratingData(D)
    drawnow
end

D.sessNumGratings = ones(size(D.GratingOnsets));
D.sessNumSpikes = ones(size(D.spikeTimes));

pupil = D.eyePos(:,3);
nanix = isnan(pupil);
pad = filtfilt(ones(50,1)/50, 1, double(nanix));
pupil(pad>0) = nan;
D.eyePos(:,3) = pupil;

function plotGratingData(D)
    figure(111); clf
    subplot(3,1,1)
    plot.raster(D.spikeTimes, D.spikeIds, 1); hold on
    t = D.GratingOnsets;
    plot((t*[1 1])', (ones(numel(t),1)*ylim)', 'r', 'Linewidth', 2); hold on
    t = D.GratingOffsets;
    plot((t*[1 1])', (ones(numel(t),1)*ylim)', 'g', 'Linewidth', 2)
    title('Spikes + Gratings')
    
    subplot(3,1,2)
    plot(D.eyeTime, D.eyePos(:,1), 'k')
    hold on
    plot(D.eyeTime, D.eyePos(:,2), 'Color', .5*[1 1 1])
    plot(D.eyeTime(D.eyeLabels==2), D.eyePos(D.eyeLabels==2, 1), 'r.')
    title('Eye Position')
    
    subplot(3,1,3)
    plot(D.treadTime, D.treadSpeed, 'k', 'Linewidth', 2); hold on
    title('Treadmill')




