function Stat = do_spike_count_eye_analyses(D, opts)
% Stat = do_spike_count_analyses(D, opts)
% 
% defopts.winstart = 0.035;
% defopts.prewin = .2;
% defopts.postwin = .1;
% defopts.binsize = .01;
% defopts.run_thresh = 3;
% defopts.debug = false;
% defopts.nboot = 100;
% defopts.spike_rate_thresh = 1;

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


Stat = struct();
% reproduce some the metrics from the Allen Institute analyses
Stat.pvis = nan(NC, 3); % p-value from anova1 (spikes vs condition), where condition is either all or direction only
Stat.cid = nan(NC, 1);  % unit id
Stat.nori = nan(NC,3);  % number of trials per orientation
Stat.ndir = nan(NC, 3); % number of trials per direction
Stat.runmod = nan(NC,2); % running modulation (Allen institute metric)
Stat.murun = nan(NC,2); % mean running (at pref stimulus)
Stat.mustat = nan(NC,2); % mean not running (at pref stimulus)
Stat.prun = nan(NC,2); % pvalue for running (Allen Institute metric)

% bootstrapped firing rates under different conditions
Stat.frBaseR = nan(NC, 3); % baseline firing rate during running (errorbars and mean)
Stat.frBaseS = nan(NC, 3);
Stat.frStimR = nan(NC, 3); 
Stat.frStimS = nan(NC, 3);
Stat.frPrefR = nan(NC, 3);
Stat.frPrefS = nan(NC, 3);
Stat.speedR  = nan(NC, 3);
Stat.speedS  = nan(NC, 3);

% Tuning curves, PSTHS 
total_directions = unique(D.GratingDirections(~isnan(D.GratingDirections)));
total_orientations = unique(mod(total_directions, 180));
num_directions = numel(total_directions);
num_orientations = numel(total_orientations);

Stat.meanrate = nan(NC,1); % mean firing rate
Stat.baselinerate = nan(NC,1);

Stat.directions = total_directions;

Stat.dratemarg = nan(NC,num_directions, 3); % direction tuning curve marginalized over TF and SF
Stat.dratemargR = nan(NC,num_directions, 3); % direction tuning curve marginalized over TF and SF
Stat.dratemargS = nan(NC,num_directions, 3); % direction tuning curve marginalized over TF and SF

Stat.dratebest = nan(NC,num_directions, 3); % direction tuning curve at best SF
Stat.dratebestR = nan(NC,num_directions, 3); % direction tuning curve at best SF
Stat.dratebestS = nan(NC,num_directions, 3); % direction tuning curve at best SF

Stat.drateweight = nan(NC,num_directions, 3); % direction tuning curve weighted by TF / SF tuning
Stat.drateweightR = nan(NC,num_directions, 3); % direction tuning curve weighted by TF / SF tuning
Stat.drateweightS = nan(NC,num_directions, 3); % direction tuning curve weighted by TF / SF tuning

Stat.orientations = total_orientations; % orientations
Stat.oratemarg = nan(NC,num_orientations, 3); % orientation
Stat.oratemargR = nan(NC,num_orientations, 3); % orientation
Stat.oratemargS = nan(NC,num_orientations, 3); % orientation

Stat.oratebest = nan(NC,num_orientations, 3);
Stat.oratebestR = nan(NC,num_orientations, 3);
Stat.oratebestS = nan(NC,num_orientations, 3);

Stat.orateweight = nan(NC,num_orientations, 3);
Stat.orateweightR = nan(NC,num_orientations, 3);
Stat.orateweightS = nan(NC,num_orientations, 3);

Stat.corrSacFRRho = nan(NC,1);
Stat.corrVarXYPval= nan(NC,1); 

Stat.corrSacFRRho= nan(NC,1); Stat.corrSacFRPval= nan(NC,1); 
Stat.corrSacMagFRRho= nan(NC,1); Stat.corrSacmagFRPval= nan(NC,1); 
Stat.corrSizeFRRho= nan(NC,1); Stat.corrSizeFRPval= nan(NC,1); 
Stat.corrVarXFRRho= nan(NC,1); Stat.corrVarXFRPval= nan(NC,1); 
Stat.corrVarYFRRho= nan(NC,1); Stat.corrVarYFRPval= nan(NC,1); 
Stat.corrVarXYFRRho= nan(NC,1); Stat.corrVarXYFRPval= nan(NC,1);

Stat.corrSacFRRhoSt= nan(NC,1); Stat.corrSacFRPvalSt= nan(NC,1); 
Stat.corrSacMagFRRhoSt= nan(NC,1); Stat.corrSacmagFRPvalSt= nan(NC,1); 
Stat.corrSizeFRRhoSt= nan(NC,1); Stat.corrSizeFRPvalSt= nan(NC,1); 
Stat.corrVarXFRRhoSt= nan(NC,1); Stat.corrVarXFRPvalSt= nan(NC,1); 
Stat.corrVarYFRRhoSt= nan(NC,1); Stat.corrVarYFRPvalSt= nan(NC,1); 
Stat.corrVarXYFRRhoSt= nan(NC,1); Stat.corrVarXYFRPvalSt= nan(NC,1);

Stat.corrSacFRRhoB= nan(NC,1); Stat.corrSacFRPvalB= nan(NC,1); 
Stat.corrSacMagFRRhoB= nan(NC,1); Stat.corrSacmagFRPvalB= nan(NC,1); 
Stat.corrSizeFRRhoB= nan(NC,1); Stat.corrSizeFRPvalB= nan(NC,1); 
Stat.corrVarXFRRhoB= nan(NC,1); Stat.corrVarXFRPvalB= nan(NC,1); 
Stat.corrVarYFRRhoB= nan(NC,1); Stat.corrVarYFRPvalB= nan(NC,1); 
Stat.corrVarXYFRRhoB= nan(NC,1); Stat.corrVarXYFRPvalB= nan(NC,1); 

Stat.corrSacRho= nan(NC,1); Stat.corrSacPval= nan(NC,1); 
Stat.corrSizeRho= nan(NC,1); Stat.corrSizePval= nan(NC,1); 
Stat.corrVarXRho= nan(NC,1); Stat.corrVarXPval= nan(NC,1); 
Stat.corrVarYRho= nan(NC,1); Stat.corrVarYPval= nan(NC,1); 
Stat.corrVarXYRho= nan(NC,1);Stat.corrVarXYPval= nan(NC,1);

Stat.corrSacMagRhoB= nan(NC,1); Stat.corrSacmagPvalB= nan(NC,1); 
Stat.corrSacMagRhoSt= nan(NC,1); Stat.corrSacMagPvalSt= nan(NC,1); 
Stat.corrSacMagRho= nan(NC,1); Stat.corrSacmagPval= nan(NC,1); 

Stat.regRnSac_B= nan(NC,2);  Stat.regRnSac_STATS= nan(NC,4); 
Stat.regRnSpd_B= nan(NC,2);  Stat.regRnSpd_STATS= nan(NC,4); 
Stat.regRnSacSpd_B= nan(NC,3); Stat.regRnSacSpd_STATS= nan(NC,4); 


Stat.SacHzR= nan(NC,3); Stat.SacHzS= nan(NC,3); Stat.SacHzAll= nan(NC,3); 
Stat.SacMgR= nan(NC,3); Stat.SacMgS= nan(NC,3); Stat.SacMgAll= nan(NC,3); 
Stat.PupSzR= nan(NC,3); Stat.PupSzS= nan(NC,3); Stat.PupSzAll= nan(NC,3); 
Stat.EyeVarR= nan(NC,3); Stat.EyeVarS= nan(NC,3); Stat.EyeVarAll= nan(NC,3);

Stat.robs = cell(NC, 1);
Stat.runningspeed = cell(NC, 1);

% correlation with running
Stat.runrho = nan(NC,1);
Stat.runrhop = nan(NC,1);

for cc = 1:NC
%
    fprintf('%d/%d\n', cc, NC)
    cid = cids(cc);
    Stat.cid(cc) = cid;

    unitix = D.spikeIds == cid;
    sessix = unique(D.sessNumSpikes(unitix));

    gtix = find(ismember(D.sessNumGratings, sessix));
    gtix(isnan(D.GratingDirections(gtix))) = [];

    onsets = D.GratingOnsets(gtix);
    winsize = median(D.GratingOffsets(gtix) - D.GratingOnsets(gtix));

    t0 = min(onsets) - 2*winsize;
    st = D.spikeTimes(unitix) - t0;
    onsets = onsets - t0;

    st(st < min(onsets)) = [];
    sp = struct('st', st, 'clu', ones(numel(st),1));

    % quick bin spikes during gratings
    R = binNeuronSpikeTimesFast(sp, onsets, winsize);
    R = R ./ winsize;

%     SKIP IF SPIKE RATE DURING GRATINGS BELOW THRESHOLD
    if mean(R) < opts.spike_rate_thresh % skip if spike rate is less that 
        continue
    end

    % --- Get binned spikes, stimulus, behavior
    [stimconds, robs, behavior, unitopts] = bin_ssunit(D, cid, 'win', [-opts.prewin opts.postwin], 'plot', false, 'binsize', opts.binsize);
    
%     SKIP NO GOOD EYETRACE
    if nnz(~isnan(behavior{3}))==0 % 
        continue
    end

    % bin spikes at 50ms and run an anova
    binmat = (unitopts.lags(:) < 0.05:.05:1.05) - (unitopts.lags(:) < 0:.05:1);
    binmat(:,sum(binmat)==0) = [];
    binmat = binmat ./ sum(binmat);
    
    brobs = robs*binmat;
    
    group = reshape(repmat(stimconds{1}, 1, size(brobs,2)), [], 1);
    Stat.pvis(cc,3) = anova1(brobs(:), group, 'off');

    %time index of stimuli on
    tixSt = unitopts.lags > 0 & unitopts.lags < (unitopts.lags(end) - opts.postwin);
    R = sum(robs(:,tixSt),2) ./ (sum(tixSt)*unitopts.binsize);
    Rall = sum(robs,2)./ (size(robs,2)*unitopts.binsize);

    %time index of baseline (prestim)
    tixB = unitopts.lags < 0 ; 
    frbase = sum(robs(:,tixB),2) ./ (sum(tixB)*unitopts.binsize);
    baseline = mean(frbase);
    Stat.baselinerate(cc) = baseline;


    runspeed = nanmean(behavior{1},2); %#ok<*NANMEAN> 

        % pull out eyedata, moved from main_import_script_3
        goodix = 1:size(robs,2);%getStableRange(sum(robs,2), 'plot', false);
        timebintotal=unitopts.NLags*unitopts.binsize;
    
        eyeFlag = behavior{4}(:,goodix);
        nSac=sum((diff(eyeFlag'))==1)';
        eyeSize = mean(behavior{3}(:,goodix),2);
        eyePosX = behavior{5}(:,goodix);
        eyePosY = behavior{6}(:,goodix);
        varX= var(eyePosX,[],2);
        varY= var(eyePosY,[],2);
        varXY= hypot(varX,varY);


        % Separate trial into pre and during stimulus
        eyeFlagSt = behavior{4}(:,tixSt); %Flag for eye state (saccade/fix/lost etc)
        nSacSt=sum((diff(eyeFlagSt'))==1)';
        eyeFlagB = behavior{4}(:,tixB); %Flag for eye state (saccade/fix/lost etc)
        nSacB=sum((diff(eyeFlagB'))==1)';

        eyeSizeSt = mean(behavior{3}(:,tixSt),2);
        eyeSizeB = mean(behavior{3}(:,tixB),2);
        
        eyePosXSt = behavior{5}(:,tixSt);
        eyePosXB = behavior{5}(:,tixB);
        eyePosYSt = behavior{6}(:,tixSt);
        eyePosYB = behavior{6}(:,tixB);
        
        varXSt= var(eyePosXSt,[],2);
        varYSt= var(eyePosYSt,[],2);
        varXYSt= hypot(varXSt,varYSt);

        varXB= var(eyePosXB,[],2);
        varYB= var(eyePosYB,[],2);
        varXYB= hypot(varXB,varYB);

        runSpd = runspeed;%behavior{1}(goodIx,:);
        spd = mean(runSpd,2);
    
        %% Entire trial
        % For each saccade find magnitude
        sacflags=diff(eyeFlag');
        nsc=nnz(sacflags==1);
        %transpose to loop over time (rows now cols) in order
        saconset=find(sacflags==1);
        sacoffset=find(sacflags==-1);
        [onrw,oncol]=ind2sub(size(sacflags),saconset);
        [offrw,offcol]=ind2sub(size(sacflags),sacoffset);
        %
        clear sacmag sacmagspd sacmagtri
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
        
        [Stat.corrSacMagRho(cc), Stat.corrSacMagPval(cc)] = corr(sacmagspd(~isnan(sacmagspd)), sacmag(~isnan(sacmagspd)), 'type', 'Spearman');
        nt=size(spd,1);
        sacmagtri=sacmagtri(~isnan(sacmagspd));
        for ii=1:nt
            if isempty(sacmag(sacmagtri==ii))
                MSac(ii,1)=0; % no saccade
            else
                MSac(ii,1)=mean(sacmag(sacmagtri==ii));
            end
        end

        
        %% During stim
        % For each saccade find magnitude
        sacflagsSt=diff(eyeFlagSt');
        nscSt=nnz(sacflagsSt==1);
        %transpose to loop over time (rows now cols) in order
        saconset=find(sacflagsSt==1);
        sacoffset=find(sacflagsSt==-1);
        [onrw,oncol]=ind2sub(size(sacflagsSt),saconset);
        [offrw,offcol]=ind2sub(size(sacflagsSt),sacoffset);
        %
        clear sacmagSt sacmagspdSt sacmagtriSt
        jj=1;
        for ii=1:nscSt
            trialmatch=offrw(offcol==oncol(ii)); %same column (trials), find rows(times) when saccade stopped
            lastsind=find(trialmatch>onrw(ii),1);
            if ~isempty(lastsind)
            %     presac onrw(ii) oncl(ii)
            %     offsac onrw(ii) trialmatch(lastsind)
                %now need to switch rows and columns back for eyepos
                xpre=eyePosXSt(oncol(ii),max([onrw(ii)-1 1]));%position before saccade
                ypre=eyePosYSt(oncol(ii),max([onrw(ii)-1 1]));
                xpost=eyePosXSt(oncol(ii),trialmatch(lastsind)+1);
                ypost=eyePosXSt(oncol(ii),trialmatch(lastsind)+1);
        
                sacmagSt(jj,1)=hypot(xpost-xpre,ypost-ypre);
                sacmagspdSt(jj,1)=spd(oncol(ii));
                sacmagtriSt(jj,1)=(oncol(ii));
                jj=jj+1;
            end
        end
        
        [Stat.corrSacMagRhoSt(cc), Stat.corrSacMagPvalSt(cc)] = corr(sacmagspdSt(~isnan(sacmagspdSt)), sacmagSt(~isnan(sacmagspdSt)), 'type', 'Spearman');
        nt=size(spd,1);
         sacmagtriSt=sacmagtriSt(~isnan(sacmagspdSt));
        for ii=1:nt
            if isempty(sacmagSt(sacmagtriSt==ii))
                MSacSt(ii,1)=0; % no saccade
            else
                MSacSt(ii,1)=mean(sacmagSt(sacmagtriSt==ii));
            end
        end
        frstim =R;

        %% Baseline
        % For each saccade find magnitude
        sacflagsB=diff(eyeFlagB');
        nscB=nnz(sacflagsB==1);
        %transpose to loop over time (rows now cols) in order
        saconset=find(sacflagsB==1);
        sacoffset=find(sacflagsB==-1);
        [onrw,oncol]=ind2sub(size(sacflagsB),saconset);
        [offrw,offcol]=ind2sub(size(sacflagsB),sacoffset);
        %
        clear sacmagB sacmagspdB sacmagtriB
        jj=1;
        for ii=1:nscB
            trialmatch=offrw(offcol==oncol(ii)); %same column (trials), find rows(times) when saccade stopped
            lastsind=find(trialmatch>onrw(ii),1);
            if ~isempty(lastsind)
            %     presac onrw(ii) oncl(ii)
            %     offsac onrw(ii) trialmatch(lastsind)
                %now need to switch rows and columns back for eyepos
                xpre=eyePosXB(oncol(ii),max([onrw(ii)-1 1]));%position before saccade
                ypre=eyePosYB(oncol(ii),max([onrw(ii)-1 1]));
                xpost=eyePosXB(oncol(ii),trialmatch(lastsind)+1);
                ypost=eyePosXB(oncol(ii),trialmatch(lastsind)+1);
        
                sacmagB(jj,1)=hypot(xpost-xpre,ypost-ypre);
                sacmagspdB(jj,1)=spd(oncol(ii));
                sacmagtriB(jj,1)=(oncol(ii));
                jj=jj+1;
            end
        end
        
        [Stat.corrSacMagRhoB(cc), Stat.corrSacMagPvalB(cc)] = corr(sacmagspdB(~isnan(sacmagspdB)), sacmagB(~isnan(sacmagspdB)), 'type', 'Spearman');
        nt=size(spd,1);
        sacmagtriB=sacmagtriB(~isnan(sacmagspdB));
        for ii=1:nt
            if isempty(sacmagB(sacmagtriB==ii))
                MSacB(ii,1)=0; % no saccade
            else
                MSacB(ii,1)=mean(sacmagB(sacmagtriB==ii));
            end
        end
    



    
%         % Correlations speed with base firing rate
%         [corrRho(cc), corrPval(cc)] = corr(spd(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
%     
        %behavior correlates of overall rate
        [Stat.corrSacFRRho(cc), Stat.corrSacFRPval(cc)] = corr(nSac(~isnan(spd)), Rall(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSacMagFRRho(cc), Stat.corrSacmagFRPval(cc)] = corr(MSac(~isnan(spd)), Rall(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSizeFRRho(cc), Stat.corrSizeFRPval(cc)] = corr(eyeSize(~isnan(spd)), Rall(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXFRRho(cc), Stat.corrVarXFRPval(cc)] = corr(varX(~isnan(spd)), Rall(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarYFRRho(cc), Stat.corrVarYFRPval(cc)] = corr(varY(~isnan(spd)), Rall(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXYFRRho(cc), Stat.corrVarXYFRPval(cc)] = corr(varXY(~isnan(spd)), Rall(~isnan(spd)), 'type', 'Spearman');

        %behavior correlates of stim-driven firing
        [Stat.corrSacFRRhoSt(cc), Stat.corrSacFRPvalSt(cc)] = corr(nSacSt(~isnan(spd)), frstim(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSacMagFRRhoSt(cc), Stat.corrSacmagFRPvalSt(cc)] = corr(MSacSt(~isnan(spd)), frstim(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSizeFRRhoSt(cc), Stat.corrSizeFRPvalSt(cc)] = corr(eyeSizeSt(~isnan(spd)), frstim(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXFRRhoSt(cc), Stat.corrVarXFRPvalSt(cc)] = corr(varXSt(~isnan(spd)), frstim(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarYFRRhoSt(cc), Stat.corrVarYFRPvalSt(cc)] = corr(varYSt(~isnan(spd)), frstim(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXYFRRhoSt(cc), Stat.corrVarXYFRPvalSt(cc)] = corr(varXYSt(~isnan(spd)), frstim(~isnan(spd)), 'type', 'Spearman');

        %behavior correlates of baseline firing
        [Stat.corrSacFRRhoB(cc), Stat.corrSacFRPvalB(cc)] = corr(nSacB(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSacMagFRRhoB(cc), Stat.corrSacmagFRPvalB(cc)] = corr(MSacB(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSizeFRRhoB(cc), Stat.corrSizeFRPvalB(cc)] = corr(eyeSizeB(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXFRRhoB(cc), Stat.corrVarXFRPvalB(cc)] = corr(varXB(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarYFRRhoB(cc), Stat.corrVarYFRPvalB(cc)] = corr(varYB(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXYFRRhoB(cc), Stat.corrVarXYFRPvalB(cc)] = corr(varXYB(~isnan(spd)), frbase(~isnan(spd)), 'type', 'Spearman');
    

        %behavior correlates of running with firing per cell
        [Stat.corrSacRho(cc), Stat.corrSacPval(cc)] = corr(spd(~isnan(spd)), nSac(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrSizeRho(cc), Stat.corrSizePval(cc)] = corr(spd(~isnan(spd)), eyeSize(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXRho(cc), Stat.corrVarXPval(cc)] = corr(spd(~isnan(spd)), varX(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarYRho(cc), Stat.corrVarYPval(cc)] = corr(spd(~isnan(spd)), varY(~isnan(spd)), 'type', 'Spearman');
        [Stat.corrVarXYRho(cc), Stat.corrVarXYPval(cc)] = corr(spd(~isnan(spd)), varXY(~isnan(spd)), 'type', 'Spearman');
    
        %Regression of saccades/running/both on rate: 
        [Stat.regRnSac_B(cc,:),Stat.regRnSac_BINT{cc},Stat.regRnSac_R{cc},Stat.regRnSac_RINT{cc},Stat.regRnSac_STATS(cc,:)] = regress(R,[nSac ones(nt,1)]);
        [Stat.regRnSpd_B(cc,:),Stat.regRnSpd_BINT{cc},Stat.regRnSpd_R{cc},Stat.regRnSpd_RINT{cc},Stat.regRnSpd_STATS(cc,:)] = regress(R,[spd ones(nt,1)]);
        [Stat.regRnSacSpd_B(cc,:),Stat.regRnSacSpd_BINT{cc},Stat.regRnSacSpd_R{cc},Stat.regRnSacSpd_RINT{cc},Stat.regRnSacSpd_STATS(cc,:)] = regress(R,[nSac spd ones(nt,1)]);

        runTrials = find(spd > opts.run_thresh);
        statTrials = find((spd) < opts.run_thresh);
        mixTrials = [runTrials; statTrials];
        
        nrun = numel(runTrials);
        nstat = numel(statTrials);
        
        % 95% confidence interval
        n = min(nrun, nstat);
        if n>0
    %         frBaseR(cc,:) = prctile(mean(frbase(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
    %         frBaseS(cc,:) = prctile(mean(frbase(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);
    %         
    %         iix = opts.lags > 0.04 & opts.lags < opts.lags(end)-.15;
    %         frstim = sum(robs(:,iix),2) / (max(opts.lags(iix)) - min(opts.lags(iix)));
    %         
    %         frStimR(cc,:) = prctile(mean(frstim(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
    %         frStimS(cc,:) = prctile(mean(frstim(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);
   
        
            Stat.SacHzR(cc,:) = prctile(mean(nSac(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.SacHzS(cc,:) = prctile(mean(nSac(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.SacHzAll(cc,:) = prctile(mean(nSac((randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
        
        
            Stat.SacMgR(cc,:) = prctile(nanmean(MSac(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.SacMgS(cc,:) = prctile(nanmean(MSac(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.SacMgAll(cc,:) = prctile(nanmean(MSac((randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
        
            Stat.PupSzR(cc,:) = prctile(mean(eyeSize(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.PupSzS(cc,:) = prctile(mean(eyeSize(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.PupSzAll(cc,:) = prctile(mean(eyeSize((randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
        
            Stat.EyeVarR(cc,:) = prctile(mean(varXY(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.EyeVarS(cc,:) = prctile(mean(varXY(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;
            Stat.EyeVarAll(cc,:) = prctile(mean(varXY((randi(nstat, [n opts.nboot])))), [2.5 50 97.5])./timebintotal;

        end
        


    % get conditions
    direction = stimconds{1};
    orientation = mod(direction, 180);
    speed = stimconds{2};
    freq = stimconds{3};

    speeds = unique(speed(:))';
    freqs = unique(freq(:))';
    directions = unique(direction(:))';
    orientations = unique(orientation(:))';

    % bin all conditions
    Dmat = direction == directions;
    Omat = orientation == orientations;
    Fmat = freq == freqs;
    Smat = speed == speeds;

    nd = numel(directions);
    nf = numel(freqs);
    ns = numel(speeds);
    nt = numel(R);
    no = numel(orientations);
    Xbig = zeros(nt, nd, nf, ns);
    XbigO = zeros(nt, no, nf, ns);
    for iis = 1:ns
        for iif = 1:nf
            Xbig(:,:,iif,iis) = Dmat .* Fmat(:,iif) .* Smat(:,iis);
            XbigO(:,:,iif,iis) = Omat .* Fmat(:,iif) .* Smat(:,iis);
        end
    end




    

    runspeed = nanmean(behavior{1},2); %#ok<*NANMEAN> 
    
    Stat.meanrate(cc) = mean(R);
    Stat.robs{cc} = R;
    Stat.runningspeed{cc} = runspeed;

    ix = ~isnan(runspeed);
    [rho, pval] = corr(R(ix), runspeed(ix), 'Type', 'Spearman');
    Stat.runrho(cc) = rho;
    Stat.runrhop(cc) = pval;

    statix = runspeed < opts.run_thresh;
    runix = runspeed > opts.run_thresh;

    Xbig = reshape(Xbig, nt, []);

    % get stimulus selectivity
    [i,group] = find(Xbig>0);
    [~, ind] = sort(i);
    group = group(ind);

    pval = anova1(R, group, 'off');
    Stat.pvis(cc,1) = pval;

    % get direction selectivity
    [i,group] = find(Dmat>0);
    [~, ind] = sort(i);
    group = group(ind);

    pval = anova1(R, group, 'off');
    Stat.pvis(cc,2) = pval;
    
    % TUNING CURVES

    % --- direction marginalize over TF/SF (Running)
    direc_ix = ismember(total_directions, directions);

    X = Dmat.*runix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);
    Stat.dratemargR(cc,direc_ix,:) = ci';

    if mean(statix) < mean(runix)
        [~, mxid] = max(ci(2,:));
    end

    
    % --- direction marginalize over TF/SF (Stationary)
    X = Dmat.*statix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);
    Stat.dratemargS(cc,direc_ix,:) = ci';

    if mean(statix) > mean(runix)
        [~, mxid] = max(ci(2,:));
    end
    
%     figure(1); clf
%     plot(directions, squeeze(Stat.dratemargR(cc,direc_ix,:)), 'r'); hold on
%     plot(directions, squeeze(Stat.dratemargS(cc,direc_ix,:)), 'b');

    % --- running modulation at preferred direction (Allen Institute metric)
    Stat.ndir(cc,1) = sum(Dmat(:,mxid));

    Drun = Dmat(:,mxid).*runix;
    murun = sum(Drun.*R) ./ sum(Drun);
    Dstat = Dmat(:,mxid).*statix;
    mustat = sum(Dstat.*R) ./ sum(Dstat);
    C = sign(murun - mustat);
    rmax = max(murun, mustat);
    rmin = min(murun, mustat);
    Stat.runmod(cc,1) = C * (rmax - rmin) / abs(rmin);
    Stat.murun(cc,1) = murun;
    Stat.mustat(cc,1) = mustat;
    [~, Stat.prun(cc,1)] = ttest2(R(Drun>0), R(Dstat>0));


    runPrefTrials=Drun;
    statPrefTrials=Dstat;



    % --- direction tuning at best SF / TF
    X = Xbig;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);
    mu = ci(2,:);
    Ici = reshape(ci', nd, [], 3);
    I = reshape(mu, nd, []);
    [ii,mxid] = find(I==max(I(:)));
    if numel(mxid) > 1
        [~, jj] = max(sum(I));
        mxid = jj;
        [~, ii] = max(I(:,mxid));
    end

    Stat.dratebest(cc,direc_ix,:) = squeeze(Ici(:,mxid,:));


    X = Xbig.*runix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Ici = reshape(ci', nd, [], 3);
    
    Stat.dratebestR(cc,direc_ix,:) = squeeze(Ici(:,mxid,:));

    X = Xbig.*statix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Ici = reshape(ci', nd, [], 3);

    Stat.dratebestS(cc,direc_ix,:) = squeeze(Ici(:,mxid,:));

%     figure(1); clf
%     plot(directions, squeeze(Stat.dratebestR(cc,direc_ix,:)), 'r'); hold on
%     plot(directions, squeeze(Stat.dratebestS(cc,direc_ix,:)), 'b');

    % get running modulation at best single stimulus
    [~, id] = max(ci(2,:));
    Drun = Xbig(:,id).*runix;
    murun = sum(Drun.*R) ./ sum(Drun);
    Dstat = Xbig(:,id).*statix;
    mustat = sum(Dstat.*R) ./ sum(Dstat);
    C = sign(murun - mustat);
    rmax = max(murun, mustat);
    rmin = min(murun, mustat);
    Stat.runmod(cc,2) = C * (rmax - rmin) / abs(rmin);
    Stat.murun(cc,2) = murun;
    Stat.mustat(cc,2) = mustat;
    [~, Stat.prun(cc,2)] = ttest2(R(Drun>0), R(Dstat>0));


    n = reshape(sum(Xbig), nd, []);
    
    Stat.ndir(cc,2) = sum(n(:,mxid));
    
%     [ii,jj] = find(I==max(I(:)));
    Stat.ndir(cc,3) = n(ii,mxid); % number of trials at best stimulus

    % --- direction tuning weighted by SF / TF tuning
    w = max(I)'; w = w ./ sum(w);
    
    X = reshape(reshape(Xbig, nt*nd, [])*w, nt, nd);
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.drateweight(cc,direc_ix,:) = ci';

    X = reshape(reshape(Xbig, nt*nd, [])*w, nt, nd);
    X = X.*runix;
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.drateweightR(cc,direc_ix,:) = ci';
    
    X = reshape(reshape(Xbig, nt*nd, [])*w, nt, nd);
    X = X.*statix;
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.drateweightS(cc,direc_ix,:) = ci';

    % ORIENTATION
    XbigO = reshape(XbigO, nt, []);

    % orientation, marginalized
    X = Omat;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    ori_ix = ismember(total_orientations, orientations);
    Stat.oratemarg(cc,ori_ix,:) = ci';

    % orientation, marginalized (Running)
    X = Omat.*runix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    ori_ix = ismember(total_orientations, orientations);
    Stat.oratemargR(cc,ori_ix,:) = ci';

    % orientation, marginalized (Stationary)
    X = Omat.*statix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    ori_ix = ismember(total_orientations, orientations);
    Stat.oratemargR(cc,ori_ix,:) = ci';
    
    % --- orientation tuning, best SF / TF (Running)

    X = XbigO.*runix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Ici = reshape(ci', no, [], 3);

    Stat.oratebestR(cc,ori_ix,:) = squeeze(Ici(:,mxid,:));
    
    % orientation tuning, best SF / TF (Running)

    X = XbigO.*statix;
    x = X.*R;
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Ici = reshape(ci', no, [], 3);

    Stat.oratebestS(cc,ori_ix,:) = squeeze(Ici(:,mxid,:));

    % -- combine running and stationary
    X = XbigO;
    x = X.*R;

    
    % bootstrap errorbars
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);
    mu = ci(2,:);
    Ici = reshape(ci', no, [], 3);
    I = reshape(mu, no, []);
    Stat.oratebest(cc,ori_ix,:) = squeeze(Ici(:,mxid,:));

    n = reshape(sum(XbigO), no, []);
    
    Stat.nori(cc,2) = sum(n(:,mxid));

    [ii,jj] = find(I==max(I(:)));
    if numel(ii)>1
        [~, jj] = max(sum(I));
        [~, ii] = max(I(:,jj));
    end
    Stat.nori(cc,3) = n(ii,jj); % number of trials at best stimulus
    
    % --- direction tuning weighted by SF / TF tuning
    w = max(I)'; w = w ./ sum(w);
    
    X = reshape(reshape(XbigO, nt*no, [])*w, nt, no);
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.orateweight(cc,ori_ix,:) = ci';

    X = reshape(reshape(XbigO, nt*no, [])*w, nt, no);
    X = X.*runix;
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.orateweightR(cc,ori_ix,:) = ci';
    
    X = reshape(reshape(XbigO, nt*no, [])*w, nt, no);
    X = X.*statix;
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.orateweightS(cc,ori_ix,:) = ci';


    
    X = reshape(reshape(XbigO, nt*no, [])*w, nt, no);
    X = X.*runix;
    x = X.*R;
    bootix = randi(size(x,1), [size(x,1) opts.nboot]);
    a = squeeze(sum(reshape(x(bootix,:), [size(bootix) size(X,2)]),1)) ./ ...
        squeeze(sum(reshape(X(bootix,:), [size(bootix) size(X,2)]),1));
    ci = prctile(a, [2.5 50 97.5]);

    Stat.orateweight(cc,ori_ix,:) = ci';

    runTrials = find(runix); %find(runspeed > thresh);
    statTrials = find(statix); %find(abs(runspeed) < 1);




    nrun = numel(runTrials);
    nstat = numel(statTrials);

    n = min(nrun, nstat);
%     n = max(nrun, nstat);

%Saving mean speeds for comparison with regressor
    Stat.speedR(cc,:) = prctile(mean(spd(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
    Stat.speedS(cc,:) = prctile(mean(spd(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);

    Stat.frBaseR(cc,:) = prctile(mean(frbase(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
    Stat.frBaseS(cc,:) = prctile(mean(frbase(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);

    Stat.frStimR(cc,:) = prctile(mean(R(runTrials(randi(nrun, [n opts.nboot])))), [2.5 50 97.5]);
    Stat.frStimS(cc,:) = prctile(mean(R(statTrials(randi(nstat, [n opts.nboot])))), [2.5 50 97.5]);

    runPref = find(runPrefTrials);
    statPref = find(statPrefTrials);

    nrunPref = numel(runPref);
    nstatPref = numel(statPref);

    n = min(nrunPref, nstatPref);

    if n>0
    Stat.frPrefR(cc,:) = prctile(mean(R(runPref(randi(nrunPref, [n opts.nboot])))), [2.5 50 97.5]);
    Stat.frPrefS(cc,:) = prctile(mean(R(statPref(randi(nstatPref, [n opts.nboot])))), [2.5 50 97.5]);
    end
    
end


