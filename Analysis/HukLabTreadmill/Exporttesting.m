% Export for gain modeling,
% Structure with data binned to frame rate? Or trial by trial
% Robs, Stim (Direction, Freq, Speed), Treadmill Speed, Eye


%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory


%% load session
gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/V1eg/';
subject = 'npx_V1_calc_eg';
cd('/home/huklab/Documents/NPX_pilot/V1Locomotion/Code/')
D = load_subject(subject,gratingpath);
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')

%% Eye data
% Removes times where pupil size is lost
D.eyeLabels(D.eyePos(:,3)==0,:)=4; %track lost
D.eyePos(D.eyePos(:,3)==0,:)=nan;

[EyeStat, Export]=do_eye_analyses_for_export(D); %simpler
%[EyeStat, Export]=do_spike_count_eye_analyses_for_export(D);
%Stat.(subject) = do_spike_count_analyses(D);


%%
% %-- 12/19/2024 01:41:59 PM --%
% cd('/home/huklab/Documents/NPX_pilot/V1Locomotion/Code/')
% clear
% close all
% main_4
% subjects
% subjects0
% subjects
% nsubjs
% subject='npx_V1_calc_eg'
% cd('/home/huklab/Documents/NPX_pilot/V1Locomotion/Code/')
% D = load_subject(subject,gratingpath);
% cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
% nnz(validStim)
% invalid
% nnz(invalid)
% NC
% size(R)
% size(onsets)
% size(spd)
% NT
% range(idOn)
% onsets=StimOnset(validStim);
% size(onsets)
% min(onsets)
% for cc = 1:NC
% %
% fprintf('%d/%d\n', cc, NC)
% cid = cids(cc);
% Stat.cid(cc) = cid;
% unitix = D.spikeIds == cid;
% sessix = unique(D.sessNumSpikes(unitix));
% gtix = find(ismember(D.sessNumGratings, sessix));
% gtix(isnan(D.GratingDirections(gtix))) = [];
% onsets=StimOnset(validStim);
% winsize = median(D.GratingOffsets(gtix) - D.GratingOnsets(gtix));
% t0 = min(onsets) - 2*winsize;
% st = D.spikeTimes(unitix) - t0;
% onsets = onsets - t0;
% st(st < min(onsets)) = [];
% sp = struct('st', st, 'clu', ones(numel(st),1));
% % quick bin spikes during gratings
% R = binNeuronSpikeTimesFast(sp, onsets, winsize);
% R = R ./ winsize;
% Robs(:,cc)=R;
% end
% size(runSpeed)
% size(Robs)
% size(spd)
% size(StimDir)
% size(StimSpeed)
% unique(StimSpeed)
% unique(StimFreq)
% unique(D.GratingFrequency)
% unique(D.GratingSpeeds)
% Export.StimSF=StimFreq;
% Export.StimTF=StimSpeed.*StimFreq;
% clear
% close all
% Exporttesting
% isnan(Export.runSpeed)
% nnz(isnan(Export.runSpeed))
% nnz(isnan(Export.SacHz))
% help var
% nnz(isnan(Export.SacHz))
% nnz(isnan(Export.pupilarea))
% nnz(isnan(Export.varxy))
% nnz(isnan(Export.varXY))
% nnz(isnan(Export.pupilarea))
% nnz(isnan(runSpeed))
% nnz(isnan(runSpeed))/numel(runSpeed)
% size(valid)
% nbins
% nnz(isnan(treadSpeed))/numel(treadSpeed)
% size(isnan(runSpeed))
% size(sum(isnan(runSpeed)))
% (sum(isnan(runSpeed)))
% (sum(isnan(runSpeed)'\))
% (sum(isnan(runSpeed)'))
% (sum(isnan(runSpeed)'))>(nbins/2)
% ((sum(eyeFlag'==4))>(nbins/2))
% nnz((sum(eyeFlag'==4))>(nbins/2))
% nnz(isnan(Export.varXY))
% nnz(isnan(Export.Eyevar))
% save('Gru_V1_Calc_20220412_test','Export');
% load('Gru_V1_Calc_20220412_test.mat')
% clear
% close all

%%