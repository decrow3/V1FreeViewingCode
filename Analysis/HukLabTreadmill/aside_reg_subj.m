%% Collect up Rho and Pvals
corrSacFRRho= []; corrSacFRPval= []; 
corrSacMagFRRho= []; corrSacmagFRPval= []; 
corrSizeFRRho= []; corrSizeFRPval= []; 
corrVarXFRRho= []; corrVarXFRPval= []; 
corrVarYFRRho= []; corrVarYFRPval= []; 
corrVarXYFRRho= []; corrVarXYFRPval= [];

corrSacFRRhoSt= []; corrSacFRPvalSt= []; 
corrSacMagFRRhoSt= []; corrSacmagFRPvalSt= []; 
corrSizeFRRhoSt= []; corrSizeFRPvalSt= []; 
corrVarXFRRhoSt= []; corrVarXFRPvalSt= []; 
corrVarYFRRhoSt= []; corrVarYFRPvalSt= []; 
corrVarXYFRRhoSt= []; corrVarXYFRPvalSt= [];

corrSacFRRhoB= []; corrSacFRPvalB= []; 
corrSacMagFRRhoB= []; corrSacmagFRPvalB= []; 
corrSizeFRRhoB= []; corrSizeFRPvalB= []; 
corrVarXFRRhoB= []; corrVarXFRPvalB= []; 
corrVarYFRRhoB= []; corrVarYFRPvalB= []; 
corrVarXYFRRhoB= []; corrVarXYFRPvalB= []; 

corrSacRho= []; corrSacPval= []; 
corrSizeRho= []; corrSizePval= []; 
corrVarXRho= []; corrVarXPval= []; 
corrVarYRho= []; corrVarYPval= []; 
corrVarXYRho= [];corrVarXYPval= [];

corrSacMagRhoB= []; corrSacmagPvalB= []; 
corrSacMagRhoSt= []; corrSacMagPvalSt= []; 
corrSacMagRho= []; corrSacmagPval= []; 

regRnSac_B= [];  regRnSac_STATS= []; 
regRnSpd_B= [];  regRnSpd_STATS= []; 
regRnSacSpd_B= []; regRnSacSpd_STATS= [];

runrho=[];murun=[];mustat=[];
SacHzRall=[];SacHzSall=[];runningspeed=[];

speedR=[];
speedS=[];

for isubj = 1:nsubjs
    %
    subject = subjects{isubj};
    good = ~isnan(Stat.(subject).meanrate); %~(frStimS(:,2) < 1 | frStimR(:,2) < 1);
    runrho=[runrho; Stat.(subject).runrho(good)];
    murun=[murun; Stat.(subject).murun(good)];
    mustat=[mustat; Stat.(subject).mustat(good)];

    corrSacFRRho= [corrSacFRRho; Stat.(subject).corrSacFRRho(good)]; corrSacFRPval= [corrSacFRPval ; Stat.(subject).corrSacFRPval(good)]; 
    corrSacMagFRRho= [corrSacMagFRRho ; Stat.(subject).corrSacMagFRRho(good)]; corrSacmagFRPval= [corrSacmagFRPval ; Stat.(subject).corrSacmagFRPval(good)]; 
    corrSizeFRRho= [corrSizeFRRho ; Stat.(subject).corrSizeFRRho(good)]; corrSizeFRPval= [corrSizeFRPval ; Stat.(subject).corrSizeFRPval(good)]; 
    corrVarXFRRho= [corrVarXFRRho ; Stat.(subject).corrVarXFRRho(good)]; corrVarXFRPval= [corrVarXFRRho ; Stat.(subject).corrVarXFRRho(good)]; 
    corrVarYFRRho= [corrVarYFRRho ; Stat.(subject).corrVarYFRRho(good)]; corrVarYFRPval= [corrVarYFRPval ; Stat.(subject).corrVarYFRPval(good)]; 
    corrVarXYFRRho= [corrVarXYFRRho ; Stat.(subject).corrVarXYFRRho(good)]; corrVarXYFRPval= [corrVarXYFRPval ; Stat.(subject).corrVarXYFRPval(good)];
    
    corrSacFRRhoSt= [corrSacFRRhoSt ; Stat.(subject).corrSacFRRhoSt(good)]; corrSacFRPvalSt= [corrSacFRPvalSt ; Stat.(subject).corrSacFRPvalSt(good)]; 
    corrSacMagFRRhoSt= [corrSacMagFRRhoSt ; Stat.(subject).corrSacMagFRRhoSt(good)]; corrSacmagFRPvalSt= [corrSacmagFRPvalSt ; Stat.(subject).corrSacmagFRPvalSt(good)]; 
    corrSizeFRRhoSt= [corrSizeFRRhoSt ; Stat.(subject).corrSizeFRRhoSt(good)]; corrSizeFRPvalSt= [corrSizeFRPvalSt ; Stat.(subject).corrSizeFRPvalSt(good)]; 
    corrVarXFRRhoSt= [corrVarXFRRhoSt ; Stat.(subject).corrVarXFRRhoSt(good)]; corrVarXFRPvalSt= [corrVarXFRPvalSt ; Stat.(subject).corrVarXFRPvalSt(good)]; 
    corrVarYFRRhoSt= [corrVarYFRRhoSt ; Stat.(subject).corrVarYFRRhoSt(good)]; corrVarYFRPvalSt= [corrVarYFRPvalSt ; Stat.(subject).corrVarYFRPvalSt(good)]; 
    corrVarXYFRRhoSt= [corrVarXYFRRhoSt ; Stat.(subject).corrVarXYFRRhoSt(good)]; corrVarXYFRPvalSt= [corrVarXYFRPvalSt ; Stat.(subject).corrVarXYFRPvalSt(good)];
    
    corrSacFRRhoB= [corrSacFRRhoB ; Stat.(subject).corrSacFRRhoB(good)]; corrSacFRPvalB= [corrSacFRPvalB ; Stat.(subject). corrSacFRPvalB(good)]; 
    corrSacMagFRRhoB= [corrSacMagFRRhoB ; Stat.(subject).corrSacMagFRRhoB(good)]; corrSacmagFRPvalB= [corrSacmagFRPvalB ; Stat.(subject).corrSacmagFRPvalB(good)]; 
    corrSizeFRRhoB= [corrSizeFRRhoB ; Stat.(subject).corrSizeFRRhoB(good)]; corrSizeFRPvalB= [corrSizeFRPvalB ; Stat.(subject).corrSizeFRPvalB(good)]; 
    corrVarXFRRhoB= [corrVarXFRRhoB ; Stat.(subject).corrVarXFRRhoB(good)]; corrVarXFRPvalB= [corrVarXFRPvalB ; Stat.(subject).corrVarXFRPvalB(good)]; 
    corrVarYFRRhoB= [corrVarYFRRhoB ; Stat.(subject).corrVarYFRRhoB(good)]; corrVarYFRPvalB= [corrVarYFRPvalB ; Stat.(subject).corrVarYFRPvalB(good)]; 
    corrVarXYFRRhoB= [corrVarXYFRRhoB ; Stat.(subject).corrVarXYFRRhoB(good)]; corrVarXYFRPvalB= [corrVarXYFRPvalB ; Stat.(subject).corrVarXYFRPvalB(good)]; 
    
    corrSacRho= [corrSacRho ; Stat.(subject).corrSacRho(good)]; corrSacPval= [corrSacPval ; Stat.(subject).corrSacPval(good)]; 
    corrSizeRho= [corrSizeRho ; Stat.(subject).corrSizeRho(good)]; corrSizePval= [corrSizePval ; Stat.(subject).corrSizePval(good)]; 
    corrVarXRho= [corrVarXRho ; Stat.(subject).corrVarXRho(good)]; corrVarXPval= [corrVarXPval ; Stat.(subject).corrVarXPval(good)]; 
    corrVarYRho= [corrVarYRho ; Stat.(subject).corrVarYRho(good)]; corrVarYPval= [corrVarYPval ; Stat.(subject).corrVarYPval(good)]; 
    corrVarXYRho= [corrVarXYRho ; Stat.(subject).corrVarXYRho(good)];corrVarXYPval= [corrVarXYPval ; Stat.(subject).corrVarXYPval(good)];
    
    corrSacMagRhoB= [corrSacMagRhoB ; Stat.(subject).corrSacMagRhoB(good)]; corrSacmagPvalB= [corrSacmagPvalB ; Stat.(subject).corrSacmagPvalB(good)]; 
    corrSacMagRhoSt= [corrSacMagRhoSt ; Stat.(subject).corrSacMagRhoSt(good)]; corrSacMagPvalSt= [corrSacMagPvalSt ; Stat.(subject).corrSacMagPvalSt(good)]; 
    corrSacMagRho= [corrSacMagRho ; Stat.(subject).corrSacMagRho(good)]; corrSacmagPval= [corrSacmagPval ; Stat.(subject).corrSacmagPval(good)]; 
    
    regRnSac_B= [regRnSac_B ; Stat.(subject).regRnSac_B(good,:)];  
% regRnSac_STATS= [regRnSac_STATS ; Stat.(subject).regRnSac_STATS(good,:)]; 
     regRnSpd_B= [regRnSpd_B ; Stat.(subject).regRnSpd_B(good,:)]; % regRnSpd_STATS= [regRnSpd_STATS ; Stat.(subject).regRnSpd_STATS(good,:)]; 
%     regRnSacSpd_B= [regRnSacSpd_B ; Stat.(subject).regRnSacSpd_B(good,:)]; regRnSacSpd_STATS= [regRnSacSpd_STATS ; Stat.(subject).regRnSacSpd_STATS(good,:)];


    SacHzRall=[SacHzRall; Stat.(subject).SacHzR(good,:)];
    SacHzSall=[SacHzSall; Stat.(subject).SacHzS(good,:)];
    runningspeed=[runningspeed; Stat.(subject).runningspeed(good,:)];
    
    speedR=[speedR; Stat.(subject).speedR(good,:)];
    speedS=[speedS; Stat.(subject).speedS(good,:)];





%% Baseline firing
ntrials=numel(corrSacFRPvalB);
figure(95);clf
subplot(2,2,1);
hold off
histogram(corrSacFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade rate corr with baseline firing')
pc=nnz(corrSacFRPvalB>0.05)/numel(corrSacFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with saccade rate \n", pc, nnz(corrSacFRPvalB>0.05), ntrials)

%

subplot(2,2,2);
hold off
histogram(corrSacMagFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade mag corr with baseline firing')
pc=nnz(corrSacmagFRPvalB>0.05)/numel(corrSacmagFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with saccade magnitude \n", pc, nnz(corrSacmagFRPvalB>0.05), ntrials)


subplot(2,2,3); hold off
histogram(corrSizeFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Pupil size corr with baseline firing')
pc=nnz(corrSizeFRPvalB>0.05)/numel(corrSizeFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with pupil size \n", pc, nnz(corrSizeFRPvalB>0.05), ntrials)

%

subplot(2,2,4); hold off
histogram(corrVarXYFRRhoB,[-0.225:0.05:.225],'Normalization','probability'  )
title('Eyevar corr with baseline firing')
pc=nnz(corrVarXYFRPvalB>0.05)/numel(corrVarXYFRPvalB)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had baseline firing rates that significantly correlated with eye position variance \n", pc, nnz(corrVarXYFRPvalB>0.05), ntrials)

%% Stimuli driven
figure(96);clf
subplot(2,2,1);
hold off
histogram(corrSacFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade rate corr with stim driven firing')
pc=nnz(corrSacFRPvalSt>0.05)/numel(corrSacFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with saccade rate \n", pc, nnz(corrSacFRPvalSt>0.05), ntrials)


subplot(2,2,2);
hold off
histogram(corrSacMagFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Sacade mag corr with stim driven firing')
pc=nnz(corrSacmagFRPvalSt>0.05)/numel(corrSacmagFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with saccade magnitude \n", pc, nnz(corrSacmagFRPvalSt>0.05), ntrials)

subplot(2,2,3); hold off
histogram(corrSizeFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Pupil size corr with stim driven firing')
pc=nnz(corrSizeFRPvalSt>0.05)/numel(corrSizeFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with pupil size \n", pc, nnz(corrSizeFRPvalSt>0.05), ntrials)

%
subplot(2,2,4); hold off
histogram(corrVarXYFRRhoSt,[-0.225:0.05:.225],'Normalization','probability'  )
title('Eyevar corr with stim driven firing')
pc=nnz(corrVarXYFRPvalSt>0.05)/numel(corrVarXYFRPvalSt)*100;
fprintf("%02.2f%% (%02.0f/ %02.0f) cells had stimulus driven firing rates that significantly correlated with eye position variance \n", pc, nnz(corrVarXYFRPvalSt>0.05), ntrials)


%%
fprintf("-----------COMBINED DATA-----------\n")

fprintf("Running speed and firing rate: %.4f, r^2=%.4f \n\n",nanmean(runrho),nanmean(runrho.^2))

fprintf("Saccade rate and firing rate: %.4f, r^2=%.4f \n",(nanmean(corrSacFRRho)),nanmean(corrSacFRRho.^2))
fprintf("Saccade rate and base firing rate: %.4f, r^2=%.4f \n",(nanmean(corrSacFRRhoB)),nanmean(corrSacFRRhoB.^2))
fprintf("Saccade rate and stim-driven firing rate: %.4f, r^2=%.4f \n\n",(nanmean(corrSacFRRhoSt)),nanmean(corrSacFRRhoSt.^2))

fprintf("Pupil size and firing rate: %.4f, r^2=%.4f \n",(nanmean(corrSizeFRRho)),(nanmean(corrSizeFRRho))^2)
fprintf("Pupil size and base firing rate: %.4f, r^2=%.4f \n",(nanmean(corrSizeFRRhoB)),(nanmean(corrSizeFRRhoB))^2)
fprintf("Pupil size and stim-driven firing rate: %.4f, r^2=%.4f \n\n",(nanmean(corrSizeFRRhoSt)),nanmean(corrSizeFRRhoSt.^2))

fprintf("Saccade magnitude and firing rate: %.4f, r^2=%.4f \n",(nanmean(corrSacMagFRRho)),(nanmean(corrSacMagFRRho))^2)
fprintf("Saccade magnitude and base firing rate: %.4f, r^2=%.4f \n",(nanmean(corrSacMagFRRhoB)),(nanmean(corrSacMagFRRhoB))^2)
fprintf("Saccade magnitude and stim-driven firing rate: %.4f, r^2=%.4f \n\n",(nanmean(corrSacMagFRRhoSt)),nanmean(corrSacMagFRRhoSt.^2))



% fprintf("",(nanmean(corrSacMagFRRhoSt))
% fprintf("",(nanmean(corrSacMagFRRhoSt))^2

% Careful about outliers
baseR=frBaseRall(:,2);
baseR((baseR)>nanmean(baseR)+3*nanstd(baseR))=nan;
baseR((baseR)<nanmean(baseR)-3*nanstd(baseR))=nan;
baseS=frBaseSall(:,2);
baseS((baseS)>nanmean(baseS)+3*nanstd(baseS))=nan;
baseS((baseS)<nanmean(baseS)-3*nanstd(baseS))=nan;
baseratio=(baseR./baseS);
baseratio=baseratio(~isnan(baseratio));
baseratio=baseratio(~isnan(baseratio));
baseratio=baseratio(baseratio~=0);

Basegm=geomean(baseratio,'omitnan' );
fprintf("Geomean (running/stationary) during baseline: %.4f\n" ,Basegm)
stimratio=frStimRall(:,2)./frStimSall(:,2);
stimratio=stimratio(stimratio~=0);
Stimgm=geomean(stimratio,'omitnan' );
fprintf("Geomean (running/stationary) during stimuli: %.4f\n" ,Stimgm)
prefratio=murun(:,1)./(mustat(:,1));
prefratio=prefratio(prefratio~=0);
prefdirgm=geomean(prefratio,'omitnan' );
fprintf("Geomean at preferred direction (running/stationary): %.4f\n" ,prefdirgm)

%%

SacHzDelta=(SacHzRall-SacHzSall); %Per cell diff in saccades
ExpectedFRDelta=SacHzDelta.*regRnSac_B(:,1); %Per cell expected change in firing rate

% FRDeltaStim=frStimR(:,2)-frStimS(:,2);
% fprintf("Mean delta (running-stationary) observed during stimuli: %.4f Hz \n" ,nanmean(FRDeltaStim))
% fprintf("Expected mean delta(running-stationary) from reg. on nSaccades: %.4f Hz \n" ,nanmean(ExpectedFRDelta(:,2)))
% %fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-ExpectedFRDelta(:,2)))
% 
% clear meanspd
% for ii=1:size(runningspeed)
% if ~isempty(runningspeed{ii})
% meanspd(ii)=nanmean(runningspeed{ii});
% else
% meanspd(ii)=nan;
% end
% end
% 
% % %Mean amount of running on trial * B for running effect
% % ExpectedFRRunDelta=meanspd'.*regRnSpd_B(:,1);
% % 
% % fprintf("Expected mean delta(running-stationary) from reg. on run speed: %.4f Hz \n" ,nanmean(ExpectedFRRunDelta))
% % %fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-ExpectedFRRunDelta))
% % 
% % 
% % %fprintf("Expected running rate from geomean: mu*statrate %.4f Hz \n" ,nanmean(Stimgm*frStimS(:,2)))
% % meandelta=(frStimS(:,2).*nanmean(stimratio) - frStimS(:,2));
% % fprintf("Expected mean delta(running-stationary) from mean running effect: %.4f  Hz\n", nanmean(meandelta))
% % %fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-meandelta))
% % fprintf("(Mean only uses one number for  cells) \n \n")
% % % FRDeltaStim./frStimS;
% %  nanmean(FRDeltaStim./frStimS(:,2))

%%

% FRDeltaStim=frStimR(:,2)-frStimS(:,2);
% meangain=mean(FRDeltaStim./frStimS(:,2),'omitnan');
% fprintf("Mean gain (running-stationary)/stationary observed during stimuli: %.4f  \n" ,meangain)
% 
% SacHzDelta=(SacHzR-SacHzS); %Per cell diff in saccades
% ExpectedFRDelta=SacHzDelta.*regRnSac_B(:,1); %Per cell expected change in firing rate
% meangainsac=mean(ExpectedFRDelta(:,2)./frStimS(:,2),'omitnan');
% fprintf("Expected gain (running-stationary)/stationary from reg. on nSaccades: %.4f  \n" ,meangainsac)
% 
% meangainrun=mean((ExpectedFRRunDelta)./frStimS(:,2),'omitnan');
% fprintf("Expected gain (running-stationary)/stationary from reg. on run speed: %.4f  \n" ,meangainrun)
% 
% 
% meandelta=(frStimS(:,2).*nanmean(stimratio) - frStimS(:,2));
% meanmeangain=mean((meandelta)./frStimS(:,2),'omitnan');
% fprintf("Expected gain (running-stationary)/stationary from mean running effect: %.4f  \n", meanmeangain)
% %fprintf("Difference in deltas: %.4f Hz \n\n", nanmean(FRDeltaStim-meandelta))
% fprintf("(Mean only uses one number for  cells) \n \n")
% % FRDeltaStim./frStimS;
%% Delta redon
FRDeltaStim=frStimR(:,2)-frStimS(:,2);
meangain=mean(frStimR(:,2)./frStimS(:,2),'omitnan');
fprintf("Mean delta (running-stationary) observed during stimuli: %.4f Hz\n" ,nanmean(FRDeltaStim))

SacHzDelta=(SacHzRall-SacHzSall); %Per cell diff in saccades
ExpectedFRDelta=SacHzDelta.*regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRSacHzR=SacHzRall.*regRnSac_B(:,1)+regRnSac_B(:,2);
ExpectedFRSacHzS=SacHzSall.*regRnSac_B(:,1)+regRnSac_B(:,2);
ExpectedFRDelta2=ExpectedFRSacHzR-ExpectedFRSacHzS;
fprintf("Expected delta (running-stationary) from reg. on nSaccades: %.4f Hz \n" ,nanmean(ExpectedFRDelta2(:,2)))


speedDelta=(speedR-speedS); %Per cell diff in saccades
ExpectedFRspeedDelta=speedDelta.*regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRspeedR=speedR.*regRnSpd_B(:,1)+regRnSpd_B(:,2);
ExpectedFRspeedS=speedS.*regRnSpd_B(:,1)+regRnSpd_B(:,2);
meangainsac=mean(ExpectedFRspeedR(:,2)./ExpectedFRspeedS(:,2),'omitnan');
fprintf("Expected delta (running-stationary) from reg. on run speed: %.4f Hz \n" ,nanmean(ExpectedFRspeedDelta(:,2)))

%
FRDeltaStim=frStimR(:,2)-frStimS(:,2);
meangain=mean(frStimR(:,2)./frStimS(:,2),'omitnan');
fprintf("Mean gain (running/stationary) observed during stimuli: %.4f  \n" ,meangain)

%SacHzDelta=(SacHzR-SacHzS); %Per cell diff in saccades
ExpectedFRDelta=SacHzDelta.*regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRSacHzR=SacHzRall.*regRnSac_B(:,1)+regRnSac_B(:,2);
ExpectedFRSacHzS=SacHzSall.*regRnSac_B(:,1)+regRnSac_B(:,2);
meangainsac=mean(ExpectedFRSacHzR(:,2)./ExpectedFRSacHzS(:,2),'omitnan');
fprintf("Expected gain (running/stationary) from reg. on nSaccades: %.4f  \n" ,meangainsac)


speedDelta=(speedR-speedS); %Per cell diff in saccades
ExpectedFRDelta=speedDelta.*regRnSac_B(:,1); %Per cell expected change in firing rate
ExpectedFRspeedR=speedR.*regRnSpd_B(:,1)+regRnSpd_B(:,2);
ExpectedFRspeedS=speedS.*regRnSpd_B(:,1)+regRnSpd_B(:,2);
meangainrun=mean(ExpectedFRspeedR(:,2)./ExpectedFRspeedS(:,2),'omitnan');
fprintf("Expected gain (running/stationary) from reg. on run speed: %.4f  \n" ,meangainrun)

%% geomean
geomeangain=geomean(frStimR(:,2)./frStimS(:,2),'omitnan');
fprintf("Geomean gain (running/stationary) observed during stimuli: %.4f  \n" ,geomeangain)
geomeangainsac=geomean(ExpectedFRSacHzR(:,2)./ExpectedFRSacHzS(:,2),'omitnan');
fprintf("Expected gain (running/stationary) from reg. on nSaccades: %.4f  \n" ,geomeangainsac)
geomeangainrun=geomean(ExpectedFRspeedR(:,2)./ExpectedFRspeedS(:,2),'omitnan');
fprintf("Expected gain (running/stationary) from reg. on run speed: %.4f  \n" ,geomeangainrun)
end