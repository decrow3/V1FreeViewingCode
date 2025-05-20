% Find the running info on the simultaneous deep and shallow V1 recording sessions.

%% Session 14
gratingpath='/media/huklab/Data/NPX/HuklabTreadmill/gratings/';
subject = 'gru';
D = load_subject(subject,gratingpath);

%% 
D.eyeLabels=D.eyeLabels((D.sessNumEye==14));
D.eyePos=D.eyePos((D.sessNumEye==14),:);
D.eyeTime=D.eyeTime((D.sessNumEye==14),:);
D.sessNumEye=D.sessNumEye((D.sessNumEye==14),:);
D.spikeIds=D.spikeIds((D.sessNumSpikes==14 ),:);
D.spikeTimes=D.spikeTimes((D.sessNumSpikes==14 ),:);
D.sessNumSpikes=D.sessNumSpikes((D.sessNumSpikes==14),:);
D.treadSpeed=D.treadSpeed((D.sessNumTread==14),:);
D.treadTime=D.treadTime((D.sessNumTread==14),:);
D.sessNumTread=D.sessNumTread((D.sessNumTread==14),:);
D.GratingContrast=D.GratingContrast((D.sessNumGratings==14),:);
D.GratingDirections=D.GratingDirections((D.sessNumGratings==14),:);
D.GratingFrequency=D.GratingFrequency((D.sessNumGratings==14),:);
D.GratingOffsets=D.GratingOffsets((D.sessNumGratings==14),:);
D.GratingOnsets=D.GratingOnsets((D.sessNumGratings==14),:);
D.GratingSpeeds=D.GratingSpeeds((D.sessNumGratings==14),:);
D.sessNumGratings=D.sessNumGratings((D.sessNumGratings==14),:);

Eyestat.(subject) = do_eye_analyses(D);


nrun = numel(Eyestat.(subject).runTrials);
nstat = numel(Eyestat.(subject).statTrials);
n=min([nrun nstat]);


%% Histograms
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


%% Session 16
D = load_subject(subject,gratingpath);

%% 
D.eyeLabels=D.eyeLabels((D.sessNumEye==16));
D.eyePos=D.eyePos((D.sessNumEye==16),:);
D.eyeTime=D.eyeTime((D.sessNumEye==16),:);
D.sessNumEye=D.sessNumEye((D.sessNumEye==16),:);
D.spikeIds=D.spikeIds((D.sessNumSpikes==16 ),:);
D.spikeTimes=D.spikeTimes((D.sessNumSpikes==16 ),:);
D.sessNumSpikes=D.sessNumSpikes((D.sessNumSpikes==16),:);
D.treadSpeed=D.treadSpeed((D.sessNumTread==16),:);
D.treadTime=D.treadTime((D.sessNumTread==16),:);
D.sessNumTread=D.sessNumTread((D.sessNumTread==16),:);
D.GratingContrast=D.GratingContrast((D.sessNumGratings==16),:);
D.GratingDirections=D.GratingDirections((D.sessNumGratings==16),:);
D.GratingFrequency=D.GratingFrequency((D.sessNumGratings==16),:);
D.GratingOffsets=D.GratingOffsets((D.sessNumGratings==16),:);
D.GratingOnsets=D.GratingOnsets((D.sessNumGratings==16),:);
D.GratingSpeeds=D.GratingSpeeds((D.sessNumGratings==16),:);
D.sessNumGratings=D.sessNumGratings((D.sessNumGratings==16),:);

Eyestat.(subject) = do_eye_analyses(D);


nrun = numel(Eyestat.(subject).runTrials);
nstat = numel(Eyestat.(subject).statTrials);
n=min([nrun nstat]);


%% Histograms
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
