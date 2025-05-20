%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory

 %% Get timing strobes from spikeGLX setup
timingfilename='/mnt/NPX/Timing Tests Treadmill Rig/timing-test-20230806/TestTiming-newfirm1_g0/TestTiming-newfirm1_g0_t0.nidq.bin';
nchannels=4;
samples=[];
dataformat=[];
memorymap=1;
NiDaQout = spGLX_load_binary_file(timingfilename,nchannels,samples,dataformat,memorymap);
ttl=dec2bin(NiDaQout.Data.data(4,:));


datapixx=ttl(:,1:7)=='1';


Event.samplenumber=find(sum(datapixx,2)>0);
Event.ts=(Event.samplenumber)/30000.076569037657;% from nidaq metafile
Event.signal=ttl(Event.samplenumber,1:7);
%
imecSync=ttl(:,8)=='1';


PDanalog=NiDaQout.Data.data(1,:);
PD=PDanalog>2*10^4;
risePD=find(diff(PD)>0);

%% Received signals from imec card for syncing 
spikefilename='/mnt/NPX/Timing Tests Treadmill Rig/timing-test-20230806/TestTiming-newfirm1_g0/TestTiming-newfirm1_g0_imec0/TestTiming-newfirm1_g0_t0.imec0.ap.bin';
nchannels=385;
samples=[];
dataformat=[];
memorymap=1;
Tout = spGLX_load_binary_file(spikefilename,nchannels,samples,dataformat,memorymap);
imecSent=Tout.Data.data(385,:); %this takes forever
%Signals that were sent from imec card
riseRec=find(diff(imecSync)>0); %in samples on nidaq
riseSent=find(diff(imecSent)>40); %in samples on AP
%%
%save('/home/huklab/Documents/SpikeSorting/Output/2023-08-06/Rocky-0806_g0_t0.nidq_imecAp.mat','riseSent','riseRec','imecSent','imecSync')
% load('/mnt/NPX/Timing Tests Treadmill Rig/timing-test-20230806/Rocky-0806_g0_t0.nidq_imecAp.mat')
%%

Event.ClosestAPSamples = interp1(riseRec,riseSent,Event.samplenumber);
Event.Ap_ts=Event.ClosestAPSamples/30000.0; % from AP metafile

%risephD=find(diff(phDiode)>0); %in samples on nidaq
Event.PhDiodeClosestAPSamples = interp1(riseRec,riseSent,risePD);
Event.PD_ts=Event.PhDiodeClosestAPSamples/30000.0; % from AP metafile

%% Events recorded on NiDaq are now in terms of seconds from start of NPX
%recording

%save(['/home/huklab/Documents/SpikeSorting/Output/2023-08-06/' 'Rocky-0806_g0_t0_2.nidq' 'dp_ev_ts.mat'],'Event')
% 
% close all
% clear Tout NiDaQout
%%

%/mnt/NPX/Rocky/20230806/Rocky-v1-0806_g0/Rocky-v1-0806_g0_t0.nidq.bin'
%load(['/home/huklab/Documents/SpikeSorting/Output/2023-08-06/' 'Rocky-0806_g0_t0_2.nidq' 'dp_ev_ts.mat'],'Event')

%%
RecStrobes=(Event.signal(:,7)=='1');
nStrobes=sum(RecStrobes)
t=Event.Ap_ts(RecStrobes);

Recbits=Event.signal(:,1:6);
%Bits 2 and 4 were switched!
Recbits(:,2)=Event.signal(:,4);
Recbits(:,4)=Event.signal(:,2);


RecWords=Recbits(sum(Recbits,2)>0,:);
RecWords=bin2dec(num2str(RecWords));
RecWords_ts=Event.Ap_ts(sum(Recbits,2)>0);
%%
scatter(RecWords_ts,RecWords,'.')
title('PTB sent signals')
%% Load stimuli data from once data was recording, move calibration/setup etc to subfolder
StimPath='/mnt/NPX/Timing Tests Treadmill Rig/timing-test-20230806/';
Exp = io.basic_marmoview_import(StimPath);

%% Check number of timestamps
ptbsw = cell2mat(cellfun(@(x) x.STARTCLOCK(:)', Exp.D, 'uni', 0));
ptbew = cell2mat(cellfun(@(x) x.ENDCLOCK(:)', Exp.D, 'uni', 0));

ptbwords = reshape([ptbsw'; ptbew'], [], 1); % all strobe words

fprintf('PTB sent %d strobed words during this session\n', numel(ptbwords))
fprintf('Ephys recovered %d strobes (flip on channel 3)\n', nStrobes)

%% Prepare to synch using xcorr
ptbst = cellfun(@(x) x.STARTCLOCKTIME, Exp.D); % ptb trial starts
ptbet = cellfun(@(x) x.ENDCLOCKTIME, Exp.D); % ptb trial ends

%% Plot sent signals
figure(1);clf
hold on
scatter(ptbst,ptbsw(:,1))
scatter(ptbst,ptbsw(:,2))
scatter(ptbst,ptbsw(:,3))
scatter(ptbst,ptbsw(:,4))
scatter(ptbst,ptbsw(:,5))
scatter(ptbst,ptbsw(:,6))

title('PTB sent signals')
%% Remove unrecorded times for syncing
%
ptbt = reshape([ptbst'; ptbet'], [], 1); % all strobe times from ptb

xlim([1691103500 1691106500])
% min =1691100000;
% max =Inf;
%%

%%
figure(1); clf
plot(diff(ptbt), 'o')
hold on


%t = ts(find(diff(ts) > 2e-3) + 1);
%figure(2); clf
%T = max([t-t(1);ptbt-ptbt(1)]);
T=max([ptbt-ptbt(1)]);

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

%% Get best lag to match time series
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
bin_size = 1e-3;
bins = 0:bin_size:T;

[ecnt, ~, eid] = histcounts(t-t(1)+offset, bins);
[pcnt, ~, pid] = histcounts(ptbt-ptbt(1), bins);

figure(4); clf; set(gcf, 'Color', 'w')
plot(bins(1:end-1), pcnt); hold on
plot(bins(1:end-1), ecnt, '.')
%xlim([3000, 7000])
legend({'PTB Strobes', 'EPHYS Strobes'})
xlabel('Time from first PTB Strobe (seconds)')
title('Binned Strobe Times')


%%
 potential_matches = find(ecnt==1 & pcnt==1);
%potential_matches = find(ecnt>0 & pcnt>0);
fprintf('found %d potential matches at lag %02.3fs\n', numel(potential_matches), lags(id)*bin_size)

ematches = ismember(eid,potential_matches);
pmatches = ismember(pid, potential_matches);

assert(sum(ematches)==sum(pmatches), 'numer of matching timestamps does not match')
clf

%%

fun = synchtime.align_clocks(ptbt(pmatches), t(ematches));
plot(ptbt(pmatches), t(ematches), 'o'); hold on
plot(xlim, fun(xlim), 'k')
xlabel('PTB clock')
ylabel('Ephys clock')


sprintf('%15.15f',fun(0))

%% Plot check
figure(2)
plot(ptbst+fun(0),ptbsw,'.')
hold on
scatter(RecWords_ts,RecWords,'.')
hold off
%%
x1=-1; x2=12200;
for ll=1:10
figure(2)
plot(ptbst+fun(0),ptbsw,'.')
hold on
xlim([x1 x2])
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2])
pause(0.2)

end
%% Not looking great

%Exp.ptb2Ephys = fun;


%% Run Fixbiter.m
toffset=-1691103782.135725021362305;
sprintf('%15.15f',toffset)
w=[1 -toffset];
@(t)(t-w(2))/w(1)
myfun=@(t)(t-w(2))/w(1);

% %% Manual sync gave offset of -1686257209.140790939331055 for ptb2Ephys
% w=[1 1686257209.14079];
% @(t)(t-w(2))/w(1)
% myfun=@(t)(t-w(2))/w(1);

%% check
x1=-2000; x2=2500;
for ll=1:10
figure(2)
plot(ptbst+myfun(0),ptbsw,'.')
hold on
xlim([x1 x2]); ylim([0 64]);
pause(0.2)
scatter(RecWords_ts,RecWords,'.')
hold off
xlim([x1 x2]); ylim([0 64]);
pause(0.2)

end

%%
Exp.ptb2Ephys = myfun;


%% fake strobe times (necessary to synchronize the eye tracker)
for iTrial = 1:numel(Exp.D)

    t = Exp.ptb2Ephys(Exp.D{iTrial}.STARTCLOCKTIME);
    if t > 0
        Exp.D{iTrial}.START_EPHYS = t;
    end

    t = Exp.ptb2Ephys(Exp.D{iTrial}.ENDCLOCKTIME);
    if t > 0
        Exp.D{iTrial}.END_EPHYS = t;
    end
end

%%
firstchannel=double(Tout.Data.data(1,:));
%%
thresholdedSIG=(firstchannel>0) - (firstchannel<-50);
%% Signals that were sent from imec card
riseRec=find(diff(imecSync)>0); %in samples on nidaq
riseSent=find(diff(imecSent)>40); %in samples on AP
fallSent=find(diff(imecSent)<-40); %in samples on AP
%%
%%
subplot(211)
plot(riseSent'./30000,'.'); hold on
plot(riseRec./30003.003,'.'); hold off
title('RiseTimes(s)')
subplot(212)

%%
plot(riseSent'./30000 - riseRec./30003.003,'.');
title('Difference(s)')
%%
% Event.ClosestAPSamples = interp1(riseRec,riseSent(2:end),Event.samplenumber);
% Event.Ap_ts=Event.ClosestAPSamples/30000; % from AP metafile

pdAnalog=(NiDaQout.Data.data(1,:));
%this is not right, something is wrong with this circuit, if we want -ve to
%mean current flows then +ve needs to be connected to ground



phDiode=pdAnalog<10000;




risephD=find(diff(phDiode)>0); %in samples on nidaq
risephD(find(diff(risephD)<20)+1)=[]; %Remove repeats/flicker
Event.PhDiodeClosestAPSamples = interp1(riseRec,riseSent,risephD);
Event.PD_ts=Event.PhDiodeClosestAPSamples/30000; % from AP metafile

%% These are already on imec card time (AP time)
thresholdedUP=(firstchannel>0);
thresholdedDWN=(firstchannel<-50);
riseSIG=find(diff(thresholdedUP)>0);
fallSIG=find(diff(thresholdedDWN)>0);
%%
Event.riseSIGClosestAPSamples = riseSIG;%interp1(riseRec,riseSent,riseSIG);
Event.SIGup_ts=Event.riseSIGClosestAPSamples/30000; % from AP metafile

Event.fallSIGClosestAPSamples = fallSIG;%interp1(riseRec,riseSent,fallSIG);
Event.SIGdwn_ts=Event.fallSIGClosestAPSamples/30000; % from AP metafile

%%
subplot(111); hold off
plot(Event.PD_ts,1.1*ones(size(Event.PD_ts)),'.'); hold on
title('PD times')
plot(Event.SIGup_ts,ones(size(Event.SIGup_ts)),'.')
plot(Event.SIGdwn_ts,0.9*ones(size(Event.SIGdwn_ts)),'.'); ylim([0.8 1.2])
%for span=1:2:2500
% xlim([span span+500]); ylim([0.5 1.5]);
% pause(0.01)
% end
%%
xlim([190 250])
% xlim([2625 2650])

%% PD_ts measured on the nidaq is ~0.5s before the signal on the probe?
% How and why? Could the headstage calibration be wrong??

subplot(211)
plot(pdAnalog)
subplot(212)
plot(firstchannel)


%% Xcorr all signals
t1=Event.SIGup_ts;
t2=Event.PD_ts;

% %% Upsampling in time?
% bin_size=1/1000;
% 
% T1interp=t1(1):bin_size:t1(end);
% for ii=1:numel(T1interp)-1
% bin1(ii)=sum(t1>T1interp(ii)&t1<T1interp(ii+1));
% %val(ii)=sum(allclockagain(t1>T1interp(ii)&t1<T1interp(ii+1)));
% end
% t2interp=t2(1):bin_size:t2(end);
% for ii=1:numel(t2interp)-1
% bin2(ii)=sum(t2>t2interp(ii)&t2<t2interp(ii+1));
% %Rval(ii)=sum(RecWords(t2>t2interp(ii)&t2<t2interp(ii+1)));
% end
% 
% plot(bin1)
% plot(bin2)

%no looping, NOTE: need to be on same time series!
sampling=1000;
bin_size=1/sampling; 
t0=min(round([t1(1) t2(1)].*sampling));
tend=max(round([t1(end) t2(end)].*sampling));
tbins= t0:tend;
%t1bins=round(t1(1)*1000):round(t1(end)*1000);
rounded1=round(t1*sampling);
[C,IA,IB1] = intersect(rounded1,tbins);
indchain=zeros(size(tbins));
indchain(IB1)=1;
bin1=indchain;

%t2bins=round(t2(1)*1000):round(t2(end)*1000);
rounded2=round(t2*sampling);
[C,IA,IB2] = intersect(rounded2,tbins);
indchain2=zeros(size(tbins));
indchain2(IB2)=1;
bin2=indchain2;
%%
figure(3); clf;
subplot(111);
plot(tbins./sampling,bin2.*1.1,'.');hold on
 plot(tbins./sampling,bin1,'.');hold off

ylim([0.9 1.2])
% xlim([2625 2650])
 xlim([150 200]);
%Still a huge difference?? How is that possible?


%% Get best lag to match time series
%to convert back to time
[out, lags] = xcorr(bin1, bin2, 'None');
[~, id] = max(out);

figure(4); clf; set(gcf, 'Color', 'w')
plot(lags*bin_size, out); hold on
xlabel('Lag (seconds)')
ylabel('Cross Correlation')
plot(lags(id)*bin_size, out(id), 'o')
legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')

offset = lags(id)*bin_size;
fprintf('Best guess for initial offset: %02.2f (s)\n', offset)


%% Close enough?

figure(2); subplot(111)
hold off
plot(t1-t1(1)-offset,ones(size(t1)),'.')
pause(0.5)
hold on
plot(t2-t2(1),1.1*ones(size(t2)),'.'); ylim([.9 1.2]); %xlim([750 800]);

%% Just local peak
middleout=round(length(out)/2 - 40000):round(length(out)/2 + 40000);
plot(lags(middleout)*bin_size, out(middleout));
[~, midid] = max(out(middleout));
midlags=lags(middleout);
midout=out(middleout);
plot(midlags(midid)*bin_size, midout(midid), 'o')
legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')

midoffset=midlags(midid)*bin_size
fprintf('Best guess for initial offset: %02.3f (s)\n', midoffset)

%%
xlim([-20 20])

%% Marmoview sent times
% Flash = cellfun(@(x) x.PR.Flashtime, Exp.D, 'UniformOutput', false); % ptb trial starts
nTrials=length(Exp.D);
FlashTimes=cell(1,nTrials);
figure(1);clf;
hold on
for ii=1:nTrials
    if isfield(Exp.D{ii},'PR')
        if isfield(Exp.D{ii}.PR,'Flashtime')
            FlashTimes{ii}=Exp.D{ii}.PR.Flashtime;
        end
    end
    
end
hold off
plot(cell2mat(FlashTimes)','.');

%%
PTBPDtime=cell2mat(FlashTimes);
PTBPDtimeV=PTBPDtime(:);

SentPDtime=Exp.ptb2Ephys(PTBPDtimeV);


%%
subplot(111); hold off
plot(Event.PD_ts,1.1*ones(size(Event.PD_ts)),'.'); hold on
title('PD times')
plot(Event.SIGup_ts,ones(size(Event.SIGup_ts)),'.')
plot(Event.SIGdwn_ts,0.9*ones(size(Event.SIGdwn_ts)),'.'); ylim([0.8 1.2])
%for span=1:2:2500
% xlim([span span+500]); ylim([0.5 1.5]);
% pause(0.01)
% end

%% Compare marmoview and probe

%% Xcorr all signals
t1=Event.SIGup_ts;
t2=SentPDtime;


%no looping, NOTE: need to be on same time series!
sampling=1000;
bin_size=1/sampling; 
t0=min(round([t1(1) t2(1)].*sampling));
tend=max(round([t1(end) t2(end)].*sampling));
tbins= t0:tend;
%t1bins=round(t1(1)*1000):round(t1(end)*1000);
rounded1=round(t1*sampling);
[C,IA,IB1] = intersect(rounded1,tbins);
indchain=zeros(size(tbins));
indchain(IB1)=1;
bin1=indchain;

%t2bins=round(t2(1)*1000):round(t2(end)*1000);
rounded2=round(t2*sampling);
[C,IA,IB2] = intersect(rounded2,tbins);
indchain2=zeros(size(tbins));
indchain2(IB2)=1;
bin2=indchain2;
%%
figure(3); clf;
subplot(111);
plot(tbins./sampling,bin2.*1.1,'.');hold on
 plot(tbins./sampling,bin1,'.');hold off

ylim([0.9 1.2])
% xlim([2625 2650])
 xlim([150 200]);
%Still a huge difference?? How is that possible?


%% Get best lag to match time series
%to convert back to time
[out, lags] = xcorr(bin1, bin2, 'None');
[~, id] = max(out);

figure(4); clf; set(gcf, 'Color', 'w')
plot(lags*bin_size, out); hold on
xlabel('Lag (seconds)')
ylabel('Cross Correlation')
plot(lags(id)*bin_size, out(id), 'o')
legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')

offset = lags(id)*bin_size;
fprintf('Best guess for initial offset: %02.2f (s)\n', offset)


%% Close enough?

figure(2); subplot(111)
hold off
plot(t1-offset,ones(size(t1)),'.')
pause(0.5)
hold on
plot(t2 ,1.1*ones(size(t2)),'.'); ylim([.9 1.2]); %xlim([750 800]);

%% Just local peak
middleout=round(length(out)/2 - 40000):round(length(out)/2 + 40000);
plot(lags(middleout)*bin_size, out(middleout));
[~, midid] = max(out(middleout));
midlags=lags(middleout);
midout=out(middleout);
plot(midlags(midid)*bin_size, midout(midid), 'o')
legend({'X Corr', 'Best Guess for Offset'}, 'Location', 'Best')

midoffset=midlags(midid)*bin_size
fprintf('Best guess for initial offset: %02.3f (s)\n', midoffset)

%%
xlim([-20 20])