% Combing files for full V2/V1 session
% Attempt 2, january 2023
% Brie V2-V1 from 2022-06-25

%Probe 2, Neuropix-PXI-101.0, was in V1
% Three sessions where recorded
%'/mnt/NPX/Brie/2022-08-09/2022-08-09_17-31-35/Record Node 103/experiment1/recording1/continuous/Neuropix-PXI-101.0/'


%
Path0='/mnt/NPX/';
Subject= 'Brie';
Date        = '06-25';
savepath = ['/home/huklab/Documents/SpikeSorting/Combining/' Subject '_2022-' Date '_Probe2/'];
%savepath= '/home/huklab/Documents/SpikeSorting/Output/Brie_2022-06-25_re/'
startsession=2;
endsession=3;

%disp(({Dates.name}));


Dates= dir(fullfile(Path0, Subject));
Location    = dir(fullfile(Path0, Subject,['*' Date]));




subfolders=dir(fullfile(Location.folder, Location.name));
for ii=1:length(subfolders)
    valid(ii)=subfolders(ii).name(1)~='.' && subfolders(ii).isdir && ~contains(subfolders(ii).name,'Combined') && ~strcmp(subfolders(ii).name,'unused');
end


nSessions=sum(valid);
SessionNames={subfolders(valid).name};

if isempty(endsession)
    endsession=nSessions;
end

for session=startsession:endsession
    SessionNames2{1+session-startsession}=SessionNames{session};
end
SessionNames=SessionNames2
nSessions=length(SessionNames);

 %% Check for subfolders, recording1 recording2 etc
% for session=1:nSessions
%     RecordingNames{session}=dir(fullfile(Location.folder, Location.name, SessionNames{ii}, 'Record Node 103/experiment1/recording*'));
% 
% end

% TODO FIX THIS UP

%%
%Is this always the same? Sometimes there are more recordings
%Hardcoding for now
Probe1subpath = ['Record Node 103/experiment1/recording1/continuous/' ...
    'Neuropix-PXI-101.0/'];
Probe2subpath = ['Record Node 103/experiment1/recording1/continuous/' ...
    'Neuropix-PXI-101.3/'];

%%

% if ~isfolder(fullfile(Location.folder, Location.name, 'Combined','Probe1'))
%     mkdir(fullfile(Location.folder, Location.name,'Combined','Probe1'))
% end
% 
% if ~isfolder(fullfile(Location.folder, Location.name, 'Combined','Probe2'))
%     mkdir(fullfile(Location.folder, Location.name,'Combined','Probe2'))
% end


%savepath =fullfile(Location.folder, Location.name,'Combined','Probe1');
if ~isfolder(savepath)
    mkdir(savepath)
end

numChannels=384;

%% Pull and write to files
for ii=1:nSessions
    if ii==1
        fileID = fopen(fullfile(savepath, 'continuous.dat'),'w');
    else
        fileID = fopen(fullfile(savepath, 'continuous.dat'),'a');
    end

    filepath {ii}       = fullfile(Location.folder, Location.name, SessionNames{ii}, Probe2subpath);
    cont    {ii}        = dir(fullfile(filepath{ii}, 'continuous.dat'));
    sych_timestamps{ii} = dir(fullfile(filepath{ii}, 'synchronized_timestamps.npy'));
    timestamps{ii}      = dir(fullfile(filepath{ii}, 'timestamps.npy'));

    size_cont(ii)               = cont{ii}.bytes;
    size_sych_timestamps(ii)    = sych_timestamps{ii}.bytes;
    size_timestamps(ii)         = timestamps{ii}.bytes;

    m = memmapfile(fullfile(cont{ii}.folder,cont{ii}.name), 'Format', 'int16')
    Data{ii}=reshape(m.Data,length(m.Data)/384,384);
    COUNT = fwrite(fileID,Data{ii},'int16') 
    fclose(fileID)
end


CombinedCont=dir(fullfile(savepath, 'continuous.dat'))
sum(size_cont)-CombinedCont.bytes


%% Notes on sorting, don't run phy template-gui params.py directly
% source Physort/bin/activate
% cd ''
% phy extract-waveforms params.py
% env QTWEBENGINE_CHROMIUM_FLAGS="--single-process" phy template-gui ./params.py


%% Timestamps and sychs
for ii=1:nSessions

    filepath {ii}       = fullfile(Location.folder, Location.name, SessionNames{ii}, Probe2subpath);
    cont    {ii}        = dir(fullfile(filepath{ii}, 'continuous.dat'));
    sych_timestamps{ii} = dir(fullfile(filepath{ii}, 'synchronized_timestamps.npy')); %Just a vector of -1s
    timestamps{ii}      = dir(fullfile(filepath{ii}, 'timestamps.npy')); % Just a vector incrementing on 30,000Hz clock time

    size_cont(ii)               = cont{ii}.bytes;
    size_sych_timestamps(ii)    = sych_timestamps{ii}.bytes;
    size_timestamps(ii)         = timestamps{ii}.bytes;
    
    times{ii}=readNPY(fullfile(timestamps{ii}.folder,timestamps{ii}.name));

    %After concatenation, time between files is destroyed
    if ii==1
        breaks(ii)=0;
    else
        breaks(ii)=times{ii}(1)-times{ii-1}(end);
    end
end

sp_timestamps=cat(1,times{:});
save(fullfile(savepath,'sp_timestamps.mat'),'sp_timestamps','-v7.3')

%% Continuous FPGA timestamps match the sp_timestamps, good sign that the clocks are the same
FPGApath='Record Node 103/experiment1/recording1/continuous/Rhythm_FPGA-102.0/';
for ii=1:nSessions

    filepath {ii}       = fullfile(Location.folder, Location.name, SessionNames{ii}, Probe2subpath);
    cont    {ii}        = dir(fullfile(filepath{ii}, 'continuous.dat'));
    sych_timestamps{ii} = dir(fullfile(filepath{ii}, 'synchronized_timestamps.npy')); %Just a vector of -1s
    timestamps{ii}      = dir(fullfile(filepath{ii}, 'timestamps.npy')); % Just a vector incrementing on 30,000Hz clock time

    size_cont(ii)               = cont{ii}.bytes;
    size_sych_timestamps(ii)    = sych_timestamps{ii}.bytes;
    size_timestamps(ii)         = timestamps{ii}.bytes;
    
    FPGA_ctimes{ii}=readNPY(fullfile(timestamps{ii}.folder,timestamps{ii}.name));
end

FPGA_ctimestamps=cat(1,FPGA_ctimes{:});
save(fullfile(savepath,'FPGAc_sp_timestamps.mat'),'FPGA_ctimestamps','-v7.3')



%% FPGA and strobe timings
FPGApath='Record Node 103/experiment1/recording1/events/Rhythm_FPGA-102.0/TTL_1/';
sync='Record Node 103/experiment1/recording1/sync_messages.txt';

for ii=1:nSessions

    FPGAfilepath {ii}       = fullfile(Location.folder, Location.name, SessionNames{ii}, FPGApath);
    
    timingfileloc           = [FPGAfilepath{ii}  '/'];
    FPGAchannel_states{ii}      = readNPY([timingfileloc 'channel_states.npy']);
    FPGAchannels {ii}           = readNPY([timingfileloc 'channels.npy']);
    FPGAfull_words {ii}         = readNPY([timingfileloc 'full_words.npy']);
    FPGAtimestamps {ii}         = readNPY([timingfileloc 'timestamps.npy']);
    
    FPGAtimestamps {ii}         = FPGAtimestamps {ii}  - breaks(ii);

    synctext{ii} = importdata(fullfile(Location.folder, Location.name, SessionNames{ii}, sync));
    
    
    for jj=1:length(synctext{ii})
        FPGAtextind(jj)=strcmp(synctext{ii}{jj}(12:22),'Rhythm FPGA');
        if length(synctext{ii}{jj})>47
        NPXtextind(jj)=strcmp(synctext{ii}{jj}(12:47),'Neuropix-PXI Id: 101 subProcessor: 0');
        end
    end
    
    %This has the starttime for the FPGA Rhythm on the same(hopefully)
    %30000Hz clock
    FPGAtext=synctext{ii}{FPGAtextind};
    %NPXtext=synctext{ii}{NPXtextind};
    % Grab digits in text file up to the @ character:
    % regexp(FPGAsynctext{ii}{find(FPGAtext)},'\w*@','match')
    FPGAsynctexttime{ii}=str2num(FPGAtext(regexp(FPGAtext,'\w*@'):(regexp(FPGAtext,'@')-1)));
    %NPXsynctexttime{ii}=str2num(NPXtext(regexp(NPXtext,'\w*@'):(regexp(NPXtext,'@')-1)));
end


% These timestamps seem to match the NPX clock!!
FPGAfull_words_all  = cat(1,FPGAfull_words{:});
FPGAchannels_all    = cat(1,FPGAchannels{:});
FPGAtimestamps_all  = cat(1,FPGAtimestamps{:});
FPGAchannel_states_all = cat(1,FPGAchannel_states{:});

%% Sync drift across FPGA and NPX, NO! SHOULD NOT NEED TO
% FPGA_clockstarts = cat(1,FPGAsynctexttime{:});
% NPX_clockstarts = cat(1,NPXsynctexttime{:});
% FPGA2NPX=fit(FPGA_clockstarts,NPX_clockstarts,'poly1');
% 
% %Putting PTB signals recorded on FPGA in terms of NPX clock
% synced_timestamps_all=FPGA2NPX(double(FPGAtimestamps_all));

%% Saving for consistency with previous analysis

ts.timestamps = FPGAtimestamps_all;%synced_timestamps_all;
ts.fullWords = FPGAfull_words_all;
ts.eventVirCh = FPGAchannels_all;
ts.eventVirChStates = FPGAchannel_states_all;
d_ts = double(FPGAtimestamps_all);
d_ts_s = d_ts/30000;
ts.timestamps_s = d_ts_s;
save(fullfile(savepath,'FPGAtimestamps.mat'), 'ts','d_ts_s','-v7.3');

plot(d_ts_s,ts.eventVirChStates,'.')
plot(d_ts_s,ts.fullWords,'.')


