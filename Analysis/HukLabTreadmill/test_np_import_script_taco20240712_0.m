%% Step 0: set your paths
% The FREEVIEWING codebase uses matlab preferences to manage paths (so
% different users can have different paths)
cd('/home/huklab/Documents/NPX_pilot/V1FreeViewingCode')
addFreeViewingPaths('ucla1') % switch to your user
addpath Analysis/HukLabTreadmill/ % the code will always assume you're running from the FreeViewing base directory


%% Try to synchronize using strobed values
% OpenEphys board timing
dataPath = '/mnt/MGS/Ephys/Raw/taco_20240712/Record Node 101/experiment1/recording1/events/Acquisition_Board-114.Rhythm Data/TTL/';
timestampsOE = readNPY(fullfile(dataPath,'timestamps.npy'));
statesOE      =readNPY(fullfile(dataPath,'states.npy'));
full_wordsOE  =readNPY(fullfile(dataPath,'full_words.npy'));
sample_numbersOE =readNPY(fullfile(dataPath,'sample_numbers.npy'));

%% NPX imec timing
dataPath = '/mnt/MGS/Ephys/Raw/taco_20240712/Record Node 101/experiment1/recording1/events/Neuropix-PXI-100.ProbeA-AP/TTL/';
timestampsNPX = readNPY(fullfile(dataPath,'timestamps.npy'));
statesNPX      =readNPY(fullfile(dataPath,'states.npy'));
full_wordsNPX  =readNPY(fullfile(dataPath,'full_words.npy'));
sample_numbersNPX =readNPY(fullfile(dataPath,'sample_numbers.npy'));

%% First sync imec to OE
riseSent=sample_numbersNPX(statesNPX==1);
fallSent=sample_numbersNPX(statesNPX==-1);

%% Analog signals
% %Continous data on OE board saved as 32channels of
% %headstage and 8 channels of ADC for 40 channels x two bytes per channel,
% %80 bytes per sample. 30000samples/sec,2.4e6bytes per sec
% OEdaqdata=memmapfile('/mnt/MGS/Ephys/Raw/taco_20240712/Record Node 101/experiment1/recording1/continuous/Acquisition_Board-114.Rhythm Data/continuous.dat');
% %% test plots of channels to find the 1Hz imec signal
% hold off
% for ii=33:40
%     start=2.4e6*3600; %1hr into recording
%     Analogtest=double(nidaqdata.Data(start+((2*ii-1):80:2.4e7)))*1+double(nidaqdata.Data(start+((2*ii):80:2.4e7)))*2^8;plot(Analogtest)
%     title(num2str(ii))
%     pause(1)
%     hold on
% end
% %not seeing it...
% hold off
%%
riseRec=sample_numbersOE(statesOE==8);
fallRec=sample_numbersOE(statesOE==-8);

plot(riseSent,riseRec,'.')
%% Remove signals from before the start of recording, (only needed if using time in seconds rather than samples)
riseSent=riseSent(riseRec~=-1);
riseRec=riseRec(riseRec~=-1);

fallSent=fallSent(fallRec~=-1);
fallRec=fallRec(fallRec~=-1);

plot(riseSent,riseRec,'.')
%maximum drift 
max(double(abs(riseSent-riseRec)))/30000

%% Convert OE timestamps into matched ephys time
sample_numbersOE_on_NPXclock = interp1(double(riseRec),double(riseSent),double(sample_numbersOE));

%% Putting it into equivalent .kwe format
eventId = double(sign(statesOE)==1);
time_samples = sample_numbersOE_on_NPXclock; %this might be expecting samples not time
event_channels = double(abs(statesOE)-1);

%% remove last channel (7) that was used for 1Hz syncing
eventId=eventId(event_channels~=7);
time_samples=time_samples(event_channels~=7);
event_channels=event_channels(event_channels~=7);

%% save out digital events
%load events only
timestamps = time_samples;%hdf5read(filename, '/event_types/TTL/events/time_samples');
highlow = eventId;%hdf5read(filename, '/event_types/TTL/events/user_data/eventID');
bitNumber = event_channels;%hdf5read(filename, '/event_types/TTL/events/user_data/event_channels');

%no longer have a strobe signal, so need to fake it
strobeSet=find(highlow==1); %any bit high
strobeUnset=find(highlow==0); %any bit goes low
strobeUnset=[1; strobeUnset];

value=nan(size(strobeSet));
for iStrobe=1:length(strobeSet)
     ts=timestamps <= timestamps(strobeSet(iStrobe)) & timestamps >= timestamps(strobeUnset(iStrobe)) & bitNumber~=7;
     value(iStrobe)=sum(2.^bitNumber(ts) .* highlow(ts));
end

eventTimes=double(timestamps(strobeSet));
eventValues = value;
flagBits = nan(size(value));
flagData = value;
invertedBits = false;
switchbits = false;

%eventSamples depends on recording number.
location='/recordings';
info=h5info([inDir filesep name],location);  %where file was filename=[inDir filesep name];
nRecs=length(info.Groups);

st_index        = double(sample_numbersNPX(1));%strcmp('start_time',{info.Groups(iRec).Attributes.Name});
recStartTime    = double(timestampsNPX(1));%double(info.Groups(iRec).Attributes(st_index).Value);
eventSamples    = eventTimes-recStartTime;

%save([outDir filesep 'events_' exptDate '.mat'], 'eventTimes', 'eventSamples', 'eventValues', 'flagBits', 'flagData', 'invertedBits', 'switchbits','recStartTime');



%% Loadup PTB signals sent through datapixx
load('/home/huklab/Documents/MGS/pdsdata')

PTBword=cellfun(@(X) X.unique_number,PDS.data,'UniformOutput',false)
for ii=1:length(PTBword)
PTBword{ii}(1)=PTBword{ii}(1)-2000;
end
ptbwords=cell2mat(PTBword);
ptbt=cell2mat(PTBts')';

