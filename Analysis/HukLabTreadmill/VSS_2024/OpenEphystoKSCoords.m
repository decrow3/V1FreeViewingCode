Path0='/mnt/NPX/';
Subject= 'Brie';
Dates= dir(fullfile(Path0, Subject));


%disp(({Dates.name}));

Date        = '2022-06-25';
Location    = dir(fullfile(Path0, Subject,['*' Date]));


subfolders=dir(fullfile(Location.folder, Location.name));
for ii=1:length(subfolders)
    valid(ii)=subfolders(ii).name(1)~='.' && subfolders(ii).isdir && ~strcmp(subfolders(ii).name,'Combined') && ~strcmp(subfolders(ii).name,'unused');
end


nSessions=sum(valid);
SessionNames={subfolders(valid).name};

%% Probe2
session=2;
Settingssubpath = ['Record Node 103/settings.xml'];
filepath        = fullfile(Location.folder, Location.name, SessionNames{session}, Settingssubpath);
% OpenEphys to KSCoords
A=readtable(filepath);

for jj=1:384
ii=jj-1;
try
Xloc{jj}=A.(['CH' num2str(ii) '_1Attribute']);
Yloc{jj}=A.(['CH' num2str(ii) '_2Attribute']);
xpos_probe1(jj)=Xloc{jj}(1,1);
xpos_probe2(jj)=Xloc{jj}(1,2);
ypos_probe1(jj)=Yloc{jj}(1,1);
ypos_probe2(jj)=Yloc{jj}(1,2);
catch
Yloc{jj}=-1; %CH191, the 192nd channel is a ground
xpos_probe1(jj)=-1;
ypos_probe1(jj)=-1;
xpos_probe2(jj)=-1;
ypos_probe2(jj)=-1;

end
end

%%
name=[Date '_' Subject '_Probe1'];
fname=['/home/huklab/Documents/MATLAB/' name];

shankInd=zeros(size(ypos_probe1));
nchan=384;
xcoords=xpos_probe1';
ycoords=ypos_probe1';
% xcoords=xpos_probe2';
% ycoords=ypos_probe2';
connected=ones(size(ypos_probe1));
connected(192)=0;
connected=ypos_probe1>0;
%Copying for SGLXMetaToCoords
        newName = [fname,'_kilosortChanMap.mat'];
        
        chanMap = (1:nchan)';
        chanMap0ind = chanMap - 1;
        connected = logical(connected)';
        %xcoords = shankInd*shankPitch + xCoord;   %KS2 not yet using kcoords, so x coord includes shank sep
        %ycoords = yCoord; % variable names need to  match KS standard
        kcoords = shankInd' + 1;     %KS1 uses kcoords to force templates to be on one shank
        %name = fname;
        save( newName, 'chanMap', 'chanMap0ind', 'connected', 'name', 'xcoords', 'ycoords', 'kcoords' );

          scatter(xpos_probe1,ypos_probe1)
%         scatter(xpos_probe2,ypos_probe2)