function varargout = dataFactoryTreadmill(sessionId, varargin)
% DATAFACTORY is a big switch statement. You pass in the ID of the sesion
% you want and it returns the processed data.
%
% Input:
%   Session [Number or ID string]
% Output: 
%   Exp [struct]
%
% Example:
%   Exp = io.dataFactory(5); % load the 5th session
ip = inputParser();
ip.addParameter('spike_sorting', []) % pass in a spike sorting argument
ip.addParameter('abort_if_missing', false)
ip.parse(varargin{:});

dataPath = getpref('FREEVIEWING', 'HUKLAB_DATASHARE');
addpath(fullfile(dataPath, 'PLDAPStools'))

meta_file = fullfile(dataPath, 'datasets.xls');

data = readtable(meta_file);
nSessions = size(data,1);

if nargin < 1
    sessionId = [];
end

if nargin >=1 && iscell(sessionId) % assume that inputs are condition combinations
    
    argConds = sessionId; % AND conditioned arguments (as pairs)
    
    [~, idx] = io.get_experiments_and(data, argConds{:});
    idx = find(idx);
    for i = idx(:)'
        fprintf('%d) %s\n', i, data.Tag{i})
    end
    
    varargout{1} = data.Tag(idx);
    return
end

if ischar(sessionId)
    sessionId = find(strcmp(data.Tag, sessionId));
end

if isempty(sessionId)
    fprintf('No valid session id passed in. You can call dataFactory with an id number or string:\n')
    for i = 1:nSessions
        fprintf('%d) %s\n', i, data.Tag{i})
    end
    varargout{1} = data.Tag(1:nSessions)';
    return
end
    
if isnumeric(sessionId)
    thisSession = data(sessionId,:);
end
    
field_list = fields(data);
if any(contains(field_list, 'CalibMat'))
    S.CalibMat = [data.CalibMat_1(sessionId), ...
        data.CalibMat_2(sessionId), ...
        data.CalibMat_3(sessionId), ...
        data.CalibMat_4(sessionId), ...
        data.CalibMat_5(sessionId)];
else
    S.CalibMat = [];
end

S.Latency = 8.3e-3; % delay from PTB time to actual monitor refresh (even at 240Hz there's about an 8ms latency)
S.rect = [-20 -60 50 10]; % default gaze-centered ROI (pixels)

S.processedFileName = [thisSession.Tag{1} '.mat'];

rootDir = strrep(thisSession.RootDirectory{1}, '/', filesep);
rawDir = strrep(thisSession.Directory{1}, '/', filesep);
serverDir = getpref('FREEVIEWING', 'HUKLAB_DATASHARE');

S.rootDir = rootDir;
S.rawDir = rawDir;
S.serverDir = serverDir;
S.rawFilePath = fullfile(serverDir, rootDir, rawDir);
S.spikeSorting = ip.Results.spike_sorting;

% try loading the file
fname = fullfile(dataPath, S.processedFileName);
if exist(fname, 'file')
    fprintf('Loading [%s]\n', S.processedFileName)
    Exp = load(fname);
    
    if ~isempty(ip.Results.spike_sorting)
        spfname = fullfile(dataPath, 'spikes', sprintf('%s_%s.mat', strrep(Exp.FileTag, '.mat', ''), ip.Results.spike_sorting));
        if exist(spfname, 'file')
            sp = load(spfname);
            Exp.osp = sp;
        else
            error('dataFactory: requested spike sorting does not exist')
        end
        if contains(dataPath, 'Gabe') %temporary for Gabe
            Exp.osp = sp.sp;
        end
    end
    
    fprintf('Done\n')
elseif ip.Results.abort_if_missing
    Exp = [];
    
else % try importing the file
       
    S.importFun = str2func(['io.' thisSession.ImportFun{1}]);
    fprintf('Could not find [%s]\n', fname)
    fprintf('Trying to import the data from [%s]\n', S.rawFilePath)
    
    % incase the import function changes to match changes in the
    % dataformat, this is a function handle that can be specific to certain
    % functions
    Exp = S.importFun(S); 
    
    % --- get spikes from concatenated session file
    fpath = fileparts(S.rawFilePath);
    subj = S.processedFileName(1:strfind(S.processedFileName, '_')-1);
    spikesfname = fullfile(fpath, [subj '_All_cat.mat']);
    fprintf('Loading spikes from [%s]\n', spikesfname)
    D = load(spikesfname);
    
    id = strfind(S.processedFileName,'_');
    subj = S.processedFileName(1:id);
    dat = datenum(S.processedFileName(id+(1:8)), 'yyyymmdd');
    dateStr = [subj datestr(dat, 'yyyy-mm-dd')];
    st = [];
    clu = [];
    
    exname = strrep(S.processedFileName, '.mat', '');
    figDir = fullfile(dataPath, 'imported_sessions_qa', exname);
    if ~exist(figDir, 'dir')
        mkdir(figDir)
    end
    
    recId = find(contains(D.z.RecId, dateStr));
    for rId = recId(:)'
        
        figure(rId); clf
      
        unitlist = find(cellfun(@(x) ~isempty(x), D.z.Times{rId}));
        nUnits = numel(unitlist);
        sx = ceil(sqrt(nUnits));
        sy = round(sqrt(nUnits));
        ax = plot.tight_subplot(sx, sy, 0.05);
        
        for iunit = 1:nUnits
            kunit = unitlist(iunit);
            
            set(gcf, 'currentaxes', ax(iunit))
            plot(squeeze(D.z.Shapes(:,kunit,:))', 'Linewidth', 2)
            title(sprintf('Unit: %d', kunit))
            
            stmp = double(D.z.Times{rId}{kunit}) / D.z.Sampling;
            st = [st; stmp];
            clu = [clu; kunit*ones(numel(stmp),1)];
        end
        
    end
    
    plot.suplabel('Waveforms', 't');
    plot.fixfigure(gcf, 10, [sy sx]*2)
    saveas(gcf, fullfile(figDir, 'waveforms.pdf'))
    
    
    
    Exp.osp = struct('st', st, 'clu', clu, 'cids', unitlist);
    fprintf('Done\n')
%     % some more meta data
%     % some more meta data
%     ops = io.loadOps(S.rawFilePath);
%     if numel(ops)>1
%         ops = ops(1);
%     end
    
%     ops = io.convertOpsToNewDirectory(ops, S.rawFilePath);
%     
%     chMap = load(ops.chanMap);
%     Exp.chmap = chMap;
%     
%     S.numChan = sum([ops.Nchan]);
%     S.numSU = sum(Exp.osp.isiV<.2);
%     S.numU = numel(Exp.osp.isiV);
    
    % --- Save Exp Struct to data folder for future use (PROC folder)
    % should consider later how we want to store things (per protocol,
    % but off hand it seems feasible given the generic structure of
    % of the D struct maybe we could concatenate all sessions in one
    Exp.ProcDataFolder = dataPath;
    Exp.DataFolder = getpref('FREEVIEWING', 'SERVER_DATA_DIR');
    Exp.FileTag = S.processedFileName;
    
    save(fname,'-v7.3', '-struct', 'Exp');
    fprintf('Exp struct saved to %s\n',fname);
    fprintf('Mat File %s\n',S.processedFileName);
    
    % --- save meta data to dataset table
    fprintf('Updating Dataset Table\n')
    
    if strcmp(thisSession.ImportFun, 'importFreeViewingHuklab')
        disp('Updating protocol list')
        protocols =  {'Grating', 'Gabor', 'Dots', 'BackImage', ...
            'FixRsvpStim', ...
            'FixCalib', ...
            'ForageStaticLines', ...
            'DriftingGrating'};
        
        prot_list = [];
        for iprot = 1:numel(protocols)
            vt = io.getValidTrials(Exp, protocols{iprot});
            if ~isempty(vt)
                if isempty(prot_list)
                    prot_list = protocols{iprot};
                else
                    prot_list = [prot_list ',' protocols{iprot}];
                end
            end
        end
        
        thisSession.StimulusProtocols = {prot_list};
    end
%     thisSession.numChannels = S.numChan;
%     thisSession.numUnits = S.numU;
%     thisSession.numSingleUnit = S.numSU;
    
    data(sessionId,:) = thisSession;
    writetable(data, meta_file);
    fprintf('Done\n')
end

varargout{1} = Exp;

if nargout > 1
%     if ~isfield(S, 'cids')
%         S.cids = Exp.osp.cids;
%     end

    varargout{2} = S;
    
    if nargout > 2 
        error("dataFactoryTreadmill: more than two outputs requested. I don't know what to do with that")
    end
        
end
    