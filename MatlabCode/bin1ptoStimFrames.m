function [Y, binfun] = bin1ptoStimFrames(dF, eventTimes, imgfr_ts, binSize, keep_sparse)
% Creating analogous function to binspiketimesfast but for widefield
% imaging, -> convert dF at imaging rate to dF at stimulus rate, Y will be
% dF in the format nStimulusframes x nPixels, as in Rate(frame) x unit

% BIN SPIKE TIMES THE FASTEST WAY POSSIBLE IN MATLAB
% 
% Inputs:
%   sp [struct]: Kilosort output struct
%   has fields:
%       st [T x 1]: spike times
%       clu [T x 1]: unit id
% Outpus:
%   Y [nBins x nNeurons]
%
% Example Call:
%   Y = binNeuronSpikeTimesFast(Exp.osp, eventTimes, binsize)

% NT = numel(eventTimes);
% NC = max(sp.cids);

% bs = min(diff(eventTimes)) - eps;
% if bs < binSize
%     warning('maximum bin size is the minimum distance between events')
%     binSize = bs;
% end

% eventTimes = [eventTimes(:)'; eventTimes(:)' + binSize];
% eventTimes = eventTimes(:);
% % eventTimes = [eventTimes(:); eventTimes(end) + binSize];
% Y = zeros(NT, NC);
% for cc = sp.cids(:)'
%     cnt = histcounts(sp.st(sp.clu==cc), eventTimes);
% %     cnt(diff(eventTimes) > .1) = 0;
%     Y(:,cc) = cnt(1:2:end);
% end


% return

if nargin < 5
    keep_sparse = false;
end

%Vectorise video to pixels x time if needed
if length(size(dF))>2
    n_img_frames=size(dF,3);
    dF=reshape(dF,[],n_img_frames);
end
n_img_frames=size(dF,2);
npix=size(dF,1);

% conversion from time to bins
binfun = @(x) (x==0) + ceil(x / binSize);

% bin imaging frame times, into stim framerate, placeholder
bst = binfun(imgfr_ts); % Puts img frames into stimframes

% bin frame times, eventTimes are trialstarts in seconds
bft = binfun(eventTimes);

% % create binned spike times
% Y = sparse(bst, double(sp.clu), ones(numel(bst), 1), max(max(bst),max(bft)),double(max(sp.clu)));
% 
% % index in with binned frame times
% Y = Y(bft,:);

%Surely not the fastest way to do this
Y=zeros(length(eventTimes),npix);
for ii = 1:npix
    Y(:,ii) = interp1(bst,dF(ii,:),bft); %simple lookup table for ft to imaging time
end


if ~keep_sparse
    Y = full(Y);
end




