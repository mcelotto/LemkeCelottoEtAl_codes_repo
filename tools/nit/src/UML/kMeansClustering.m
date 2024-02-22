function [sPredicted, S, accuracy, varargout] = kMeansClustering(R, S, kMeansOpts, varargin)
%%% *function [sPredicted, accuracy, varargout] = kMeansClustering(R, S, kMeansOpts, varargin)*
%%%
%%% ### Description
%%% k-means clustering algorithm for multidimensional neural acitvity. Allows for
%%% computation of mutual information between the estimated clusters and
%%% the original stimuli.
%%%
%%% ### Inputs:
%%% - *R*: *nDimensions X nTrials* array including response data.
%%% - *S*: *1 X nTrials* stimulus array.
%%% - *kMeansOpts* (optional): options structure for kMeans (see further notes section for more
%%%           details. 
%%% - *MIopts* (optional): options structure for computing MI (see further notes section for more
%%%           details. Only used if kMeansOpts.MI is TRUE
%%% ### Outputs:
%%% - *outputs*: vector of same length as *outputsList* returning the specified outputs.
%%%
%%% ### Further notes:
%%% #### The options structure
%%% The kMeansOpts structure can include any the following fields:
%%% - kMeansOpts.MI: TRUE for computing MI between predicted clusters and
%%% original stimuli in S. Default is FALSE
%%% - kMeansOpts.distance: Distance metric, in p-dimensional space, used
%%% for minimization, specified as the comma-separated pair consisting of
%%% 'Distance' and 'sqeuclidean', 'cityblock', 'cosine', 'correlation', or
%%% 'hamming'. Default 'sqeuclidean'. Further info in 
%%% https://www.mathworks.com/help/stats/kmeans.html#namevaluepairarguments
%%% - kMeansOpts.parallel: TRUE for using parallel pool (only if parallel
%%% computing toolbox installed). Default FALSE
%%% The opts structure can include any the following fields:
%%% - opts.nbins: number of bins array (can be a scalar or a vector of `size(nb) = size(R,1)`). If a scalar is specified then all dimensions are binned using the same number of bins. Alternatively, each dimension is binned using its specific number of bins. 
%%% - opts.bin_method: binning method (can be a scalar or a cell array of `size(method) = size(R,1)`). If a scalar is specified then all dimensions are binned using the same method. Alternatively, each dimension is binned using its specific method.
%%% Possible options for binning methods are:


assert(nargin<5,...
    'Number of input arguments larger than expected')
assert(nargin>1,...
    'Not enough input arguments')

% flags TRUE or FALSE
if exist('kMeansOpts','var')
    if isfield(kMeansOpts,'MI')
        MIflag = kMeansOpts.MI;
    else
        MIflag = 0;
    end
    if ~isfield(kMeansOpts,'parallel')
        kMeansOpts.parallel = 0;
    else
        if kMeansOpts.parallel == 1
            try
                s = eval("ver('parallel')");
                clearvars s
            catch
                kMeansOpts.parallel = 0;
                warning('You set kMeansOpts.parallel TRUE but Parallel Computing Toolbox is not installed')
            end
        end
    end
    if ~isfield(kMeansOpts,'distance')
        kMeansOpts.distance = 'sqeuclidean';
    end
else
    MIflag = 0;
    kMeansOpts.parallel = 0;
    kMeansOpts.distance = 'sqeuclidean';
end

% adapt data from NIT format to MATLAB k-means
R = R';
nBins = length(unique(S));

sPredicted = kmeans(R,nBins,'distance',kMeansOpts.distance,...
    'Options',statset('UseParallel',kMeansOpts.parallel));
sPredicted = sPredicted';

nClusters = length(unique(S));
uniqueS = unique(S);
newUniqueS = 1:nClusters;
newS = nan(size(S));
for i = 1:nClusters
    idx = S==uniqueS(i);
    newS(idx) = i;
end
S = newS;
uniqueS = unique(S);
clearvars newS idx
    
permsS = perms(uniqueS);
% [~, index] = ismember(uniqueS,permsS,'rows');

accuracy = 0;
newSpredicted = [];
for i = 1:size(permsS,1)
    currentPerm = permsS(i,:);
    sPredictedTemp = nan(size(S));
    for j = 1:nClusters
        idx = sPredicted==j;
        sPredictedTemp(idx) = currentPerm(j);
    end
    accuracyTemp = sum(sPredictedTemp == S)/length(S);
    if accuracyTemp > accuracy
       accuracy=accuracyTemp;
       newSpredicted = sPredictedTemp;
    end
end
sPredicted = newSpredicted;
clearvars accuracyTemp newSpredicted currentPerm sPredictedTemp
end

