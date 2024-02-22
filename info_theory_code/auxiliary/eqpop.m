function binnedX = eqpop(X, nb, varargin)
%%% *function binnedX = eqpop(X, nb, varargin)*
%%%
%%% ### Description
%%% Builds edges for equipopulated binning. Unlike the other binning functions provided, that bin data according to the values of X, equally-populated binning bins according to the count. For this reason, this function returns already the binned version of X and **not** the bin edges.
%%%
%%% ### Inputs:
%%% - *X*: array of values to be discretized.
%%% - *nb*: number of bins.
%%%
%%% ### Outputs:
%%% - *binnedX*: unlike the other binning functions, that bin.
%%%

uniqueX = unique(X(:));

Nu = length(uniqueX);
N = length(X);

if Nu<nb
    error("Too many bins for the selected data.");
end

binSize = (N+1)/nb;

if sum(X == mode(X)) > binSize
    warning('eqpop:tooFrequentResponse',"Using equally populated bins with a response array showing at least one response that is more frequent than nSamples/nBins. The same response will be binned differently in different trials: this may artificially increase the entropy of X. Consider changing the binning strategy.")
end

[~ , sortedIndxs] = sort(X);

edges = round(linspace(1,length(X)+1,nb+1));
edges = edges-1;
edges(1) = 1;

% eqpop returns edges as indices, not as values of the response
edges = repelem(edges,2);
edges(1) = [];
edges(end) = [];
edges(3:2:end-1) = edges(3:2:end-1) + 1;
edges = reshape(edges,2,nb);
binnedX = nan(size(X));
for tb = 1:nb
    st = edges(1,tb);
    en = edges(2,tb);
    binnedX(sortedIndxs(st:en)) = tb;
end
