function [prdir, varargout] = prcomp(data, standardize)
%%% *function [prdir, prcomp, prcompvariance, tsq, exp_var] = prcomp(data)*
%%%
%%% ### Description
%%% PCA runs a Principal Components Analysis on the data provided.
%%%
%%% ### Inputs:
%%% - *data*: n_samples x n_features data matrix.
%%% - *standardize*: logical value, true if data should be standardized
%%% (i.e. scaled to have zero-mean and one standard deviation), false otherwise.
%%%
%%% ### Outputs:
%%% - *prdir*: directions of principal components.
%%% - *prcomp*: principal components of the input data.
%%% - *prcompvariance*: principal component variances (eigenvalues of
%%% covariance matrix of data).
%%% - *tsq*: Hotelling’s T-Squared Statistic (sum of squares of the
%%% standardized scores for each observation).
%%% - *exp_var*: explained variance for each of the PC's.

if(nargout > 5)
    error("Too many output arguments.")
end

n_samples = length(data(:,1));
n_features = length(data(1,:));

% standardize data to zero-mean and one std if required
if standardize
    data = data-mean(data,1);
    data = data./std(data);
    data(isnan(data)) = 0;
end

% calculate eigenvalues of covariance matrix only if n_samples >= n_features
if n_samples >= n_features
    cov = (data'*data)./(n_samples-1);
    [prdir, evalues] = eig(cov,'nobalance');
    [varargout{2},i] = sort(diag(evalues),'descend');
    prdir = prdir(:,i);
% else work on Gram matrix for computational performance
else
    % max number of principal components is n_samples - 1 if n_samples < n_features
    n_pc = n_samples-1;
    gram = (data*data')./(n_samples-1);
    [evecs, evalues] = eig(gram,'nobalance');
    [evalues,i] = sort(diag(evalues),'descend');
    evecs = evecs(:,i);
    % principal axes are related to eigenvectors of covariance matrix by:
    prdir = data'*evecs./((n_samples-1)*evalues').^0.5;
    prdir = prdir(:,1:n_pc);
    if nargout >= 3
        varargout{2} = evalues(1:n_pc);
    end
end

% principal components
if nargout >= 3
    varargout{1} = data*prdir;
end
% Hotelling’s T-Squared Statistic
if nargout >= 4
    varargout{3} = varargout{1}-mean(varargout{1},1);
    varargout{3} = varargout{3}./std(varargout{3});
    varargout{3} = sum(varargout{3}.^2,2);
end
% Explained variance
if nargout >= 5
    varargout{4} = (var(varargout{1})./sum(var(varargout{1}))*100)';
end
end