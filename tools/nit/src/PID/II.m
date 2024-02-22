function [II_biased, varargout] = II(S, R, C, varargin)
%%% *function [II_biased, varargout] = II(S, R, C, varargin)*
%%%
%%% ### Description
%%% This function computes the intersection information II, either naive or bias-corrected, given input data.
%%%
%%% ### Inputs:
%%% - *S*: must be an array of *nDimsS X nTrials* elements representing the discrete value of the stimulus presented in each trial.
%%% - *R*: must be an array of *nDimsR X nTrials* response matrix describing the response of each of the *nDims* dimensions for each trial.
%%% - *C*: must be an array of *nDimsC X nTrials* elements representing the discrete value of the choice made by the subject in each trial.
%%% - *opts*: options used to calculate II (see further notes).
%%%
%%% ### Outputs:
%%% - *II_biased*: naive, biased estimation of intersection information $II = min[SI(S:{R,C}), SI(C:{R,S})]$ for the input data.
%%% - *II_unbiased*: unbiased estimate of intersection information.  This is returned as second argument only if `opts.bias` is not `'naive'`.
%%% - *II_null*: null distribution if opts.nullhyp = 1
%%% - *pval*: p value for the statistical test if opts.nullhyp = 1
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field                                | description                                                                                                                                                                                                                                                 | allowed values                                                                                                                                                                                                                                     | default   |
%%% |--------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
%%% | opts.bias                            | specifies the bias correction method                                                                                                                                                                                                                        | `'naive'` (no bias correction)<br>`'le'` (linear extrapolation)<br>`'qe'` (quadratic extrapolation)                                                                                                                                                | `'naive'` |
%%% | opts.max_draws_per_split_number      | specifies the maximum number of draws on which to calculate the unbiased estimate of II. E.g. `opts.max_draws_per_split_number = 10`, the II is calculated as average on 20 trials for 2 splits and 10 trials for 1 split. (ignored if opts.bias = 'naive') | int > 0                                                                                                                                                                                                                                        | 50        |
%%% | opts.bin_method_stimulus             | specifies the binning method for the stimulus                                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'none'` |
%%% | opts.bin_method_activity             | specifies the binning method for the activity                                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_choice               | specifies the binning method for the behavior                                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'none'` |
%%% | opts.n_bins_stimulus                 | number of bins to be used to reduce the dimensionality of the stimulus  | int > 1 | 3         |
%%% | opts.n_bins_activity                 | number of bins to be used to reduce the dimensionality of the activity  | int > 1 | 3         |
%%% | opts.n_bins_choice                   | number of bins to be used to reduce the dimensionality of the choice    | int > 1 | 3         |
%%% | opts.nullhyp                         | If == 1 test the II value against the null hypothesis of Y being conditionally independent from X once S is given, via a permutation test. If == 0 don't test significance
%%% | opts.nh_perm                         | number of permutations used to build the null hypothesis distribution, by default 200
%%%
%
%  This source code is part of:
%  NIT - Neuroscience Information Toolbox
%  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.

% check inputs dimensions
assert(length(R(1,:)) == length(S(1,:)) && length(R(1,:)) == length(C(1,:)),...
    "Number of trials are differing in S, R and C");

% check for nans
if any(isnan(S))
    error("S contains NaNs. Aborting.")
end
if any(isnan(R))
    error("R contains NaNs. Aborting.")
end
if any(isnan(C))
    error("C contains NaNs. Aborting.")
end

% set default options
if nargin == 3
    opts.bias = 'naive';
    opts.bin_method_activity = 'eqpop';
    opts.bin_method_stimulus = 'none';
    opts.bin_method_choice = 'none';
    opts.n_bins_stimulus = 3;
    opts.n_bins_activity = 3;
    opts.n_bins_choice = 3;
    opts.max_draws_per_split_number = 200;
    opts.nullhyp = 0;
    opts.nh_perm = 0;
else
    opts = varargin{1};
    % check bias reduction method and set default if not provided
    if ~isfield(opts,'bias')
        opts.bias = 'naive';
    else
        assert(strcmp(opts.bias,'naive') ||...
            strcmp(opts.bias, 'le') ||...
            strcmp(opts.bias, 'qe'),...
            ['opts.bias argument can be only ''naive'', ''le'' or ''qe''. ',...
            'Specified value ''%s'''], opts.bias);
    end
    % check draws_per_split_number and set default if not provided
    if ~isfield(opts,'max_draws_per_split_number')
        opts.max_draws_per_split_number = 50;
    else
        assert(opts.max_draws_per_split_number > 0,...
            'opts.max_draws_per_split_number should be > 0');
    end   
    % check binning method and set default if not provided
    if ~isfield(opts,'bin_method_activity')
        opts.bin_method_activity = 'none';
    else
        assert(strcmp(opts.bin_method_activity,'none') ||...
            strcmp(opts.bin_method_activity,'eqpop') ||...
            strcmp(opts.bin_method_activity, 'eqspace') ||...
            strcmp(opts.bin_method_activity, 'ceqspace') ||...
            strcmp(opts.bin_method_activity, 'gseqspace'),...
            ['opts.bin_method_activity argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_activity);
    end
    if ~isfield(opts,'bin_method_stimulus')
        opts.bin_method_stimulus = 'none';
    else
        assert(strcmp(opts.bin_method_stimulus,'none') ||...
            strcmp(opts.bin_method_stimulus,'eqpop') ||...
            strcmp(opts.bin_method_stimulus, 'eqspace') ||...
            strcmp(opts.bin_method_stimulus, 'ceqspace') ||...
            strcmp(opts.bin_method_stimulus, 'gseqspace'),...
            ['opts.bin_method_stimulus argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_stimulus);
    end
    if ~isfield(opts,'bin_method_choice')
        opts.bin_method_choice = 'none';
    else
        assert(strcmp(opts.bin_method_choice,'none') ||...
            strcmp(opts.bin_method_choice,'eqpop') ||...
            strcmp(opts.bin_method_choice, 'eqspace') ||...
            strcmp(opts.bin_method_choice, 'ceqspace') ||...
            strcmp(opts.bin_method_choice, 'gseqspace'),...
            ['opts.bin_method_choice argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_choice);
    end
    % check n_bins and set default if not provided
    if ~isfield(opts,'n_bins_stimulus')
        opts.n_bins_stimulus = 3;
    else
        assert(opts.n_bins_stimulus > 1,...
            'opts.n_bins_stimulus should be > 1');
    end
    if ~isfield(opts,'n_bins_activity')
        opts.n_bins_activity = 3;
    else
        assert(opts.n_bins_activity > 1,...
            'opts.n_bins_activity should be > 1');
    end
    if ~isfield(opts,'n_bins_choice')
        opts.n_bins_choice = 3;
    else
        assert(opts.n_bins_choice > 1,...
            'opts.n_bins_choice should be > 1');
    end
    if ~isfield(opts,'nullhyp')
        opts.nullhyp = 0;
    else
        assert((opts.nullhyp == 0) || (opts.nullhyp == 1),...
            'opts.nullhyp should 0 or 1');
    end  
    if ~isfield(opts,'nh_perm') 
        opts.nh_perm = 200;
    else
        assert(opts.nh_perm >= 0,...
            'opts.nh_perm should be >= 0');
    end  
end

% set up split factors and other method-depending vars and
% check that number of draws per number of splits is compatible
% with the size of the data (though unlikely to be otherwise)
if strcmp(opts.bias, 'le')
    split_factors = [1 0.5];
elseif strcmp(opts.bias, 'qe')
    split_factors = [1 0.5 0.25];
else
    split_factors = [1];
end

n_trials = length(S(1,:));
n_splits = numel(split_factors);

for i=2:numel(split_factors)
    sf = split_factors(i);
    n_split_trials = floor(n_trials*sf);
    max_n_draws = round(opts.max_draws_per_split_number/sf);
	warning('off','MATLAB:nchoosek:LargeCoefficient')
    assert(max_n_draws <= nchoosek(n_trials,n_split_trials), ['Number of ',...
        'draws per split number ''opts.max_draws_per_split_number'' specified ',...
        'is higher than allowed given sample size. Max allowed for '...
        'method %s and split factor %g is: %i'],...
        opts.bin_method_activity, sf, floor(nchoosek(n_trials,n_split_trials)*sf));
end

% bin R if required
if ~strcmp(opts.bin_method_activity, 'none')
    R_b = binr(R, opts.n_bins_activity, opts.bin_method_activity);
else
    R_b = R;
end
% bin S if required
if ~strcmp(opts.bin_method_stimulus, 'none')
    S_b = binr(S, opts.n_bins_stimulus, opts.bin_method_stimulus);
else
    S_b = S;
end
% bin C if required
if ~strcmp(opts.bin_method_choice, 'none')
    C_b = binr(C, opts.n_bins_choice, opts.bin_method_choice);
else
    C_b = C;
end

% if multi-dimensional response map it onto a 1D space
if length(R_b(:,1)) > 1
    R_b = map_Nd_array_to_1d(R_b);
end
if length(S_b(:,1)) > 1
    S_b = map_Nd_array_to_1d(S_b);
end
if length(C_b(:,1)) > 1
    C_b = map_Nd_array_to_1d(C_b);
end

% calculate II
II_splits = zeros(1,n_splits);

for s = 1:n_splits
    sf = split_factors(s);
    n_split_trials(s) = floor(n_trials*sf);
    if sf == 1
        max_n_draws = 1;
    else
        max_n_draws = round(opts.max_draws_per_split_number/sf);
    end
    II_tmp = zeros(1,max_n_draws);
    st_err_mean = zeros(1,max_n_draws);
    
    for d = 1:max_n_draws
        % fetch parallel results as they arrive
        II_tmp(d) = calculate_intersection_information(S_b',C_b',R_b,n_trials,n_split_trials(s));
        % calculate standard error of the mean based on number of samples so far
        st_err_mean(d) = std(II_tmp(1:d))/sqrt(length(II_tmp(1:d)));
    end
    II_splits(s) = mean(II_tmp);
end
% biased II
II_biased = II_splits(1);
% bias-corrected II
if strcmp(opts.bias, 'le')
    p = polyfit(1./n_split_trials, II_splits, 1);
    varargout{1} = p(2);
elseif strcmp(opts.bias, 'qe')
    p = polyfit(1./n_split_trials, II_splits, 2);
    varargout{1} = p(3);
end

if opts.nullhyp 
    % The null hypothesis distribution and the measure pvalue are computed
    % for the biased quantity. Indeed the surrogate distributions built to
    % compute the null hypothesis are affected, on average, by the same
    % limited-samplig bias of the original distribution.
    data = [S_b; R_b; C_b];
    [II_null_distribution, pval] = conditional_independence_nullhypothesis(data,II_biased,opts);
    varargout{3} = pval;
    varargout{2} = II_null_distribution;
end

end

