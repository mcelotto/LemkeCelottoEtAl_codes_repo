function [FIT_biased, varargout] = FIT_choice(C, hX, hY, Y, varargin)

%%% *function [FIT_biased, varargout] = FIT_choice(C, hX, hY, Y, varargin)*
%%%
%%% ### Description
%%% This function computes the feature-specific information transmission FIT for choice (or behavioral output). Either naive or bias-corrected, given input data.
%%%
%%% ### Inputs:
%%% - *C*: must be an array of *nDimsS X nTrials * elements representing the discrete value of the choice on each trial.
%%% - *hX*: must be an array of *nDimsX X nTrials* response matrix describing the past of the emitter variables on each of the *nDims* dimensions for each trial.
%%% - *hY*: must be an array of *nDimsY X nTrials* response matrix describing the past of the receiver variable on each of the *nDims* dimensions for each trial.
%%% - *Y*: must be an array of *nDimsY X nTrials* response matrix describing the response of each of the *nDims* dimensions for each trial.
%%% - *opts*: options used to calculate II (see further notes).
%%%
%%% ### Outputs:
%%% - *FIT_biased*: naive, biased estimation of feature-specific information transmission $II = min[SI(S:{R,C}), SI(C:{R,S})]$ for the input data.
%%% - *FIT_unbiased*: unbiased estimate of feature-specific information transmission.  This is returned as second argument only if `opts.bias` is not `'naive'`.
%%% - *FIT p-value*: pvalue of the measure tested against the null hypothesis of Y being conditionally independent from X given S
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field                                | description                                                                                                                                                                                                                                             | allowed values                                                                                                                                                                                                                                     | default   |
%%% |--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
%%% | opts.bias                            | specifies the bias correction method                                                                                                                                                                                                                    | `'naive'` (no bias correction)<br>`'le'` (linear extrapolation)<br>`'qe'` (quadratic extrapolation)                                                                                                                                                | `'naive'` |
%%% | opts.max_draws_per_split_number      | specifies the maximum number of draws on which to calculate the unbiased estimate of II. E.g. `opts.max_draws_per_split_number = 10`, the II is calculated as average on 20 trials for 2 splits and 10 trials for 1 split. (ignored if opts.bias = 'naive') | int > 0                                                                                                                                                                                                                                        | 50        |
%%% | opts.th_bias_corr_convergence_slope  | threshold slope of the standard error of the mean vs number of draws curve (linear fit over the last 5 draws) to judge convergence of the estimated II                                                                                                  | float > 0                                                                                                                                                                                                                                          | 0.5E-4    |
%%% | opts.bin_method_X                    | specifies the binning method for the X signals                                                                                                                                                                                                          | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_Y                    | specifies the binning method for the Y signals                                                                                                                                                                                                          | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_C                    | specifies the binning method for the stimulus signals                                                                                                                                                                                                   | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.n_binsX                         | number of bins to be used to reduce the dimensionality of the response X                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                            | 3         |
%%% | opts.n_binsY                         | number of bins to be used to reduce the dimensionality of the response Y                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                            | 3         |
%%% | opts.n_binsC                         | number of bins to be used to reduce the dimensionality of the stimulus S                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                            | 2         |
%%% | opts.nullhyp                         | If == 1 test the FIT value against the null hypothesis of Y being conditionally independent from X once S is given, via a permutation test. If == 0 don't test significance
%%% | opts.nh_perm                         | number of permutations used to build the null hypothesis distribution, by default 200
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
assert( length(C(1,:)) == length(hX(1,:)) || ...
     length(C(1,:)) == length(hY(1,:)) || ...
     length(C(1,:)) == length(hY(1,:)), ...
     'Number of trials are differing among the inputs')
    
assert(length(hX(1,:)) == length(C(1,:)),...
    "Number of columns of hX, hX(1,:), should be of size nTrials");
assert(length(hY(1,:)) == length(C(1,:)),...
    "Number of columns of hY, hY(1,:), should be of size nTrials");
assert(length(Y(1,:)) == length(C(1,:)),...
    "Number of columns of Y, Y(1,:), should be of size nTrials");
assert(length(hY(:,1)) == length(Y(:,1)), ...
    'Dimensions must match in hY and Y');

% check for nans
if any(isnan(C))
    error("S contains NaNs. Aborting.")
end
if any(isnan(hX))
    error("hX contains NaNs. Aborting.")
end
if any(isnan(hY))
    error("hY contains NaNs. Aborting.")
end
if any(isnan(Y))
    error("Y contains NaNs. Aborting.")
end

% set default options
if nargin == 4
    opts.bias = 'naive';
    opts.bin_method_X = 'eqpop';
    opts.bin_method_Y = 'eqpop';
    opts.bin_method_C = 'none';
    opts.n_binsX = 3;
    opts.n_binsY = 3;
    opts.n_binsC = 2;
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
    % check th_bias_corr_convergence_slope and set default if not provided
    if ~isfield(opts,'th_bias_corr_convergence_slope')
        opts.th_bias_corr_convergence_slope = 0.5E-8;
    else
        assert(opts.th_bias_corr_convergence_slope > 0,...
            'opts.th_bias_corr_convergence_slope should be > 0');
    end    
    % check binning method and set default if not provided
    if ~isfield(opts,'bin_method_X')
        opts.bin_method_X = 'eqpop';
    else
        assert(strcmp(opts.bin_method_X,'none') ||...
            strcmp(opts.bin_method_X,'eqpop') ||...
            strcmp(opts.bin_method_X, 'eqspace') ||...
            strcmp(opts.bin_method_X, 'ceqspace') ||...
            strcmp(opts.bin_method_X, 'gseqspace'),...
            ['opts.bin_method_X argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_X);
    end
    if ~isfield(opts,'bin_method_Y')
        opts.bin_method_Y = 'eqpop';
    else
        assert(strcmp(opts.bin_method_Y,'none') ||...
            strcmp(opts.bin_method_Y,'eqpop') ||...
            strcmp(opts.bin_method_Y, 'eqspace') ||...
            strcmp(opts.bin_method_Y, 'ceqspace') ||...
            strcmp(opts.bin_method_Y, 'gseqspace'),...
            ['opts.bin_method_Y argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_Y);
    end
    if ~isfield(opts,'bin_method_C')
        opts.bin_method_C = 'none';
    else
        assert(strcmp(opts.bin_method_C,'none') ||...
            strcmp(opts.bin_method_C,'eqpop') ||...
            strcmp(opts.bin_method_C, 'eqspace') ||...
            strcmp(opts.bin_method_C, 'ceqspace') ||...
            strcmp(opts.bin_method_C, 'gseqspace'),...
            ['opts.bin_method_C argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_C);
    end

    % check n_bins and set default if not provided
    if ~isfield(opts,'n_binsX')
        opts.n_binsX = 3;
    else
        assert(opts.n_binsX > 1,...
            'opts.n_binsX should be > 1');
    end 
    if ~isfield(opts,'n_binsY')
        opts.n_binsY = 3;
    else
        assert(opts.n_binsY > 1,...
            'opts.n_binsY should be > 1');
    end     
    if ~isfield(opts.n_binsC,'n_binsC')
        opts.n_binsC = 2;
    else
        assert(opts.n_binsC > 1,...
            'opts.n_binsC should be > 1');
    end   
    
    % check nullhyp and set default if not provided
    if ~isfield(opts,'nullhyp')
        opts.nullhyp = 0;
    else
        assert((opts.nullhyp == 0) || (opts.nullhyp == 1),...
            'opts.nullhyp should 0 or 1');
    end  
    if ~isfield(opts,'nh_perm') 
        opts.nh_perm = 0;
    else
        assert(opts.nh_perm >= 0,...
            'opts.nh_perm should be >= 0');
    end       
end
%%
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

n_trials = length(C(1,:));
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
        opts.bin_method_X, sf, floor(nchoosek(n_trials,n_split_trials)*sf));
end

if min(C) <= 0
    C = C - min(C) + 1; % avoid negative value for the probabilityDist function
end
% bin neural responses if required
if ~strcmp(opts.bin_method_C, 'none')
    C_b = binr(C, opts.n_binsC, opts.bin_method_C) + 1;
else
    C_b = C + 1;
end

if ~strcmp(opts.bin_method_X, 'none')
    hX_b = binr(hX, opts.n_binsX, opts.bin_method_X) + 1; % The '+ 1' is added because zeroes in the array would cause an error in the probabilityDist function
else
    hX_b = hX + 1;
end

if ~strcmp(opts.bin_method_Y, 'none')
    hY_b = binr(hY, opts.n_binsY, opts.bin_method_Y) + 1;
    Y_b = binr(Y, opts.n_binsY, opts.bin_method_Y) + 1;
else
    hY_b = hY + 1;
    Y_b = Y + 1;
end

%%
% if multi-dimensional response map it onto a 1D space
% if multi-dimensional response map it onto a 1D space
if length(hX_b(:,1)) > 1
    hX_b = map_Nd_resp_to_1d(hX_b);
end
if length(hY_b(:,1)) > 1
    hY_b = map_Nd_resp_to_1d(hY_b);
    Y_b = map_Nd_resp_to_1d(Y_b);
end
if length(C_b(:,1)) > 1
    C_b = map_Nd_resp_to_1d(C_b);
end

% calculate FIT
FIT_splits = zeros(1,n_splits);

for s = 1:n_splits
    sf = split_factors(s);
    n_split_trials(s) = floor(n_trials*sf);
    if sf == 1
        max_n_draws = 1;
    else
        max_n_draws = round(opts.max_draws_per_split_number/sf);
    end
    FIT_tmp = zeros(1,max_n_draws);
    st_err_mean = zeros(1,max_n_draws);
    
    for d = 1:max_n_draws
        % fetch parallel results as they arrive
        FIT_tmp(d) = calculate_choice_specific_information_transmission(C_b',hX_b,hY_b,Y_b,n_trials,n_split_trials(s));
        % calculate standard error of the mean based on number of samples so far
        st_err_mean(d) = std(FIT_tmp(1:d))/sqrt(length(FIT_tmp(1:d)));
    end
    FIT_splits(s) = mean(FIT_tmp);
end
% biased FIT
FIT_biased = zeros(1, opts.nh_perm+1); 
FIT_biased(1) = FIT_splits(1);
% bias-corrected FIT
if strcmp(opts.bias, 'le')
    p = polyfit(1./n_split_trials, FIT_splits, 1);
    varargout{1} = p(2);
elseif strcmp(opts.bias, 'qe')
    p = polyfit(1./n_split_trials, FIT_splits, 2);
    varargout{1} = p(3);
end

if opts.nullhyp 
    % The null hypothesis distribution and the measure pvalue are computed
    % for the biased quantity. Indeed the surrogate distributions built to
    % compute the null hypothesis are affected, on average, by the same
    % limited-samplig bias of the original distribution.
    
    data = [C_b; hX_b; hY_b; Y_b];
    
    [FIT_null_distribution, pval] = conditional_independence_nullhypothesis(data,FIT_biased(1),opts);
    FIT_biased(2:end) = FIT_null_distribution;
    varargout{2} = pval;
else
    varargout{2} = nan;
end
end
