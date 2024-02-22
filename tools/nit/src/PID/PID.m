function [PID, varargout] = PID(Y, X1, X2, varargin)
%%% *function [PID_biased, varargout] = PID(Y, X1, X2, varargin)*
%%%
%%% ### Description
%%% This function computes the partial information decomposition of X1 and
%%% X2 about Y
%%%
%%% ### Inputs:
%%% - *Y*: must be an array of *nDimensions X nTrials* elements representing the first signal
%%% - *X1*: must be an array of *nDimensions X nTrials* elements representing the second signal
%%% - *X2*: must be an array of *nDimensions X nTrials* elements representing the third signal
%%% - *opts*: options used to calculate PID (see further notes).
%%%
%%% ### Outputs:
%%% - *PID*: structure containing the naive, biased estimation for partial information decomposition
%%% - *PID_unbiased*: structure containing unbiased estimate for partial information decomposition. This is returned as second argument only if `opts.bias` is not `'naive'`.
%%% decomposition
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field                                | description                                                                                                                                                                                                                                             | allowed values                                                                                                                                                                                                                                     | default   |
%%% |--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
%%% | opts.bin_method_X1                    | specifies the binning method for the X1 signal                                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_X2                    | specifies the binning method for the X2 signal                                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_Y                    | specifies the binning method for the Y signal                                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.max_draws_per_split_number      | specifies the maximum number of draws on which to calculate the unbiased estimate of the PID. E.g. `opts.max_draws_per_split_number = 10`, the II is calculated as average on 20 trials for 2 splits and 10 trials for 1 split. (ignored if opts.bias = 'naive') | int > 0                                                                                                                                                                                                                                        | 50        |
%%% | opts.n_bins_X1                          | number of bins to be used to reduce the dimensionality of X1  | int > 1 | 3         |
%%% | opts.n_bins_X2                          | number of bins to be used to reduce the dimensionality of X2  | int > 1 | 3         |
%%% | opts.n_bins_Y                          | number of bins to be used to reduce the dimensionality of Y  | int > 1 | 3         |
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
assert(length(Y(1,:)) == length(X1(1,:)) && length(Y(1,:)) == length(X2(1,:)),...
    "Number of trials are differing in S, R and C");

%Y = Y';
% check for nans
if any(isnan(Y))
    error("Y contains NaNs. Aborting.")
end
if any(isnan(X1))
    error("X1 contains NaNs. Aborting.")
end
if any(isnan(X2))
    error("X2 contains NaNs. Aborting.")
end

% set default options
if nargin == 3
    opts.bin_method_X1 = 'eqpop';
    opts.bin_method_X2 = 'eqpop';
    opts.bin_method_Y = 'eqpop';
    opts.n_bins_X1 = 3;
    opts.n_bins_X2 = 3;
    opts.n_bins_Y = 3;
    opts.max_draws_per_split_number = 200;
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
    if isfield(opts,'bin_method_X1')
        assert(strcmp(opts.bin_method_X1,'none') ||...
            strcmp(opts.bin_method_X1,'eqpop') ||...
            strcmp(opts.bin_method_X1, 'eqspace') ||...
            strcmp(opts.bin_method_X1, 'ceqspace') ||...
            strcmp(opts.bin_method_X1, 'gseqspace'),...
            ['opts.bin_method_X1 argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_X1);
    else
        opts.bin_method_X1 = 'none';
    end
    if isfield(opts,'bin_method_X2')
        assert(strcmp(opts.bin_method_X2,'none') ||...
            strcmp(opts.bin_method_X2,'eqpop') ||...
            strcmp(opts.bin_method_X2, 'eqspace') ||...
            strcmp(opts.bin_method_X2, 'ceqspace') ||...
            strcmp(opts.bin_method_X2, 'gseqspace'),...
            ['opts.bin_method_X1 argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_X2);
    else
        opts.bin_method_X2 = 'none';
    end
    if isfield(opts,'bin_method_Y')
        assert(strcmp(opts.bin_method_Y,'none') ||...
            strcmp(opts.bin_method_Y,'eqpop') ||...
            strcmp(opts.bin_method_Y, 'eqspace') ||...
            strcmp(opts.bin_method_Y, 'ceqspace') ||...
            strcmp(opts.bin_method_Y, 'gseqspace'),...
            ['opts.bin_method argument can be only ''eqpop'', ''eqspace'', ''ceqspace'' or ''gseqspace''. ',...
            'Specified value ''%s'''], opts.bin_method_Y);
    else
        opts.bin_method_Y = 'none';
    end
    % check n_bins and set default if not provided
    if ~isfield(opts,'n_bins_X1')
        opts.n_bins_X1 = 3;
    else
        assert(opts.n_bins_X1 > 1,...
            'opts.n_bins_X1 should be > 1');
    end
    if ~isfield(opts,'n_bins_X2')
        opts.n_bins_X2 = 3;
    else
        assert(opts.n_bins_X2 > 1,...
            'opts.n_bins_X2 should be > 1');
    end
    if ~isfield(opts,'n_bins_Y')
        opts.n_bins_Y = 3;
    else
        assert(opts.n_bins_Y > 1,...
            'opts.n_bins_Y should be > 1');
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

n_trials = length(X1(1,:));
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
        opts.bin_method_X1, sf, floor(nchoosek(n_trials,n_split_trials)*sf));
end

% bin X if required
if ~strcmp(opts.bin_method_X1, 'none')
    X1_b = binr(X1, opts.n_bins_X1, opts.bin_method_X1);
else
    X1_b = X1;
end
if ~strcmp(opts.bin_method_X2, 'none')
    X2_b = binr(X2, opts.n_bins_X2, opts.bin_method_X2);
else
    X2_b = X2;
end
if ~strcmp(opts.bin_method_Y, 'none')
    Y_b = binr(Y, opts.n_bins_Y, opts.bin_method_Y);
else
    Y_b = Y;
end

% if multi-dimensional response map it onto a 1D space
if length(X1_b(:,1)) > 1
    X1_b = map_Nd_array_to_1d(X1_b);
end
if length(X2_b(:,1)) > 1
    X2_b = map_Nd_array_to_1d(X2_b);
end
if length(Y_b(:,1)) > 1
    Y_b = map_Nd_array_to_1d(Y_b);
end

% calculate PID
PIDsplits = cell(1,n_splits);

for s = 1:n_splits
    sf = split_factors(s);
    n_split_trials(s) = floor(n_trials*sf);
    if sf == 1
        max_n_draws = 1;
    else
        max_n_draws = round(opts.max_draws_per_split_number/sf);
    end
    PIDsi_tmp = zeros(1,max_n_draws);
    PIDuiy_tmp = zeros(1,max_n_draws);
    PIDuiz_tmp = zeros(1,max_n_draws);
    PIDci_tmp = zeros(1,max_n_draws);
    st_si_err_mean = zeros(1,max_n_draws);
    st_uiy_err_mean = zeros(1,max_n_draws);
    st_uiz_err_mean = zeros(1,max_n_draws);
    st_ci_err_mean = zeros(1,max_n_draws);
    
    for d = 1:max_n_draws
        % fetch parallel results as they arrive
        PID_tmp = calculate_pid(X2_b,Y_b,X1_b,n_trials,n_split_trials(s));
        PIDsi_tmp(d) = PID_tmp(1);
        PIDuiy_tmp(d) = PID_tmp(2);
        PIDuiz_tmp(d) = PID_tmp(3);
        PIDci_tmp(d) = PID_tmp(4);

        % calculate standard error of the mean based on number of samples so far
        st_si_err_mean(d) = std(PIDsi_tmp(1:d))/sqrt(length(PIDsi_tmp(1:d)));
        st_uiy_err_mean(d) = std(PIDuiy_tmp(1:d))/sqrt(length(PIDuiy_tmp(1:d)));
        st_uiz_err_mean(d) = std(PIDuiz_tmp(1:d))/sqrt(length(PIDuiz_tmp(1:d)));
        st_ci_err_mean(d) = std(PIDci_tmp(1:d))/sqrt(length(PIDci_tmp(1:d)));
    end
    PIDsi_splits(s) = mean(PIDsi_tmp);
    PIDuiy_splits(s) = mean(PIDuiy_tmp);
    PIDuiz_splits(s) = mean(PIDuiz_tmp);
    PIDci_splits(s) = mean(PIDci_tmp);
end

% biased II
PID.shared = mean(PIDsi_tmp);
PID.uniqueX1 = mean(PIDuiy_tmp);
PID.uniqueX2 = mean(PIDuiz_tmp);
PID.complementary = mean(PIDci_tmp);
% bias-corrected II
if strcmp(opts.bias, 'le')
    p = polyfit(1./n_split_trials, PIDsi_splits, 1);
    PID_unbiased.shared = p(2);
    p = polyfit(1./n_split_trials, PIDuiy_splits, 1);
    PID_unbiased.uniqueX1 = p(2);
    p = polyfit(1./n_split_trials, PIDuiz_splits, 1);
    PID_unbiased.uniqueX2 = p(2);
    p = polyfit(1./n_split_trials, PIDci_splits, 1);
    PID_unbiased.complementary = p(2);
    varargout{1} = PID_unbiased;
elseif strcmp(opts.bias, 'qe')
    p = polyfit(1./n_split_trials, PIDsi_splits, 2);
    PID_unbiased.shared = p(3);
    p = polyfit(1./n_split_trials, PIDuiy_splits, 2);
    PID_unbiased.uniqueX1 = p(3);
    p = polyfit(1./n_split_trials, PIDuiz_splits, 2);
    PID_unbiased.uniqueX2 = p(3);
    p = polyfit(1./n_split_trials, PIDci_splits, 2);
    PID_unbiased.complementary = p(3);
    varargout{1} = PID_unbiased;
end

end
