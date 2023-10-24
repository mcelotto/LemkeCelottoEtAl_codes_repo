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
%%% - *FIT p-value*: pvalue of the measure tested against the null hypothesis of Y being conditionally independent from X given S
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field                                | description                                                                                                                                                                                                                                             | allowed values                                                                                                                                                                                                                                     | default   |
%%% |--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
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
    opts.bin_method_X = 'eqpop';
    opts.bin_method_Y = 'eqpop';
    opts.bin_method_C = 'none';
    opts.n_binsX = 3;
    opts.n_binsY = 3;
    opts.n_binsC = 2;
    opts.nullhyp = 0;
    opts.nh_perm = 0;
else
    opts = varargin{1};
 
    % check binning method and set default if not provided
    if ~isfield(opts,'bin_method_X')
        opts.bin_method_X = 'eqpop';
    else
        assert(strcmp(opts.bin_method_X,'none') ||...
            strcmp(opts.bin_method_X,'eqpop'),...
            ['opts.bin_method_X argument can be only ''eqpop'' or ''none''',...
            'Specified value ''%s'''], opts.bin_method_X);
    end
    if ~isfield(opts,'bin_method_Y')
        opts.bin_method_Y = 'eqpop';
    else
        assert(strcmp(opts.bin_method_Y,'none') ||...
            strcmp(opts.bin_method_Y,'eqpop'),...
            ['opts.bin_method_Y argument can be only ''eqpop''or ''none''',...
            'Specified value ''%s'''], opts.bin_method_Y);
    end
    if ~isfield(opts,'bin_method_C')
        opts.bin_method_C = 'none';
    else
        assert(strcmp(opts.bin_method_C,'none') ||...
            strcmp(opts.bin_method_C,'eqpop'),...
            ['opts.bin_method_C argument can be only ''eqpop'' or ''none''',...
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
    if ~isfield(opts,'n_binsC')
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

n_trials = length(C(1,:));

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
% biased FIT
FIT_biased = zeros(1, opts.nh_perm+1); 
FIT_biased(1) = calculate_choice_specific_information_transmission(C_b',hX_b,hY_b,Y_b,n_trials,length(C(1,:)));

if opts.nullhyp 
    % The null hypothesis distribution and the measure pvalue are computed
    % for the biased quantity. Indeed the surrogate distributions built to
    % compute the null hypothesis are affected, on average, by the same
    % limited-samplig bias of the original distribution.
    
    data = [C_b; hX_b; hY_b; Y_b];
    
    [FIT_null_distribution, pval] = conditional_independence_nullhypothesis(data,FIT_biased(1),opts);
    FIT_biased(2:end) = FIT_null_distribution;
    varargout{1} = pval;
end
end
