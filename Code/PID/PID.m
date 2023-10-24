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
%%% - *PID_btsp*: cell structure containing naive, biased estimate for partial information decomposition, calculated according to the options `opts.btsp_type` and `opts.btsp_variables`. This is returned as third argument only if bootstrapping is requested. The argument is returned as a cell array, where each component of the cell array returns the results of shuffling according to each component of `opts.btsp_variables` (it is possible to perform multiple bootstrapping operations at once to test for significance of specific atoms).
%%% decomposition
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field                                | description                                                                                                                                                                                                                                                                      | allowed values                                                                                                                                                                                                                                                              | default   |
%%% |--------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
%%% | opts.bin_method_X1                   | specifies the binning method for the X1 signal                                                                                                                                                                                                                                   | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_X2                   | specifies the binning method for the X2 signal                                                                                                                                                                                                                                   | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_method_Y                    | specifies the binning method for the Y signal                                                                                                                                                                                                                                    | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.n_bins_X1                       | number of bins to be used to reduce the dimensionality of X1                                                                                                                                                                                                                     | int > 1                                                                                                                                                                                                                                                                     | 3         |
%%% | opts.n_bins_X2                       | number of bins to be used to reduce the dimensionality of X2                                                                                                                                                                                                                     | int > 1                                                                                                                                                                                                                                                                     | 3         |
%%% | opts.n_bins_Y                        | number of bins to be used to reduce the dimensionality of Y                                                                                                                                                                                                                      | int > 1                                                                                                                                                                                                                                                                     | 3         |
%%% | opts.btsp                            | number of bootstrap operations to perform for significance testing (those will be performed independently on each of the variables listed in `opts.btsp_variables`)                                                                                                              | int >= 0                                                                                                                                                                                                                                                                    | 0         |
%%% | opts.btsp_variables                  | list of variables to be bootstrapped, specified as a cell array of strings. If multiple bootstrapping operations are requested, `opts.btsp_variables` can contain multiple variables. In this case a `opts.btsp_type` option should be specified for each `opts.btsp_variables`  | cell array containing one (or more) strings (`"X1"`, `"X2"` or `"Y"`) corresponding to all variables to be bootstrapped                                                                                                                                                     | N/A       |
%%% | opts.btsp_type                       | type of bootstrapping to be applied for each variable (`'all'` shuffles all values of the corresponding variable in `opts.btsp_variables` across trials, while `'$VAR$conditioned'` shuffles trials by conditioning on values of the variable specified in the substring `$VAR$`)| cell array of same size of opts.btsp_variables, each element contains the type of shuffling to be performed for the variable (either "all"` or `"X1conditioned"`, `"X2conditioned"`, or `"Yconditioned"`)                                                                   | N/A       |
%%%
%%% The bootstrapping options allow the user to tailor the bootstrap
%%% operation to her/his specific needs. One should bear in mind that
%%% unconditioned bootstrapping (`opts.btsp_type = {'all'}`) of one
%%% variable (e.g. `Y`), destroys all MI between `(X1;Y)` as well as the MI
%%% between `(X2;Y)`. The use of a conditioned shuffling (`opts.btsp_type =
%%% {'conditioned'}`), instead, preserves MI values but destroys the
%%% correlations between the other two variables (`X1` and `X2` in the
%%% given example of shuffling `Y`). This can be useful when some PID atoms
%%% can be due to correlations between two of the three variables.
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
    "Number of trials are differing in X1, X2 and Y");

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
    opts.btsp = 0;
else
    opts = varargin{1};
    % check binning method and set default if not provided
    if isfield(opts,'bin_method_X1')
        assert(strcmp(opts.bin_method_X1,'none') ||...
            strcmp(opts.bin_method_X1,'eqpop'),...
            ['opts.bin_method_X1 argument can be only ''eqpop'' or ''none''. ',...
            'Specified value ''%s'''], opts.bin_method_X1);
    else
        opts.bin_method_X1 = 'none';
    end
    if isfield(opts,'bin_method_X2')
        assert(strcmp(opts.bin_method_X2,'none') ||...
            strcmp(opts.bin_method_X2,'eqpop'),...
            ['opts.bin_method_X1 argument can be only ''eqpop'' or ''none''. ',...
            'Specified value ''%s'''], opts.bin_method_X2);
    else
        opts.bin_method_X2 = 'none';
    end
    if isfield(opts,'bin_method_Y')
        assert(strcmp(opts.bin_method_Y,'none') ||...
            strcmp(opts.bin_method_Y,'eqpop'),...
            ['opts.bin_method argument can be only ''eqpop'' or ''none''. ',...
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
    % check bootstrapping options
    if ~isfield(opts,'btsp')
        opts.btsp = 0;
        opts.btsp_variables = {};
        opts.btsp_type = {};
    else
        assert(opts.btsp > 0 && mod(opts.btsp,1) == 0,...
            'opts.btsp should be an integer >= 0')
    end
    if opts.btsp > 0
        assert(isfield(opts,'btsp_variables'),...
            'opts.btsp option requires opts.btsp_variables');
        assert(isfield(opts,'btsp_type'),...
            'opts.btsp option requires opts.btsp_type');
        assert(iscell(opts.btsp_variables),...
            'opts.btsp_variables option should be a cell array');
        assert(iscell(opts.btsp_type),...
            'opts.btsp_type option should be a cell array');
        assert(length(opts.btsp_variables) == length(opts.btsp_type),...
            'cell arrays opts.btsp_variables and opts.btsp_type should have the same size')
        for i = 1:length(opts.btsp_variables)
            assert(strcmp(opts.btsp_variables{i},'Y') ||...
                strcmp(opts.btsp_variables{i},'X1') ||...
                strcmp(opts.btsp_variables{i}, 'X2'),...
                ['Elements in opts.btsp_variables{:} should be one of ''X1'', ''X2'', or ''Y''. ',...
                'Specified value ''%s'''], opts.btsp_variables{i});
            assert(strcmp(opts.btsp_type{i},'all') ||...
                strcmp(opts.btsp_type{i},'X1conditioned') ||...
                strcmp(opts.btsp_type{i},'X2conditioned') ||...
                strcmp(opts.btsp_type{i},'Yconditioned'),...
                ['Elements in opts.btsp_type{:} should be one of ''all'', ''X1conditioned'', ''X2conditioned'', or ''Yconditioned''. ',...
                'Specified value ''%s'''], opts.btsp_type{i});
            splitstring = split(opts.btsp_type{i},'conditioned');
            conditioning_vars{i} = splitstring{1};
            assert(~strcmp(opts.btsp_variables{i},conditioning_vars{i}),...
                ['Shuffle on variable ''%s'' cannot be performed by conditioning on values of the same variable. ', ...
                'The operation has no meaning.'], opts.btsp_variables{i});
        end
    end
end

n_trials = length(X1(1,:));

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

% calculate PID on unshuffled data
PID = calculate_PID(Y_b,X1_b,X2_b,n_trials);

% Calculate PID on shuffled data if required
for s = 1:length(opts.btsp_variables)
    % initialize arrays for btsp results
    btsp_atoms{s}.uniqueX1                  = NaN(1,opts.btsp);
    btsp_atoms{s}.uniqueX2                  = NaN(1,opts.btsp);
    btsp_atoms{s}.shared                    = NaN(1,opts.btsp);
    btsp_atoms{s}.complementary             = NaN(1,opts.btsp);
    btsp_atoms_unbiased{s}.uniqueX1         = NaN(1,opts.btsp);
    btsp_atoms_unbiased{s}.uniqueX2         = NaN(1,opts.btsp);
    btsp_atoms_unbiased{s}.shared           = NaN(1,opts.btsp);
    btsp_atoms_unbiased{s}.complementary    = NaN(1,opts.btsp);
    
    switch opts.btsp_variables{s}
        case 'Y'
            shuffled_var = Y_b;
        case 'X1'
            shuffled_var = X1_b;
        case 'X2'
            shuffled_var = X2_b;
    end
    switch conditioning_vars{s}
        case 'all'
            % un-conditioned case
            cond_var = [];
        case 'X1'
            % X1 conditioned
            cond_var = X1_b;
        case 'X2'
            % X2 conditioned
            cond_var = X2_b;
        case 'Y'
            % Y conditioned
            cond_var = Y_b;
    end
    
    for b = 1:opts.btsp
        if isempty(cond_var)
            shuffled_var = shuffled_var(randperm(length(shuffled_var)));
        else
            unique_vals = unique(cond_var);
            for i = 1:length(unique_vals)
                val_inds = find(cond_var == unique_vals(i));
                shuffled_val_inds = val_inds(randperm(length(val_inds)));
                shuffled_var(val_inds) = shuffled_var(shuffled_val_inds);
            end
        end
        
        switch opts.btsp_variables{s}
            case 'Y'
                PID_btsp = calculate_PID(shuffled_var,X1_b,X2_b,n_trials);
            case 'X1'
                PID_btsp = calculate_PID(Y_b,shuffled_var,X2_b,n_trials);
            case 'X2'
                PID_btsp = calculate_PID(Y_b,X1_b,shuffled_var,n_trials);
        end
        
        btsp_atoms{s}.uniqueX1(b)      = PID_btsp.uniqueX1;
        btsp_atoms{s}.uniqueX2(b)      = PID_btsp.uniqueX2;
        btsp_atoms{s}.shared(b)        = PID_btsp.shared;
        btsp_atoms{s}.complementary(b) = PID_btsp.complementary;
    end
end

% return results
offset_outputs = 1;
if opts.btsp > 0
    varargout{1} = btsp_atoms;
end

end
%% Function that calculated PID atoms for a single shuffling
function [PID] = calculate_PID(Y_b,X1_b,X2_b,...
    n_trials)

% calculate PID
PID_tmp = calculate_pid(X2_b,Y_b,X1_b,n_trials,n_trials);
PID.shared = PID_tmp(1);
PID.uniqueX1 = PID_tmp(2);
PID.uniqueX2 = PID_tmp(3);
PID.complementary = PID_tmp(4);
end