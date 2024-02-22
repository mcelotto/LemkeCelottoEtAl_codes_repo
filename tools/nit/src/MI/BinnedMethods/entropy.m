function entropies = entropy(X, Y, opts, outputsList)
%%% *function entropies = entropy(X, Y, opts, outputsList)*
%%%
%%% ### Description
%%% Compute entropy and entropy-like quantities.
%%%
%%% ### Inputs:
%%% - *X*: *nDimensionsX X nTrials* array including input data (typically the neural response).
%%% - *Y*: *nDimensionsY X nTrials* array including input data (typically the stimulus). In principle the values of Y can also represent a multi-dimensional neural response (if the MI between two neuronal populations is of interest). While the MI calculated is symmetrical to swapping X and Y, the user should be aware that the information breakdown quantities are calculated only on X. 
%%% - *opts*: options structure (see further notes section for more details).
%%% - *outputsList*: cell array of char arrays of strings specifying the quantities that the function is computing. These are going to be returned, in the same order as specified in this list in *entropies*. See further notes section for more details on the available options.
%%%
%%% ### Outputs:
%%% - *entropies*: cell array of same length as *outputsList* returning the specified outputs.
%%%
%%% ### Further notes
%%% #### The options structure
%%% The options structure can include any the following fields:
%%%
%%% - opts.n_binsX: number of bins array to be used on X (can be a scalar or a vector of `size(nb) = size(X,1)`). If a scalar is specified then all dimensions are binned using the same number of bins. Alternatively, each dimension is binned using its specific number of bins. 
%%% - opts.n_binsY: number of bins array to be used on Y (can be a scalar or a vector of `size(nb) = size(X,1)`). If a scalar is specified then all dimensions are binned using the same number of bins. Alternatively, each dimension is binned using its specific number of bins. 
%%% - opts.bin_methodX: binning method for X (can be a scalar or a cell array of `size(method) = size(X,1)`). If a scalar is specified then all dimensions are binned using the same method. Alternatively, each dimension is binned using its specific method. If not specified it is assumed that no binning is necessary and X is already discrete.
%%% - opts.bin_methodY: binning method for Y (can be a scalar or a cell array of `size(method) = size(Y,1)`). If a scalar is specified then all dimensions are binned using the same method. Alternatively, each dimension is binned using its specific method. If not specified it is assumed that no binning is necessary and Y is already discrete.
%%% 
%%% Possible options for binning methods are:
%%%
%%% | Option       | Description                        |
%%% |--------------|------------------------------------|
%%% | `'none'`     | No binning                         |
%%% | `'eqpop'`    | Equally populated binning          |
%%% | `'eqspace'`  | Equispaced binning                 |
%%% | `'ceqspace'` | Centered equispaced spaced binning |
%%% | `'geqspace'` | Gaussian equispaced spaced binning |
%%%
%%% see the documentation of [`binr`](tools/Binning/binr) function for more details on each of these binning strategies.
%%%
%%% - opts.method: this field specifies which estimation method to use and can be one of the following strings:
%%%
%%% | Option | Description     |
%%% |--------|-----------------|
%%% | `'dr'` | Direct method   |
%%% | `'gs'` | Gaussian method |
%%%
%%% The direct method requires X values to be discretized into non-negative integer values. This is handled internally to the code through the *opts.bin_methodX* option. It can also be done by the user prior to the call to this function but the user should be aware that the binned values should be strincly non-negative. **Failing to properly discretizing the X array will result in Matlab crashing.**
%%%
%%% - opts.bias: this field specifies the bias correction procedure. It can be one of the following strings:
%%%
%%% | Option    | Description             |
%%% |-----------|-------------------------|
%%% | `'qe'`    | Quadratic extrapolation |
%%% | `'pt'`    | Panzeri & Treves 1996   |
%%% | `'gsb'`   | Gaussian bias           |
%%% | `'naive'` | Biased naive estimates  |
%%%
%%% - opts.btsp (optional, default: *opt.btsp = 0*): this field must be a (non-negative) scalar specifying how many bootstrap estimates to compute.These estimates are useful as an additional bias correction or for performing statistics on the entropy values.
%%%       Bootstrap estimates are performed by means of pairing X and Y at random and computing the entropy quantities for these random pairings; each estimate corresponds to a different random pairing configuration.
%%%       See the examples below for additional information on how to use this option.
%%% - opt.xtrp (optional, default: *opt.xtrp = 0*): this field must be a (non-negative) scalar specifying how many  iterations repetitions of the extrapolation procedure should be performed. Extrapolations procedure (such as the quadratic extrapolation) perform bias estimation by computing entropy values on sub-groups of the available trials. These subgroups are created randomly. The xtrp option allows to average the extrapolation values over as many different random partitions as specified by the parameter.
%%% - opts.verbose (optional, default *opt.verbose = true*): if this field exists and is set to true a summary of the selected options is displayed and additional checks are performed on the input variables. No warnings are displayed unless this options is enabled. This feature is useful to check whether INFORMATION is being called correctly. It is therefore highly reccomended for new users or when first running of the program with new input options. However, keep in mind that these checks drammatically increases computation time and are thus not reccommended for computationally intensive session. If a custom bias correction function is called (see "BUILDING AND CALLING CUSTOM BIAS CORRECTION FUNCTIONS" below) tests are performed on the invoked function to check whether it satisfied the requirements of the toolbox.
%%%   
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description       |
%%% |---------|-------------------|
%%% | 'HX'    | $H(X)$            |
%%% | 'HXY'   | $H(X|Y)$          |
%%% | 'HlX'   | $H_{lin}(X)$      |
%%% | 'HiX'   | $H_{ind}(X)$      |
%%% | 'HiXY'  | $H_{ind}(X|Y)$    |
%%% | 'ChiX'  | $Chi(X)$          |
%%% | 'HshX'  | $H_{sh}(X)$       |
%%% | 'HshXY' | $H_{sh}(X|Y)$     |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
%%% IMPORTANT: Not all combinations of method, bias and output options are possible. For example, bias correction 'pt' can only be used together  with method 'dr'. The allowed combinations of method, bias and output options are summarized in the following tables: 
%%%   
%%% #### Allowed outputs/bias combinations for direct method estimation
%%% Legend:
%%% - **X**: combination available
%%% - **n**: naive estimate returned
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'gsb' | 'user defined' |
%%% |---------|---------|-------|-------|-------|-----------------
%%% | 'HX'    |    X    |   X   |   X   |   n   |       X        |
%%% | 'HXY'   |    X    |   X   |   X   |   n   |       X        |
%%% | 'HlX'   |    X    |   X   |   X   |   n   |       X        |
%%% | 'HiX'   |    X    |   X   |   n   |   n   |       n        |
%%% | 'HiXY'  |    X    |   X   |   X   |   n   |       n        |
%%% | 'ChiX'  |    X    |   X   |   n   |   n   |       n        |
%%% | 'HshX'  |    X    |   X   |   X   |   n   |       X        |
%%% | 'HshXY' |    X    |   X   |   X   |   n   |       X        |
%%% 
%%% #### Allowed outputs/bias combinations for Gaussian method estimation
%%% Legend:  
%%% - **X**: combination available
%%% - **n**: naive estimate returned
%%% - **0** : zero returned
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'gsb' | 'user defined' |
%%% |---------|---------|-------|-------|-------|-----------------
%%% | 'HX'    |    X    |   X   |   n   |   X   |       n        |
%%% | 'HXY'   |    X    |   X   |   n   |   X   |       n        |
%%% | 'HlX'   |    X    |   X   |   n   |   X   |       n        |
%%% | 'HiX'   |    0    |   0   |   0   |   0   |       n        |
%%% | 'HiXY'  |    X    |   X   |   n   |   X   |       n        |
%%% | 'ChiX'  |    0    |   0   |   0   |   0   |       n        |
%%% | 'HshX'  |    X    |   X   |   n   |   X   |       n        |
%%% | 'HshXY' |    X    |   X   |   n   |   X   |       n        |
%%%
%%% #### Output options with bootstrap
%%% Bootstrap estimates make sense (and can thus be computed) only for
%%% the following output quantities: 'HXY', 'HiX', 'HiXY', 'ChiX' and
%%% HshXs'.
%%%
%%% If the user only specifies the number of bootstrap repetitions
%%% (through the BTSP parameter) bootstrap estimates are computed for
%%% any of the above quantities which appears in the output list.
%%%
%%% However users are given the opportunit to precisely select which
%%% quantities to compute bootstrap for by means of appending 'bs' to
%%% the output string name of an output quantitie as follows:
%%%
%%% | Quantity | Bootstrapped quantity  |
%%% |----------|------------------------|
%%% | 'HiX'    | 'HiXbs'                |
%%% | 'HiXY'   | 'HiXYbs'               |
%%% | 'ChiX'   | 'ChiXbs'               |
%%% | 'HshXY'  | 'HshXYbs'              |
%%%
%%%  Carefully selecting for which quantities to compute bootstrap can
%%%  greately decrease computation times. For example, bootstrap
%%%  estimates of HiX and ChiX are often useless although their
%%%  computation can be very time consuming.
%%%
%%%  Bootstrap estimates are returned to the user in the form of an
%%%  array of length Opt.btsp concatenated to the actual entropy
%%%  estimate. For example, suppose that HXY is computed with Opt.btsp =
%%%  20. In this case the output corresponding to HXY will be an array
%%%  of length 21: The first element is the actual entropy estimate
%%%  while the remaining 20 elements correspond to 20 distinct bootstrap
%%%  estimates (see Fig. 2).
%%%
%%%
%%%                Actual | Bootstrap
%%%              Estimate | Estimates
%%%                       |
%%%           Index:    1 | 2   3   4      Opt.btsp+1 = 21
%%%                   ----|------------- ... -----
%%%                   | x | x | x | x |      | x |
%%%                   ----|------------- ... -----
%%%                       |
%%%
%%% ### Building and calling custom bias correction functions
%%%
%%% In the direct method each value, P(x), of a probability distribution is
%%% estimated using the normalized count, C(x), of occurrence of x as 
%%% $P(x) = C(x)/N$.
%%%
%%% C contains all the information regarding the probability distribution
%%% P(x) and thus also all the parameters necessary for estimation of the
%%% bias of H(X). The same applies to the estimation of the other Y-
%%% unconditional probabilities, such as P_sh(x), and also to the Y-
%%% conditional probabilities: in the Y-conditional case bias
%%% correction is performed for each of the Ny probability distributions
%%% P(x|y) and then averaged across the values of Y.
%%%
%%% Users can define their custom bias correction methods as follows: let
%%% `custom_bias_corr_func` be the name of the user-defined correction routine
%%% (any valid MATLAB function name can be used, see Matlab documentation).
%%% This function must receive as input the array C of size *nDimensionsX x
%%% 1* described above and return the (positive) value to be **added** to the
%%% plugin entropy estimate:
%%%       
%%%         bias = custom_bias_corr_func(C)
%%%
%%% To call the custom bias correction functon simply pass the name of the
%%% function as a string in the Opt.bias field. For the above example we
%%% have
%%%       
%%%       opts.bias = 'custom_bias_corr_func';
%%%
%%% ### Remarks
%%% - Field-names in the option structure are case-sensitive
%%% - Ouput options are case INsensitive
%%% - It is more efficient to call INFORMATION with several output options rather than calling the function repeatedly. For example:
%%% ```
%%% [X, Y] = information(X, Y, opts, {'I', 'Ish'});
%%% ```
%%% is faster than
%%% ```
%%% X = information(X, Y, opts, {'I'});
%%% Y = information(X, Y, opts, {'Ish'});
%%%```
%%%
%   Copyright (C) 2011 Cesare Magri
%   Version 5
%
% -------
% LICENSE
% -------
% This software is distributed free under the condition that:
% 1. it shall not be incorporated in software that is subsequently sold;
% 2. the authorship of the software shall be acknowledged and the following article shall be properly cited in any publication that uses results generated by the software:
%
%      Magri C, Whittingstall K, Singh V, Logothetis NK, Panzeri S: A
%      toolbox for the fast information analysis of multiple-site LFP, EEG
%      and spike train recordings. BMC Neuroscience 2009 10(1):81;
%
% 3.  this notice shall remain in place in each source file.
%
% ----------
% DISCLAIMER
% ----------
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% Check inputs and set default values
if nargin < 4
    error("entropy.m: not enough input arguments. See `help entropy` for usage info");
end

allowed_outputs_list = ["HX","HXY","HlX","HiX","HiXY","ChiX","HshX","HshXY"];
for i = 1:length(outputsList)
    if ~any(matches(allowed_outputs_list,string(outputsList{i})))
        error("Specified output quantity: `" + string(outputsList{i})...
            + "` not recognized. Valid quantities are: '"...
            + join(allowed_outputs_list, "' '") + "'")
    end
end

if any(isnan(X))
    error("X contains NaNs. Aborting.")
end
if any(isnan(Y))
    error("Y contains NaNs. Aborting.")
end

% Default verbose value:
if ~isfield(opts, 'verbose')
    opts.verbose = true;
end

if ~any(matches(["qe", "pt", "naive", "gsb"],string(opts.bias)))
    % if specified bias correction method is not any of the default ones
    % check for existence of user-defined function with name opts.bias
    if ~exist(opts.bias)
        error("User defined bias correction method: `" + opts.bias...
            + "` cannot be found in the current MATLAB search path.")
    end
    if exist(opts.bias) ~= 2
        error("User defined bias correction method: `" + opts.bias...
            + "` could be found in the current MATLAB search path but does not appear to be a MATLAB function.")
    end    
end

if ~isfield(opts, 'bin_methodX')
    if  opts.verbose
        warning("No opts.bin_methodX specified. X is assumed to be already binned. Proceeding without binning.")
    end
    opts.bin_methodX = 'none';
end
if ~isfield(opts, 'bin_methodY')
    if opts.verbose
        warning("No opts.bin_methodY specified. Y is assumed to be already binned. Proceeding without binning.")
    end
    opts.bin_methodY = 'none';
end
if strcmp(opts.bin_methodX, 'none')
    opts.n_binsX = 0;
end
if strcmp(opts.bin_methodY, 'none')
    opts.n_binsY = 0;
end
if ~isfield(opts, 'n_binsX')
    error("Please specify number of bins for binning X through opts.n_binsX variable")
end
if ~isfield(opts, 'n_binsY')
    error("Please specify number of bins for binning Y through opts.n_binsY variable")
end
str_bin_methodX = string(opts.bin_methodX);
for i = 1:length(str_bin_methodX)
    if ~any(matches(["none", "eqpop", "eqspace", "ceqspace", "geqspace"],...
            str_bin_methodX(i)))
        % if specified binning method is not any of the default ones
        % check for existence of user-defined function with name opts.bin_methodX
        if ~exist(str_bin_methodX(i))
            error("User defined X binning method: `" + str_bin_methodX(i)...
                + "` cannot be found in the current MATLAB search path.")
        end
        if exist(str_bin_methodX(i)) ~= 2
            error("User defined X binning method: `" + str_bin_methodX(i)...
                + "` could be found in the current MATLAB search path but does not appear to be a MATLAB function.")
        end
    end
end
str_bin_methodY = string(opts.bin_methodY);
for i = 1:length(str_bin_methodY)
    if ~any(matches(["none", "eqpop", "eqspace", "ceqspace", "geqspace"],...
            str_bin_methodY(i)))
        % if specified binning method is not any of the default ones
        % check for existence of user-defined function with name opts.bin_methodY
        if ~exist(str_bin_methodY(i))
            error("User defined Y binning method: `" + str_bin_methodY(i)...
                + "` cannot be found in the current MATLAB search path.")
        end
        if exist(str_bin_methodY(i)) ~= 2
            error("User defined Y binning method: `" + str_bin_methodY(i)...
                + "` could be found in the current MATLAB search path but does not appear to be a MATLAB function.")
        end
    end
end
assert(length(X(1,:)) == length(Y(1,:)),...
    "Dimensions of X should be nDimensionsX X nTrials, while Y should be of shape nDimensionsY X nTrials");
if ~strcmp(opts.bin_methodX, 'none')
    X = binr(X, opts.n_binsX, opts.bin_methodX);
end
if ~strcmp(opts.bin_methodY, 'none')
    Y = binr(X, opts.n_binsY, opts.bin_methodY);
end
% if Y is multi-dimensional, project it to a 1d space
if length(Y(:,1)) > 1
    Y = map_Nd_array_to_1d(Y);
end
% build 3D X matrix
[X, nt] = buildx(Y, X);
opts.nt = nt;

% build params structure
XMatrixName = inputname(1);
if isfield(opts, 'pars')
    pars = opts.pars;
else
    pars = build_parameters_structure(X, opts, XMatrixName, outputsList);
end


% COMPUTING ENTROPIES =====================================================
HXY   = zeros(pars.Ns, (pars.btsp * pars.doHXYbs  )+1);
HlXY  = zeros(pars.Ns, (pars.btsp * pars.doHlXYbs )+1);
HiX   = zeros(1      , (pars.btsp * pars.doHiXbs  )+1);
ChiX  = zeros(1      , (pars.btsp * pars.doChiXbs )+1);
HshXY = zeros(pars.Ns, (pars.btsp * pars.doHshXYbs)+1);

pars.isBtsp  = 0;
if pars.biasCorrNum==1
    [HX, HXY(:,1), HlX, HlXY(:,1), HiX(1), HiXY, ChiX(1), HshX, HshXY(:,1)] = xtrploop(X, pars);
else
    [HX, HXY(:,1), HlX, HlXY(:,1), HiX(1), HiXY, ChiX(1), HshX, HshXY(:,1)] = pars.methodFunc(X, pars);
end

Py = pars.Nt ./ pars.totNt;
% Bootstrap
if pars.btsp>0
    
    % Setting bootstrap "do options". No bootstrap is computed for HX, HlX,
    % HiXY (this quantity is only computed for test purposes) and HshX:
    pars.doHX    = false;
    pars.doHXY   = pars.doHXYbs;
    pars.doHlX   = false;
    pars.doHlXY  = pars.doHlXYbs;
    pars.doHiX   = pars.doHiXbs;
    pars.doHiXY  = false;
    pars.doChiX  = pars.doChiXbs;
    pars.doHshX  = false;
    pars.doHshXY = pars.doHshXYbs;
    
    % Signaling that the X matrix needs to be bootstrapped:
    pars.isBtsp  = 1;

    % In the bootstrap loop we use the following trick: suppose we need to
    % compute bootstrap for HXY, then "pars.HXYbs = true" and
    %
    %   k:k*pars.doHXYbs = k:k*1 = k
    %
    % If bootstrap need not be computed we have "pars.HXYbs = false" and
    %
    %   k:k*pars.doHXYbs = k:k*0 = k:0 = []
    %
    % Thus, since the instruction "x([]) = value" does not produce any
    % effect, we have that for "pars.doHXYbs = true" the instruction
    %
    %   HXY(:,k:k*pars.doHXYbs) = value;
    %
    % will store the value on the right in HXY(:,k) while for "pars.doHXYbs
    % = false" the same instruction will produce no effect.
    
    for k=2:pars.btsp+1
        
        % For the direct method (methodNum==1) the bootstrap permutation
        % are currently performed outside the method function:
        if pars.methodNum~=3
            % Randperm (inlining for speed):
            [ignore, randIndxes] = sort(rand(pars.totNt,1));
            % Randomly assigning trials to stimuli as defined by NT:
            X(:, pars.trials_indxes(randIndxes)) = X(:, pars.trials_indxes);
        end
        
        if pars.biasCorrNum==1
            [ignore2, HXY(:,k:k*pars.doHXYbs), ignore3, HlXY(:,k:k*pars.doHlXYbs), HiX(:,k:k*pars.doHiXbs), ignore4, ChiX(:,k:k*pars.doChiXbs), ignore5, HshXY(:,k:k*pars.doHshXYbs)] = ...
                xtrploop(X, pars); %#ok<ASGLU>
        else
            [ignore2, HXY(:,k:k*pars.doHXYbs), ignore3, HlXY(:,k:k*pars.doHlXYbs), HiX(:,k:k*pars.doHiXbs), ignore4, ChiX(:,k:k*pars.doChiXbs), ignore5, HshXY(:,k:k*pars.doHshXYbs)] = ...
                pars.methodFunc(X, pars); %#ok<ASGLU>
        end
    end
    
    Py_bootstrap = Py(:, ones(pars.btsp+1,1));
    
end % ---------------------------------------------------------------------



% ASSIGNING OUTPUTS =======================================================
entropies = cell(pars.Noutput,1);

%Assigning HX output
entropies(pars.whereHX) = {HX};

% Assigning HXY output
if pars.doHXYbs && pars.btsp>0
    % Bootstrap case:
    entropies(pars.whereHXY)   = {sum(  HXY .* Py_bootstrap, 1)};
else
    % NON_bootstrap case:
    entropies(pars.whereHXY)   = {sum(  HXY(:,1) .* Py, 1)};
end

if pars.Nc==1
    entropies(pars.whereHlX)   = {HX};
    entropies(pars.whereHshX)  = {HX};
    
    if pars.doHlXYbs && pars.btsp>0
        entropies(pars.whereHlXY)  = {sum(  HXY .* Py_bootstrap, 1)};
    else
        entropies(pars.whereHlXY)  = {sum(  HXY(:,1) .* Py, 1)};
    end
    
    entropies(pars.whereHiXY)  = {sum(  HXY(:,1) .* Py, 1)};
    
    if pars.doHshXYbs && pars.btsp>0
        entropies(pars.whereHshXY) = {sum(  HXY .* Py_bootstrap, 1)};
    else
        entropies(pars.whereHshXY) = {sum(  HXY(:,1) .* Py, 1)};
    end
    
    switch pars.methodNum
        case 1
            entropies(pars.whereHiX)  = {HX};
            entropies(pars.whereChiX) = {HX};
            
        % In the gaussian method case we return NaN for consistency with
        % the NC>1 case:
        case 2
            entropies(pars.whereHiX)  = {NaN};
            entropies(pars.whereChiX) = {NaN};
    end
else 
    entropies(pars.whereHlX)   = {HlX};
    
    if pars.doHlXYbs && pars.btsp>0
        entropies(pars.whereHlXY)  = {sum( HlXY .* Py_bootstrap, 1)};
    else
        entropies(pars.whereHlXY)  = {sum( HlXY .* Py(:, ones((pars.btsp * pars.doHlXYbs) + 1, 1)), 1)};
    end
    entropies(pars.whereHiX)   = {HiX};
    
    % HiXY is simply a test quantity, we don't consider bootstrap for it.
    entropies(pars.whereHiXY)  = {sum( HiXY .* Py, 1)};
    
    entropies(pars.whereChiX)  = {ChiX};
    entropies(pars.whereHshX)  = {HshX};
    
    if pars.doHshXYbs && pars.btsp>0
        entropies(pars.whereHshXY) = {sum(HshXY .* Py_bootstrap, 1)};
    else
        entropies(pars.whereHshXY) = {sum(HshXY .* Py, 1)};
    end
end
