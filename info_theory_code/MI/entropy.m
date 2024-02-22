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
%%%
%%% see the documentation of [`binr`](tools/Binning/binr) function for more details on each of these binning strategies.
%%%
%%% - opts.btsp (optional, default: *opt.btsp = 0*): this field must be a (non-negative) scalar specifying how many bootstrap estimates to compute.These estimates are useful as an additional bias correction or for performing statistics on the entropy values.
%%%       Bootstrap estimates are performed by means of pairing X and Y at random and computing the entropy quantities for these random pairings; each estimate corresponds to a different random pairing configuration.
%%%       See the examples below for additional information on how to use this option.
%%% - opts.verbose (optional, default *opt.verbose = true*): if this field exists and is set to true a summary of the selected options is displayed and additional checks are performed on the input variables. No warnings are displayed unless this options is enabled. This feature is useful to check whether INFORMATION is being called correctly. It is therefore highly reccomended for new users or when first running of the program with new input options. However, keep in mind that these checks drammatically increases computation time and are thus not reccommended for computationally intensive session. If a custom bias correction function is called (see "BUILDING AND CALLING CUSTOM BIAS CORRECTION FUNCTIONS" below) tests are performed on the invoked function to check whether it satisfied the requirements of the toolbox.
%%%   
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description       |
%%% |---------|-------------------|
%%% | 'HX'    | $H(X)$            |
%%% | 'HXY'   | $H(X|Y)$          |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
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

allowed_outputs_list = ["HX","HXY"];
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
str_bin_methodY = string(opts.bin_methodY);

assert(length(X(1,:)) == length(Y(1,:)),...
    "Dimensions of X should be nDimensionsX X nTrials, while Y should be of shape nDimensionsY X nTrials");
if ~strcmp(opts.bin_methodX, 'none')
    X = binr(X, opts.n_binsX, opts.bin_methodX);
end
if ~strcmp(opts.bin_methodY, 'none')
    Y = binr(Y, opts.n_binsY, opts.bin_methodY);
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

pars.isBtsp  = 0;
[HX, HXY(:,1)] = direct_method(X, pars);

Py = pars.Nt ./ pars.totNt;
% Bootstrap
if pars.btsp>0
    
    % Setting bootstrap "do options". No bootstrap is computed for HX, HlX,
    % HiXY (this quantity is only computed for test purposes) and HshX:
    pars.doHX    = false;
    pars.doHXY   = pars.doHXYbs;
    
    % Signaling that the X matrix needs to be bootstrapped:
    pars.isBtsp  = 1;

    for k=2:pars.btsp+1

        % Randperm (inlining for speed):
        [ignore, randIndxes] = sort(rand(pars.totNt,1));
        % Randomly assigning trials to stimuli as defined by NT:
        X(:, pars.trials_indxes(randIndxes)) = X(:, pars.trials_indxes);

        [ignore2, HXY(:,k:k*pars.doHXYbs), ignore3, ignore4, ignore5, ignore6, ignore7, ignore8, ignore9] = ...
            direct_method(X, pars);
        
    end
    
    Py_bootstrap = Py(:, ones(pars.btsp+1,1));
    
end

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
