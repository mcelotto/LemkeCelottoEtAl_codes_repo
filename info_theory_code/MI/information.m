function outputs = information(X, Y, opts, outputsList)
%%% *function outputs = information(X, Y, opts, outputsList)*
%%%
%%% ### Description
%%% Compute mutual information.
%%%
%%% ### Inputs:
%%% - *X*: *nDimensionsX X nTrials* array including input data (typically the neural response).
%%% - *Y*: *nDimensionsY X nTrials* array including input data (typically the stimulus).
%%% - *opts*: options structure (see further notes section for more details).
%%% - *outputs*: cell array of char arrays specifying the quantities that the function is computing. These are going to be returned, in the same order as specified in this list in *outputs*. See further notes section for more details on the available options.
%%%
%%% ### Outputs:
%%% - *outputs*: cell array of same length as *outputsList* returning the specified outputs.
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
%%% - opts.btsp (optional): this field must be a (non-negative) scalar specifying how many bootstrap estimates to compute.These estimates are useful as an additional bias correction or for performing statistics on the entropy values.
%%%       Bootstrap estimates are performed by means of pairing stimuli and responses at random and computing the entropy quantities for these random pairings; each estimate corresponds to a different random pairing configuration.
%%%       See the examples below for additional information on how to use this option.
%%%       DEFAULT: 0.
%%% - opts.verbose (optional): if this field exists and is set to true a summary of the selected options is displayed and additional checks are performed on the input variables. No warnings are displayed unless this options is enabled. This feature is useful to check whether INFORMATION is being called correctly. It is therefore highly reccomended for new users or when first running of the program with new input options. However, keep in mind that these checks drammatically increases computation time and are thus not reccommended for computationally intensive session.
%%%       DEFAULT: true.
%%%   
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description       | Expression (in terms of ENTROPY output options)|
%%% |---------|-------------------|------------------------------------------------|
%%% | 'I'     | I(X;Y)            | I     = HX - HXY                               |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
%%%
%%% #### Output options with bootstrap
%%% Bootstrap estimates are returned to the user in the form of an array of length Opt.btsp concatenated to the actual information estimate. For example, suppose that I is computed with Opt.btsp =
%%% 20. In this case the output corresponding to I will be an array of length 21: The first element is the actual entropy estimate while the remaining 20 elements correspond to 20 distinct bootstrap estimates (see Fig. 2).
%%%
%%%          Actual | Bootstrap
%%%        Estimate | Estimates
%%%                 |
%%%     Index:    1 | 2   3   4      Opt.btsp+1 = 21
%%%             ----|------------- ... -----
%%%             | x | x | x | x |      | x |
%%%             ----|------------- ... -----
%%%                 |
%%%
%   Copyright (C) 2009 Cesare Magri
%   Version: 1.0.5
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

if nargin < 4
    error("information.m: not enough input arguments. See `help information` for usage info");
end

allowed_outputs_list = ["I"];
for i = 1:length(outputsList)
    if ~any(matches(allowed_outputs_list,string(outputsList{i})))
        error("Specified output quantity: `" + string(outputsList{i})...
            + "` not recognized. Valid quantities are: '"...
            + join(allowed_outputs_list, "' '") + "'")
    end
end

whereI     = strcmpi(outputsList, 'i');
doI     = any(whereI);

% Checks ------------------------------------------------------------------
if any(isnan(X))
    error("X contains NaNs. Aborting.")
end
if any(isnan(Y))
    error("Y contains NaNs. Aborting.")
end

specifiedOutputOptsVec = [doI];
NspecifiedOutputOpts = sum(specifiedOutputOptsVec);
if NspecifiedOutputOpts~=length(outputsList)
    msg = 'Unknown selection or repeated option in output list.';
    error(msg);
end

% Default verbose value:
if ~isfield(opts, 'verbose')
    opts.verbose = true;
end

allOutputOpts = {'HX' 'HXY'};
positionInOuputOptsList = 0;

% What needs to be computed by ENTTROPY -----------------------------------
% Need to compute H(X)? 
doHX = false;
if doI
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHX = positionInOuputOptsList;
    doHX = true;
end

% Need to compute H(X|Y)?
doHXY = false;
if doI
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHXY = positionInOuputOptsList;
    doHXY = true;
end


% Computing information theoretic quantities ------------------------------
outputOptsList = allOutputOpts([doHX doHXY]);

H = entropy(X, Y, opts, outputOptsList);

% Assigning output --------------------------------------------------------
outputs = cell(length(outputsList),1);

% I = HX - HXY
if doI
    outputs(whereI) = {H{whereHX} - H{whereHXY}};
end
