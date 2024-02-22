function outputs = information(X, Y, opts, outputsList)
%%% *function outputs = information(X, Y, opts, outputsList)*
%%%
%%% ### Description
%%% Compute mutual information and information-like quantities.
%%%
%%% ### Inputs:
%%% - *X*: *nDimensionsX X nTrials* array including input data (typically the neural response).
%%% - *Y*: *nDimensionsY X nTrials* array including input data (typically the stimulus). In principle the values of Y can also represent a multi-dimensional neural response (if the MI between two neuronal populations is of interest). While the MI calculated is symmetrical to swapping X and Y, the user should be aware that the information breakdown quantities are calculated only on X. 
%%% - *opts*: options structure (see further notes section for more details).
%%% - *outputsList*: cell array of char arrays specifying the quantities that the function is computing. These are going to be returned, in the same order as specified in this list in *outputs*. See further notes section for more details on the available options.
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
%%% | `'eqspace'`  | Equispaced binning                 |
%%% | `'ceqspace'` | Centered equispaced spaced binning |
%%% | `'geqspace'` | Gaussian equispaced spaced binning |
%%%
%%% see the documentation of [binr](tools/Binning/binr) function for more details on each of these binning strategies.
%%%
%%% - opts.method: this field specifies which estimation method to use and can be one of the following strings:
%%%
%%% | Option | Description     |
%%% |--------|-----------------|
%%% | `'dr'` | Direct method   |
%%% | `'gs'` | Gaussian method |
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
%%% - opts.btsp (optional): this field must be a (non-negative) scalar specifying how many bootstrap estimates to compute.These estimates are useful as an additional bias correction or for performing statistics on the entropy values.
%%%       Bootstrap estimates are performed by means of pairing stimuli and responses at random and computing the entropy quantities for these random pairings; each estimate corresponds to a different random pairing configuration.
%%%       See the examples below for additional information on how to use this option.
%%%       DEFAULT: 0.
%%% - opt.xtrp (optional, default: Opt.xtrp = 0): this field must be a (non-negative) scalar specifying how many  iterations repetitions of the extrapolation procedure should be performed. Extrapolations procedure (such as the quadratic extrapolation) perform bias estimation by computing entropy values on sub-groups of the available trials. These subgroups are created randomly. The xtrp option allows to average the extrapolation values over as many different random partitions as specified by the parameter.
%%% - opts.verbose (optional): if this field exists and is set to true a summary of the selected options is displayed and additional checks are performed on the input variables. No warnings are displayed unless this options is enabled. This feature is useful to check whether INFORMATION is being called correctly. It is therefore highly reccomended for new users or when first running of the program with new input options. However, keep in mind that these checks drammatically increases computation time and are thus not reccommended for computationally intensive session.
%%%       DEFAULT: true.
%%%   
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description       | Expression (in terms of ENTROPY output options)|
%%% |---------|-------------------|------------------------------------------------|
%%% | 'I'     | I(X;Y)            | I     = HX - HXY                               |
%%% | 'Ish'   | I(X;Y) shuffle    | Ish   = HX - HiXY + HshXY - HXY                |
%%% | 'IX'    | I(X)              | IX    = HlX - HX                               |
%%% | 'ILIN'  | I_lin(Y;X)        | ILIN  = HlX - HiXY                             |
%%% | 'SYN'   | Syn               | SYN   = HX - HXY - HlX + HiXY                  |
%%% | 'SYNsh' | Syn shuffle       | SYNsh = HX + HshXY - HXY - HlX                 |
%%% | 'ISS'   | I_sig_sim         | ISS   = HiX - HlX                              |
%%% | 'IC'    | I_cor             | IC    = HX - HXY + HiXY - HiX                  |
%%% | 'ICsh'  | I_cor shuffle     | ICsh  = HX + HshXY - HXY - HiX                 |
%%% | 'ICI'   | I_cor_ind         | ICI   = ChiX - HiX                             |
%%% | 'ICD'   | I_cor_dep         | ICD   = HX - HXY - ChiX + HiXY                 |
%%% | 'ICDsh' | I_cor_dep shuffle | ICDsh = HX + HshXY - HXY - ChiX                |
%%% | 'ILB1'  | I_LB1             | ILB1  = HX - HiXY                              |
%%% | 'ILB2'  | I_LB2             | ILB2  = ChiX - HiXY                            |
%%% | 'G1'    | I(X1,X2;Y) - I(X1;Y)                                               |      
%%% | 'G2'    | I(X1,X2;Y) - I(X1;Y)                                               |
%%% | 'GMAX'  | I(X1,X2;Y) - max(I(X1;Y), I(X2;Y))                                 |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
%%%  The function also computes IXY and IXYsh according to the following equations:
%%%
%%% | Option  | Description       | Expression (in terms of ENTROPY output options) |
%%% |---------|-------------------|-------------------------------------------------|
%%% | 'IXY'   | I(X|Y)            | IXY   = HiXY - HXY                              |
%%% | 'IXYsh' | I(X|Y) shuffle    | IXYsh = HshXY - HXY                             |
%%%
%%% However this quantity is meaningful only for *nDimensions=2*. 
%%% IMPORTANT: Not all combinations of method, bias and output options are possible. For example, bias correction 'pt' can only be used together  with method 'dr'. The allowed combinations of method, bias and output options are summarized in the following tables:
%%%
%%% #### Allowed outputs/bias combinations for direct method estimation
%%%
%%% Legend:
%%% - **X**: combination available
%%% - **-**: combination NOT permitted
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%%% |---------|---------|-------|-------|-------|
%%% | 'I'     |    X    |   X   |   X   |   -   |
%%% | 'Ish'   |    X    |   X   |   X   |   -   |
%%% | 'IX'    |    X    |   X   |   X   |   -   |
%%% | 'ILIN'  |    X    |   X   |   X   |   -   |
%%% | 'SYN'   |    X    |   X   |   X   |   -   |
%%% | 'SYNsh' |    X    |   X   |   X   |   -   |
%%% | 'ISS'   |    X    |   X   |   -   |   -   |
%%% | 'IC'    |    X    |   X   |   -   |   -   |
%%% | 'ICsh'  |    X    |   X   |   -   |   -   |
%%% | 'ICI'   |    X    |   X   |   -   |   -   |
%%% | 'ICD'   |    X    |   X   |   -   |   -   |
%%% | 'ICDsh' |    X    |   X   |   -   |   -   |
%%% | 'ILB1'  |    X    |   X   |   X   |   -   |
%%% | 'ILB2'  |    X    |   X   |   -   |   -   |
%%%
%%% #### Allowed outputs/bias combinations for Gaussian method estimation
%%% 
%%% Legend:  
%%% - **X**: combination available
%%% - **-**: combination NOT permitted
%%% - **n.r.**: combination available but not recommended
%%% - **NaN** : NaN returned
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%%% |---------|---------|-------|-------|-------|
%%% | 'I'     |    X    |  n.r. |   -   |   X   |
%%% | 'Ish'   |    X    |  n.r. |   -   |   X   |
%%% | 'IX'    |    X    |  n.r. |   -   |   X   |
%%% | 'ILIN'  |    X    |  n.r. |   -   |   X   |
%%% | 'SYN'   |    X    |  n.r. |   -   |   X   |
%%% | 'SYNsh' |    X    |  n.r. |   -   |   X   |
%%% | 'ISS'   |   NaN   |  NaN  |   -   |  NaN  |
%%% | 'IC'    |   NaN   |  NaN  |   -   |  NaN  |
%%% | 'ICsh'  |   NaN   |  NaN  |   -   |  NaN  |
%%% | 'ICI'   |   NaN   |  NaN  |   -   |  NaN  |
%%% | 'ICD'   |   NaN   |  NaN  |   -   |  NaN  |
%%% | 'ICDsh' |   NaN   |  NaN  |   -   |  NaN  |
%%% | 'ILB1'  |    X    |  n.r. |   -   |   X   |
%%% | 'ILB2'  |   NaN   |  NaN  |   -   |  NaN  |
%%%
%%% #### Allowed IXY/bias combinations for direct method estimation
%%%
%%% Legend:
%%% - **X**: combination available
%%% - **-**: combination NOT permitted
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%%% |---------|---------|-------|-------|-------|
%%% | 'IXY'   |    X    |   X   |   X   |   -   |
%%% | 'IXYsh' |    X    |   X   |   X   |   -   |
%%%
%%% #### Allowed IXY/bias combinations for Gaussian method estimation
%%%
%%% Legend:
%%% - **X**: combination available
%%% - **-**: combination NOT permitted
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'gsb' |
%%% |---------|---------|-------|-------|-------|
%%% | 'IXY'   |    X    |  n.r. |   -   |   X   |
%%% | 'IXYsh' |    X    |  n.r. |   -   |   X   |
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
%%% [I, Ish] = information(X, Y, opts, {'I', 'Ish'});
%%% ```
%%% is faster than
%%% ```
%%% I = information(X, Y, opts, {'I'});
%%% Ish = information(X, Y, opts, {'Ish'});
%%%```
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

allowed_outputs_list = ["I","Ish","IX","ILIN","SYN","SYNsh","ISS","IC",...
    "ICsh","ICI","ICD","ICDsh","ILB1","ILB2"];
for i = 1:length(outputsList)
    if ~any(matches(allowed_outputs_list,string(outputsList{i})))
        error("Specified output quantity: `" + string(outputsList{i})...
            + "` not recognized. Valid quantities are: '"...
            + join(allowed_outputs_list, "' '") + "'")
    end
end

whereI     = strcmpi(outputsList, 'i');
whereIsh   = strcmpi(outputsList, 'ish');
whereIX    = strcmpi(outputsList, 'ix');
whereIXY   = strcmpi(outputsList, 'ixy');
whereIXYsh = strcmpi(outputsList, 'ixysh');
whereILIN  = strcmpi(outputsList, 'ilin');
whereSYN   = strcmpi(outputsList, 'syn');
whereSYNsh = strcmpi(outputsList, 'synsh');
whereISS   = strcmpi(outputsList, 'iss');
whereIC    = strcmpi(outputsList, 'ic');
whereICsh  = strcmpi(outputsList, 'icsh');
whereICI   = strcmpi(outputsList, 'ici');
whereICD   = strcmpi(outputsList, 'icd');
whereICDsh = strcmpi(outputsList, 'icdsh');
whereILB1  = strcmpi(outputsList, 'ilb1');
whereILB2  = strcmpi(outputsList, 'ilb2');
    
doI     = any(whereI);
doIsh   = any(whereIsh);
doIX    = any(whereIX);
doIXY   = any(whereIXY);
doIXYsh = any(whereIXYsh);
doILIN  = any(whereILIN);
doSYN   = any(whereSYN);
doSYNsh = any(whereSYNsh);
doISS   = any(whereISS);
doIC    = any(whereIC);
doICsh  = any(whereICsh);
doICI   = any(whereICI);
doICD   = any(whereICD);
doICDsh = any(whereICDsh);
doILB1  = any(whereILB1);
doILB2  = any(whereILB2);

% Checks ------------------------------------------------------------------
if any(isnan(X))
    error("X contains NaNs. Aborting.")
end
if any(isnan(Y))
    error("Y contains NaNs. Aborting.")
end

specifiedOutputOptsVec = ...
    [doI doIsh doIX doIXY doIXYsh doILIN doSYN doSYNsh doISS doIC doICsh doICI doICD doICDsh doILB1 doILB2];
NspecifiedOutputOpts = sum(specifiedOutputOptsVec);
if NspecifiedOutputOpts~=length(outputsList)
    msg = 'Unknown selection or repeated option in output list.';
    error(msg);
end

% Restrictions on possible combinantions ----------------------------------
if strcmpi(opts.method, 'dr')
    % Can't apply bias-correction gsb with method dr:
    if strcmpi(opts.bias, 'gsb')
        msg = 'Bias correction ''gsb'' can only be used in conjunction with method ''gs''.';
        error('Information:drMethodAndGsbBias', msg);
    end

    % Can't compute ISS, IC, ICsh, ICI, ICD, ICDsh or ILB2 for
    % bias-correction pt
    if strcmpi(opts.bias, 'pt') && (doISS || doIC || doICsh || doICI || doICD || doICDsh || doILB2)
        msg = 'One or more of the selected output options are not available for bias correction ''pt''.';
        error(msg);
    end
end

% Default verbose value:
if ~isfield(opts, 'verbose')
    opts.verbose = false;
end

isGaussianMethod = false;
if strcmpi(opts.method, 'gs')
    isGaussianMethod = true;
    
    % Gaussian and QE is not recommended:
    if opts.verbose && strcmpi(opts.bias, 'qe')
        msg = 'Usage of bias correction ''qe'' in conjunction with gaussian method is not recommended.';
        warning(msg);
    end
    
    % Gaussian and PT is not allowed:
    if strcmpi(opts.bias, 'pt')
        msg = 'Bias correction ''pt'' can only be used in conjunction with method ''dr''.';
        error(msg);
    end
end

allOutputOpts = {'HX' 'HXY' 'HlX' 'HiX' 'HiXY' 'ChiX' 'HshXY'};
positionInOuputOptsList = 0;

% What needs to be computed by ENTTROPY -----------------------------------
% Need to compute H(X)? 
doHX = false;
if doI || doIsh || doIX || doSYN || doSYNsh || doIC || doICsh || doICD || doICDsh || doILB1
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHX = positionInOuputOptsList;
    doHX = true;
end

% Need to compute H(X|Y)?
doHXY = false;
if doI || doIsh || doIXY || doIXYsh || doSYN || doSYNsh || doIC || doICsh || doICD || doICDsh
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHXY = positionInOuputOptsList;
    doHXY = true;
end

% Need to compute H_lin(X)?
doHlX = false;
if doIX || doSYN || doILIN || doSYNsh || doISS
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHlX = positionInOuputOptsList;
    doHlX = true;
end

% Need to compute H_ind(X)?
doHiX = false;
if (doISS || doIC || doICsh || doICI) && ~isGaussianMethod
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHiX = positionInOuputOptsList;
    doHiX = true;
end

% Need to compute H_ind(X|Y)?
doHiXY = false;
if doIsh || doIXY || doILIN || doSYN || doIC || doICD || doILB1 || doILB2
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHiXY = positionInOuputOptsList;
    doHiXY = true;
end

% Need to compute Chi(X)?
doChiX = false;
if (doICI || doICD || doICDsh || doILB2) && ~isGaussianMethod
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereChiX = positionInOuputOptsList;
    doChiX = true;
end

% Need to compute H_sh(X|Y)?
doHshXY = false;
if doIsh || doIXYsh || doSYNsh || doICsh || doICDsh
    positionInOuputOptsList = positionInOuputOptsList + 1;
    whereHshXY = positionInOuputOptsList;
    doHshXY = true;
end

% Computing information theoretic quantities ------------------------------
outputOptsList = allOutputOpts([doHX doHXY doHlX doHiX doHiXY doChiX doHshXY]);

% H = cell(positionInOuputOptsList, 1);
H = entropy(X, Y, opts, outputOptsList);

% Assigning output --------------------------------------------------------
outputs = cell(length(outputsList),1);

% I = HX - HXY
if doI
    outputs(whereI) = {H{whereHX} - H{whereHXY}};
end

% Ish = HX - HiXY + HshXY - HXY
if doIsh
    outputs(whereIsh) = {H{whereHX} - H{whereHiXY} + H{whereHshXY} - H{whereHXY}};
end

% IX = HlX - HX
if doIX
    outputs(whereIX) = {H{whereHlX} - H{whereHX}};
end

% IXY = HiXY - HXY
if doIXY
    outputs(whereIXY) = {H{whereHiXY} - H{whereHXY}};
end

% IXYsh = HshXY - HXY
if doIXYsh
    outputs(whereIXYsh) = {H{whereHshXY} - H{whereHXY}};
end

% ILIN = HlX - HiXY
if doILIN
    outputs(whereILIN) = {H{whereHlX} - H{whereHiXY}};
end

% SYN = HX - HXY - HlX + HiXY
if doSYN
    outputs(whereSYN) = {H{whereHX} - H{whereHXY} - H{whereHlX} + H{whereHiXY}};
end

% SYNsh = HX + HshXY - HXY - HlX
if doSYNsh
    outputs(whereSYNsh) = {H{whereHX} + H{whereHshXY} - H{whereHXY} - H{whereHlX}};
end

% ISS = HiX - HlX
if doISS
    if ~isGaussianMethod
        outputs(whereISS) = {H{whereHiX} - H{whereHlX}};
    else
        outputs(whereISS) = {NaN};
    end
end

% IC = HX - HXY + HiXY - HiX
if doIC
    if ~isGaussianMethod
        outputs(whereIC) = {H{whereHX} - H{whereHXY} + H{whereHiXY} - H{whereHiX}};
    else
        outputs(whereIC) = {NaN};
    end
end

% ICsh  = HX + HshXY - HXY - HiX
if doICsh
    if ~isGaussianMethod
        outputs(whereICsh) = {H{whereHX} + H{whereHshXY} - H{whereHXY} - H{whereHiX}};
    else
        outputs(whereICsh) = {NaN};
    end
end

% ICI= ChiX - HiX
if doICI
    if ~isGaussianMethod
        outputs(whereICI) = {H{whereChiX} - H{whereHiX}};
    else
        outputs(whereICI) = {NaN};
    end
end

% ICD = HX - HXY - ChiX + HiXY
if doICD
    if ~isGaussianMethod
        outputs(whereICD) = {H{whereHX} - H{whereHXY} - H{whereChiX} + H{whereHiXY}};
    else
        outputs(whereICD) = {NaN};
    end
end

% ICDsh = HX + HshXY - HXY - ChiX
if doICDsh
    if ~isGaussianMethod
        outputs(whereICDsh) = {H{whereHX} + H{whereHshXY} - H{whereHXY} - H{whereChiX}};
    else
        outputs(whereICDsh) = {NaN};
    end
end

% ILB1 = HX - HiXY
if doILB1
    outputs(whereILB1) = {H{whereHX} - H{whereHiXY}};
end

% ILB2 = ChiX - HiXY
if doILB2
    if ~isGaussianMethod
        outputs(whereILB2) = {H{whereChiX} - H{whereHiXY}};
    else
        outputs(whereILB2) = {NaN};
    end
end
