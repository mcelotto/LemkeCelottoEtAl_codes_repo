function transferentropies = transferentropy(X, Y, opts, outputsList)
%%% *function transferentropies = transferentropy(X, Y, opts, outputsList)*
%%%
%%% ### Description
%%% Computes transfer entropy from signal $X$ to signal $Y$.
%%%
%%% ### Inputs:
%%% - *X*: *nTimePts X nTrials* matrix containing the "causing" signal data matrix.
%%% - *Y*: *nTimePts X nTrials* matrix containing the "caused" signal data matrix.
%%% - *opts*: options structure.
%%% - *outputsList*: cell array of char arrays of strings specifying the quantities that the function is computing. These are going to be returned, in the same order as specified in this list in *entropies*. See further notes section for more details on the available options.
%%%
%%% ### Outputs:
%%% - *transferentropies*: cell array of same length as *outputsList* returning the specified outputs.
%%%
%%% ### Further notes
%%% ### The data matrices
%%%
%%% The function requires that the user specifies two types of data: the
%%% "causing" signal X and the "caused" signal Y. This toolbox allows to
%%% estimate transfer entropy, i.e. how much the knowledge of the past of X
%%% allows to incrase the predictability over the present of Y.
%%% 
%%% Both X and Y must be matrices of size *nTimePts X nTrials* where **nTimePts* is the
%%% number of points recorded in each trial while *nTrials* is the number of
%%% trials. Trials can be different recordings or the two signals or they
%%% can just reflect the way the user wishes to break the data into blocks
%%% for statistical testing (see [Bootstrapping](#bootstrapping) section.
%%% 
%%% It is important to remind that, in order to compute transfer entropy
%%% trials of X and corresponding trials of X must be simultaneously
%%% recorded.
%%%
%%% #### The options structure
%%% The function `transferentropy` is essentially a wrapper around function [`entropy`](MI/BinnedMethods/entropy).
%%% Consequently, `transferentropy` options include many of the options to
%%% `entropy` (exceptions are represented by the number-of-trials per
%%% stimulus opts.nt array and by the usage of the bootstrap option opts.btsp
%%% option, see [Bootstrapping](#bootstrapping) section.
%%%
%%% Below is a list of the required and optional option parameters: for a description of
%%% the inherited input options the user can refer to the
%%% [documentation](MI/BinnedMethods/entropy) of function `entropy`.
%%%
%%% ##### `transferentropy`-specific options
%%% - opts.taux (*required*): lags considered for the causing signal. Lags are specified as arrays of **strictly negative** integers.
%%% - opts.tuay (*required*): lags considered for the caused signal. Lags are specified as arrays of **strictly negative** integers.
%%% - opts.trperm (*optional*): number of random trial permutations. **Cannot be specified** in conjunction with *opts.btsp* option.
%%%
%%% ##### Options inherited from function `entropy`
%%% - opts.method (*required*): probability estimation method.
%%% - opts.bias (*required*): bias correction method.
%%% - opts.btsp (*optional*): number of bootstrap iterations. **Cannot be specified** in conjunction with *opts.trperm* option.
%%% - opts.verbose (*optional*): perform extensive checks on inputs.
%%%
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description                                                   |
%%% |---------|---------------------------------------------------------------|
%%% | 'TE'    | $X \to Y$ transfer entropy                                    |
%%% | 'TEsh'  | $X \to Y$ transfer entropy with shuffle correction            |
%%% | 'NTE'   | normalized $X \to Y$ transfer entropy                         |
%%% | 'NTEsh' | normalized $X \to Y$ transfer entropy with shuffle correction |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
%%%
%%% #### Trial-permutation and bootstrapped estimates
%%% For `opts.trperm > 0`, *trial-permutation* estimates are computed after 
%%% trials in X and Y have been randomly permuted. This can be useful for 
%%% bias performing statistics on the transfer-entropy values or to estimate
%%% the effect of a causing signal which is common to both X and Y. During
%%% trial-permutation the time-consistency of *X* and *Y* is kept, while
%%% trials are shuffled.
%%%
%%% *Bootstrapped estimates*, instead, destroy the time-consistency within
%%% trials by randomly shuffling time points within the same time lag
%%% considered. Bootstrap estimates are returned if `opts.btsp > 0`.
%%% 
%%% Trial-permutation and bootstrapping are *mutually exclusive*.
%%% Trial-permuted and bootstrap estimates can be returned for all
%%% [outputs](the-output-list) of the function.
%%%
%%% Bootstrap (or trial-permuted) estimates are returned to the user in
%%% the form of an array of length `opts.trperm` (or `opts.btsp`) concatenated
%%% to the actual transfer-entropy estimate. For example, suppose that
%%% TE is computed with Opt.btsp = 20. In this case the output
%%% corresponding to TE will be an array of length 21: the first
%%% element is the actual entropy estimate while the remaining 20
%%% elements correspond to 20 distinct bootstrap estimates.
%%%
%%%
%%%                Actual | bootstrap (trial-permuted)
%%%              Estimate | estimates
%%%                       |
%%%           Index:    1 | 2   3   4      Opt.btsp+1 = 21
%%%                   ----|------------- ... -----
%%%                   | x | x | x | x |      | x |
%%%                   ----|------------- ... -----
%%%                       |
%%%
%   -----------------------------------------------------------------------
%   TECHNICAL NOTE 1 - HOW TRANSFER-ENTROPY IS COMPUTED IN THIS ROUTINE
%   -----------------------------------------------------------------------
%   By the chain rule (see the "chain rule for information" in [1]) we have
%
%           I(Y1 Y2; X) = I(Y1; X) + I(Y2; X | Y1)
%
%   from which
%
%           I(Y2; X | Y1) = I(Y1 Y2; X) - I(Y1; X).
%
%   If we substitute Y1 = Ypast, Y2 = Xpast and X = Ypres we obtain
%
%           TE = I(Ypres ; Xpast | Ypast) =
%              = I(Ypast Xpast; Ypres) - I(Ypast; Ypres).
%
%   which expresses transfer entropy as the information about Ypres carried
%   by Xpast beyond that carried by Ypast. Here Ypast and Xpast can be
%   multi-dimensional arrays, e.g.
%
%           Ypres = Y(t);
%           Ypast = [Y(t-ty_1), ..., Y(t-ty_ny)];
%           Xpast = [X(t-tx_1), ..., X(t-tx_nx)];
%
%   where ty_1, ..., ty_ny are the ny delays for signal Y and tx_1, ...,
%   tx_nx are the nx delays for signal X.
% 
%   This formulation of the transfer entropy is also particularly
%   convenient for computation since Ypres is always one-dimensional while
%   Xpast and [Xpast Ypast] can be multi-dimensional responses.
%   Consequently:
%   - we can use the shuffling correction in the best possible way;
%   - bias increases due to increases in the dimensionality of Ypast (while
%     keeping the dimensionality of Xpast constant) will cancel each other
%     in the when taking the difference between the two information terms.
%
%   -----------------------------------------------------------------------
%   TECHNICAL NOTE 2 - HOW GOUREVITCH AND EGGERMONT NORMALIZATION IS 
%   COMPUTED IN THIS ROUTINE
%   -----------------------------------------------------------------------
%   The Normalized Transfer Entropy [2] is defined as
%
%           NTE = TE / H(Ypres | Ypast)
%
%   Using the chain rule for entropy we obtain:
%
%           H(X1 X2) = H(X1) + H(X1|X2)   -->   H(X2|X1) = H(X1 X2) - H(X1)
%
%   Sustituting X2 = Ypres and X1 = Ypast we have:
%
%           H(Ypres | Ypast) = H(Ypres Ypast) - H(Ypres)
%
%   REFERENCES
%   ----------
%   [1] TM Cover, JA Thomas, "Elements of Information Theory - Second
%       Edition", Wiley Interscience
%
%   [2] B Gourevitch, JJ Eggermont, "Evaluating Information Transfer between
%       Auditory Cortical Neurons", Journal of Neurophysiology, 97, 2533-2543
%
%   Copyright (C) 2010 Cesare Magri
%   Version: 1.1.0
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

if nargin < 4
    error("transferentropy.m: not enough input arguments. See `help transferentropy` for usage info");
end
if any(isnan(X))
    error("R contains NaNs. Aborting.")
end
if any(isnan(Y))
    error("S contains NaNs. Aborting.")
end
if ~isfield(opts, 'bin_method')
    warning("No opts.bin_method specified. R is assumed to be already binned. Proceeding without binning response matrix.")
    opts.bin_method = 'none';
    opts.n_bins = 1;
end
if ~isfield(opts, 'n_bins')
    error("Please specify number of bins for response binning through opts.n_bins variable")
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
if ~isfield(opts, 'bin_method')
    warning("No opts.bin_method specified. R is assumed to be already binned. Proceeding without binning response matrix.")
    opts.bin_method = 'none';
    opts.n_bins = 1;
end
if ~isfield(opts, 'n_bins')
    error("Please specify number of bins for response binning through opts.n_bins variable")
end
if ~isfield(opts, 'taux')
    error("Please specify a vector of (strictly negative integer) lags for variable X")
end
if ~isfield(opts, 'tauy')
    error("Please specify a vector of (strictly negative integer) lags for variable Y")
end
if ~isfield(opts, 'method')
    error("Please specify a probability estimation method")
end
if ~isfield(opts, 'bias')
    error("Please specify a bias correction strategy")
end
if ~isfield(opts, 'trperm')
    opts.trperm = 0;
end 

xtau = opts.taux;
ytau = opts.tauy;

if any(xtau>=0 | ytau>=0)
    error('Only negative delays can be specified.');
end

[nPntX, nTrlX] = size(X);
[nPntY, nTrlY] = size(Y);
if nPntX~=nPntY && nTrlX~=nTrlY
    error('Size of the two input matrices must be the same.');
else
    nTrl = nTrlX;
    nPnt = nPntX;
end

whereTE    = strcmpi(outputsList, 'te'   );
whereTEsh  = strcmpi(outputsList, 'tesh' );
whereNTE   = strcmpi(outputsList, 'nte'  );
whereNTEsh = strcmpi(outputsList, 'ntesh');

doTE    = true; % we always need to compute TE
doTEsh  = any(whereTEsh);
doNTE   = any(whereNTE);
doNTEsh = any(whereNTEsh);

nOut = any(whereTE) + doTEsh + doNTE + doNTEsh;

% Making btsp and trperm options mutually exclusive
if (isfield(opts, 'btsp') && opts.btsp>0) && (isfield(opts, 'trperm') && opts.trperm>0)
    error("Options `btsp` and `trperm` cannot be specified together.");
elseif (isfield(opts, 'btsp') && opts.btsp>0)
    nTestRep = opts.btsp;
elseif (isfield(opts, 'trperm') && opts.trperm>0)
    % Are there trials at all?
    if nTrl==1
        error("No trials to permute.");
    end
    nTestRep = opts.trperm;
else
    nTestRep = 0;
end

% Creating matrices with delayed copies of the signals 
outLen = nPnt + min([xtau(:); ytau(:)]); % Using "+" since all tau are negative!

x_with_delays = subTimeShiftedReplica(X, xtau(:), min([xtau(:); ytau(:)]), outLen);

% For Y we include also delay zero:
ytau_with_zero = [0; ytau(:)];
y_with_delays = subTimeShiftedReplica(Y, ytau_with_zero, min([xtau(:); ytau(:)]), outLen);

% Computing transfer entropy values 
TEsh  = zeros(nTestRep+1, 1);
NTE   = zeros(nTestRep+1, 1);
NTEsh = zeros(nTestRep+1, 1);

% Computing transfer-entropy quantities:
[TE, TEsh, NTE, NTEsh] = subTeEngine(x_with_delays(:,:), y_with_delays(:,:), opts, doTE, doTEsh, doNTE, doNTEsh);

% Random trial pairing
for testRepInd=1:opts.trperm
    % Randomly permuting trials in XMAT:
    x_with_delays = x_with_delays(:, :, randperm(nTrl));
    
    [TE(testRepInd+1), TEsh(testRepInd+1), NTE(testRepInd+1), NTEsh(testRepInd+1)] = ...
        subTeEngine(x_with_delays(:,:), y_with_delays(:,:), opts, doTE, doTEsh, doNTE, doNTEsh);
end

% Assigning output
transferentropies = cell(nOut,1);
transferentropies(whereTE)    = {TE};
transferentropies(whereTEsh)  = {TEsh};
transferentropies(whereNTE)   = {NTE};
transferentropies(whereNTEsh) = {NTEsh};


function [TE, TEsh, NTE, NTEsh] = subTeEngine(xmat, ymat, OptInfo, ...
                                  doTE, doTEsh, doNTE, doNTEsh)
%SUBTEENGINE Computes the transfer-entropy values
%
% USAGE
%   [TE, TESH, NTE, NTESH] = SUBTEENGINE(XMAT, YMAT, OPTINFO, DOTE, DOTESH, DONTE, DONTESH)
%
% INPUT
%   XMAT
%   YMAT
%   OPTINFO
%   DOTE
%   DOTESH
%   DONTE
%   DONTESH
%
% OUTPUT
%   

% I1 = I(Ypast; Ypres)
if doTE || doNTE
    I1 = information(ymat(2:end,:), ymat(1,:), OptInfo, {'I'});
    I1 = I1{1};
end
        
if doTEsh || doNTEsh
%     I1sh = information(R1, OptInfo, {'Ish'});
    I1sh = information(ymat(2:end,:), ymat(1,:), OptInfo, {'Ish'});
    I1sh = I1sh{1};
end


% I2 = I(Ypast Xpast; Ypres)
if doTE || doNTE
    I2 = information([ymat(2:end,:); xmat], ymat(1,:), OptInfo, {'I'});
    I2 = I2{1};
end

if doTEsh || doNTEsh
    I2sh = information([ymat(2:end,:); xmat], ymat(1,:), OptInfo, {'Ish'});
    I2sh = I2sh{1};
end

% Computing normalization factor:
if doNTE || doNTEsh
    % H1 = H(Ypres Ypast)
    %OptInfo.nt = size(ymat(:,:),2);
    OptInfo.bias = 'pt';
    HR1 = entropy(ymat(:,:), ymat(1,:), OptInfo, {'HR'});
    HR1 = HR1{1};
    
    % H2 = H(Ypres)
    %OptInfo.nt = size(ymat(:,:),2);
    HR2 = entropy(ymat(1,:), ones(size(ymat(1,:))), OptInfo, {'HR'});
    HR2 = HR2{1};
    
    normFact = HR1 - HR2;
end

% Computing transfer-entropy:
if doTE || doNTE
    TE = I2 - I1;
else
    TE = 0;
end

% Computing transfer entropy with shuffling correction:
if doTEsh || doNTEsh
    TEsh = I2sh - I1sh;
else
    TEsh = 0;
end

% Computing normalized transfer entropy:
if doNTE
    NTE = (I2 - I1) / normFact;
else
    NTE = 0;
end

% Computing normalized transfer entropy with shuffling correction:
if doNTEsh
    NTEsh = (I2sh - I1sh) / normFact;
else
    NTEsh = 0;
end


function out = subTimeShiftedReplica(x, negTauVec, minNegTau, outLen)
%SHIFTED_REPLICA Returns matrix with replica of an array shifted by
% different negative delays.
%
% USAGE
%   OUT = SHIFTED_REPLICA(X, NEGTAUVEC, OUTLEN)
%
% INPUT
%   X         - array that we wish to time-shift
%   NEGTAUVEC - array of negative time shifts
%   OUTLEN    - desired output length (must be smaller than LENGTH(X))

nTrl = size(x,2);
negTauVecLen = length(negTauVec);
% minNegTau = min(negTauVec);

out = zeros(negTauVecLen, outLen, nTrl);

for negTauInd=1:negTauVecLen

    negTau = negTauVec(negTauInd);

    posTau = negTau - minNegTau;

    if negTau<0
        xshift = subIntShift(x, negTau);
    else
        xshift = x;
    end

    if posTau>0
%         out(negTauInd, :, :) = subIntShift(xshift, posTau);
        tmp = subIntShift(xshift, posTau);
        out(negTauInd, :, :) = tmp(1:outLen, :);   
    else
        out(negTauInd, :, :) = xshift(1:outLen, :);
    end

end


function y = subIntShift(x, tau)
%INTSHIFT4TRANSFERENTROPY Shifts a signal in time.
%
% USAGE
%   Y = TIMESHIFT(X, TAU);
%   Y = TIMESHIFT(X, TAU, N);
%
% INPUT
%   X   - time series to be time-shifted. If X is an N-by-M array than each
%         row of X corresponds to a sample of the time-series while the M
%         columns are M distinct observations.
%   TAU - lag by which the signal needs to be shifted (can be negative or
%         positive).

if tau>0
    y = x(tau+1:end, :);
    
elseif tau==0
    y = x;
    
elseif tau<0
    y = x(1:end+tau, :); % Using "+" since tau negative
    
end
