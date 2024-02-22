function [MI,varargout] = copulaInformation(X, Y, opts)
%%% *function [MI,copula] = copulaInformation(X, Y, opts)*
%%%
%%% ### Description
%%% Compute mutual information using mixed-vine or non-parametric copulae.
%%%
%%% ### Inputs:
%%% - *X*: *nDimensionsX X nTrials* array including input data (typically the neural response).
%%% - *Y*: *1 X nTrials* array including input data (typically the stimulus). 
%%% - *opts*: options sctructure:
%%%       - *opts.method*: *'MVC'* for using mixed-vine copula and *'NPC'*
%%%       for non-parametric copula (see tables below for further fields depending on the method).
%%%
%%%         ## MVC
%%%
%%%     |          Field         | Optional |                                                                                                                                                              Description                                                                                                                                                             | Default Value |
%%%     |:----------------------:|:--------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|:-------------:|
%%%     |  opts.margins{i}.dist  |   False  | *'unif'* (uniform), *'norm'* (normal), *'gam'* (gamma), *'poiss'* (poisson),<br>*'bino'* (binomial), *'nbin'* (negative binomial).<br>*'norm'* and *'gam'* require opts.marginal{i}.iscont = True while *'nbin'* <br>*'bino'* *'nbin'* require opts.marginal{i}.iscont = False.<br>Note that *i* makes reference to the ith marginal |       —       |
%%%     |  opts.margins{i}.theta |   False  | *nParams x 1* array including the parameters of the respective marginal                                                                                                                                                                                                                                                              |       -       |
%%%     | opts.margins{i}.iscont |   False  | Binary variable specifying whether the marginal is continuous (*True*)<br>or discrete (*False*)                                                                                                                                                                                                                                      |       -       |
%%%     |   opts.families{i,j}   |   False  | Dependency between variables in marginals *i* and *j*. *nMarginals x nMarginals*<br>array with non-empty values in the upper triangle. <br>The values must match any of: *'gaussian' 'student' 'clayton' 'claytonrot090'<br> 'claytonrot180' 'claytonrot270'*                                                                        |       -       |
%%%     |     opts.theta{i,j}    |   False  | Parameters for each family in opts.families relating variables *i* and *j*                                                                                                                                                                                                                                                           |       0       |
%%%     |       opts.alpha       |   True   | Significance level of estimate                                                                                                                                                                                                                                                                                                       |      0.05     |
%%%     |       opts.erreps      |   True   | Maximum standard error                                                                                                                                                                                                                                                                                                               |      1e-3     |
%%%     |       opts.cases       |   True   | number of draw samples                                                                                                                                                                                                                                                                                                               |     20000     |
%%%     |      opts.verbose      |   True   | display extended output info                                                                                                                                                                                                                                                                                                         |     False     |
%%%
%%%         ## NPC
%%%
%%%     |      Field     | Optional |                                 Description                                 | Default Value |
%%%     |:--------------:|:--------:|:---------------------------------------------------------------------------:|:-------------:|
%%%     |  opts.margins  |   False  | *'cont'* and *'discrete'* for continuous or discrete marginals respectively |       —       |
%%%     |     opts.bw    |   True   |                     *'LL1'* or *'LL2'* bandwidth methods                    |     'LL1'     |
%%%     | opts.knots_fit |   True   |                  number of bins used in fitting the copula                  |      100      |
%%%     | opts.knots_est |   True   |                 number of bins used in estimating the copula                |      100      |
%%%     |  opts.parallel |   True   |               *0* for non parallel and *1* parallel computing               |       0       |
%%%     |   opts.alpha   |   True   |                   alpha for monte-carlo sampling variance                   |      0.05     |
%%%     |   opts.erreps  |   True   |    variance threshold for monte-carlo sampling for information estimation   |      1e-3     |
%%%     |    opts.iter   |   True   |                   maximum number of monte-carlo iterations                  |       50      |
%%%     |   opts.cases   |   True   |                     number of samples in each iteration                     |     20000     |
%%%     | opts.verbose   |   True   | display extended output info                                                |     False     |
%%%
%%% ### Outputs:
%%% - *MI*: mutual information estimated from the copula
%%% - *copula*: copula fitted to X and Y. The data structure will be different depending on the used copula structure (NPC Vs MVC).


%input checking
assert(nargin>=3,'Not enough input arguments')
assert(nargout<3, 'Too many output arguments')
assert(isfield(opts, 'method'),'Please provide which copula method do you want to use "MVC" or "NPC"')

assert(size(Y,2)==size(X,2),'Amount of trials is not consistent between X and Y')

if any(isnan(X))
    error("X contains NaNs. Aborting.")
end
if any(isnan(Y))
    error("Y contains NaNs. Aborting.")
end
if ~isfield(opts, 'margins')
        error('Missing mandatory argument: please specify type of marginal distributions.')
else
    for i = 1:length(opts.margins)
       if ~contains(opts.margins{i}, ["cont" "discrete"])
           error("Specified marginals do not match any of the possible options. Use 'cont' for continuous marginals and 'discrete' for discrete ones.")
       end
    end
end
if ~isfield(opts, 'verbose')
    opts.verbose = false;
end

NPCflag = 0; MVCflag = 0;
if strcmp(opts.method,'MVC')
    
    parameters.unif = [0 0];
    parameters.norm = [2 1];
    parameters.gam = [2 1];
    parameters.poiss = [1 1];
    parameters.bino = [2 1];
    parameters.nbin = [2 1];

    parameters.ind = [0 0];
    parameters.gaussian = [1 1];
    parameters.student = [2 1];
    parameters.clayton = [1 1];
    parameters.claytonrot090 = [1 1];
    parameters.claytonrot180 = [1 1];
    parameters.claytonrot270 = [1 1];
    
    MVCflag = 1;
    warn = 0;
    if ~isfield(opts, 'vinetype') || ...
            (isfield(opts, 'vinetype') && ~strcmp(opts.vinetype,'c-vine'))
        error('Please define the vine type in opts according to the documentation')
    end
    if ~isfield(opts, 'margins') || ...
            (isfield(opts, 'margins') && length(opts.margins)~=(size(X,1)+1))
        error('The amount of margins defined does not match the amount of dimensions')
    else
        for i = 1:length(opts.margins)
            
            if ~isfield(opts.margins{i}, 'dist') || (isfield(opts.margins{i}, 'dist') ...
                    && ~strcmp(opts.margins{i}.dist,'ind') && ~strcmp(opts.margins{i}.dist,'unif') ...
                    && ~strcmp(opts.margins{i}.dist,'norm') && ~strcmp(opts.margins{i}.dist,'gam') ...
                    && ~strcmp(opts.margins{i}.dist,'poiss') && ~strcmp(opts.margins{i}.dist,'bino') ...
                    && ~strcmp(opts.margins{i}.dist,'nbin'))
                error(['Marginal ',num2str(i),' type is not well defined in "dist" field'])
            end
            if ~isfield(opts.margins{i}, 'theta')
                error(['Marginal ',num2str(i),' parameters are not defined in "theta" field'])
            else
                if size(opts.margins{i}.theta) ~= parameters.(opts.margins{i}.dist)
                    error(['The amount of parameters in "theta" for marginal ',num2str(i),' is not correct'])
                end
            end
            if ~isfield(opts.margins{i}, 'iscont')
                error(['Marginal ',num2str(i),' parameters are not defined in "iscont" field'])
            else
                if opts.margins{i}.iscont~=1 && opts.margins{i}.iscont~=0
                    error(['Field "iscont" conains an incorrect value for marginal ',num2str(i)])
                end
            end
        end
    end
    
    if ~isfield(opts, 'families')
        warning('Families are not defined in "families" field');
    else
        [idx(:,1), idx(:,2)] = ...
            ind2sub(size(opts.families),find(~cellfun('isempty', opts.families)));
        check = (idx == nchoosek(1:length(opts.margins),2)) == 0;
        if any(check(:))
            error('Your families are not defined in the upper triangle of the cell array in field "families"')
        end
        for i = 1:size(idx,1)
            if ~strcmp(opts.families{idx(i,1),idx(i,2)},'ind') ...
                    && ~strcmp(opts.families{idx(i,1),idx(i,2)},'gaussian') ...
                    && ~strcmp(opts.families{idx(i,1),idx(i,2)},'student') ...
                    && ~strcmp(opts.families{idx(i,1),idx(i,2)},'clayton') ...
                    && ~strcmp(opts.families{idx(i,1),idx(i,2)},'claytonrot090') ...
                    && ~strcmp(opts.families{idx(i,1),idx(i,2)},'claytonrot180') ...
                    && ~strcmp(opts.families{idx(i,1),idx(i,2)},'claytonrot270')
                error(['Family ',num2str(i),' type is not well defined in "dist" field'])
            end
        end
        [idx(:,1), idx(:,2)] = ...
            ind2sub(size(opts.theta),find(~cellfun('isempty', opts.theta)));
        check = (idx == nchoosek(1:length(opts.margins),2)) == 0;
        if any(check(:))
            error('Your parameters are not defined in the upper triangle of the cell array in field "theta"')
        end
        for i = 1:size(idx,1)
            if size(opts.theta{idx(i,1),idx(i,2)}) ~= parameters.(opts.families{idx(i,1),idx(i,2)})
                error(['The amount of parameters for family ',num2str(i),' is not correct'])
            end
        end
    end
    
    if ~isfield(opts, 'cases')
        warn = 1;
        opts.cases = 20000;
    end
    if ~isfield(opts, 'alpha')
        warn = 1;
        opts.alpha = 0.05;
    end
    if ~isfield(opts, 'erreps')
        warn = 1;
        opts.erreps = 1e-3;
    end
    
    if warn
        warning('Not all the parameters were provided in the options structure. Using default values instead.');
    end
    
    vine.type = opts.vinetype;
    vine.margins = opts.margins;
    vine.families = opts.families;
    vine.theta = opts.theta;
    cases = opts.cases;
    
elseif strcmp(opts.method,'NPC')
    NPCflag = 1;
    warn = 0;
    if size(X,1)>1
        error('Sorry, but non-parametric copula can only be used with one-dimensional data')
    end
    if ~isfield(opts, 'margins')
        error('Please provide information about the type of marginals in opts.margins')
    end
    if ~isfield(opts, 'bw')
        warn = 1;
        opts.bw = 'LL1';
    end
    if ~isfield(opts, 'knots_fit')
        warn = 1;
        opts.knots_fit = 30;
    end
    if ~isfield(opts, 'knots_est')
        opts.knots_est = opts.knots_fit;
    end
    if ~isfield(opts, 'parallel')
        warn = 1;
        opts.parallel=0;
    end
    if ~isfield(opts, 'alpha')
        warn = 1;
        opts.alpha=0.05;
    end
    if ~isfield(opts, 'erreps')
        warn = 1;
        opts.erreps=1e-3;
    end
    if ~isfield(opts, 'iter')
        warn = 1;
        opts.iter=2;
    end
    if ~isfield(opts, 'cases')
        warn = 1;
        opts.cases=5000;
    end
    opts.plot = 0;
    if warn
        warning('Not all the parameters were provided in the options structure. Using default values instead.')
    end
end

%computing
if MVCflag
    
    if opts.verbose
    	disp('Sampling from mixed copula vine...');
    end
    x = mixedvinernd(vine,cases);
    if opts.verbose
    	disp('Fitting parameters to samples...');
    end
    % Construct argument to specify which margins are continuous
    iscont = false(length(vine.margins),1);
    for i = 1:length(vine.margins)
        iscont(i) = vine.margins{i}.iscont;
    end
    vineest = mixedvinefit(x,vine.type,iscont);
    %%Estimate entropy
    if opts.verbose
    	disp('Estimating entropy of mixed copula vine...');
    end
    alpha = 0.05; % Significance level of estimate (95% confidence)
    erreps = 1e-1; % Maximum standard error
    [h,stderr] = mixedvineentropy(vine,alpha,erreps);
    % Display results
    disp([' Estimated entropy:          ' num2str(h) ' bit']);
    disp([' Standard error of estimate: ' num2str(stderr) ' bit']);
    fprintf('\n');
    
    MI = h;
    
    if nargout > 1
        varargout{1} = vineest;
    end
    
elseif NPCflag
    
    % merge X and Y in a single X to apply transformation in case of a discrete marginal (marginals are given in the same order as inputs)
    X = [X',Y']; 
    
    % convert discrete values to continous if necessary
    for i = 1:length(opts.margins)
        if strcmp(string(opts.margins{i}), "discrete")
            X(:,i) = NPC_discrete_to_cont(X(:,i));
            opt.margins{i} = 'cont';
        end
    end
    
    % swap response and stimuli to match NPC inputs requirement
    X = [X(:,2) X(:,1)];
    
    if opts.verbose
    	disp('Building the vine structure...')
    end
    
    % Build the vine structure
    range(1:2,1)=min(X)-1e-10;
    range(1:2,2)=max(X)+1e-10;
    [vine]=NPC_prep_copula(X,opts.margins,range);
    
    if opts.verbose
    	disp('Fitting the copula bandwidths (this may take some time)...')
    end
    % Fit bandwidths of copula to data:
    [~,~,copula,~,~] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,1,0,opts.knots_fit,opts.parallel);

    if opts.verbose
    	disp('Estimating copula over grid and data points (this may take some time)...')
    end
    % Estimate copula:
    [~,~,copula,~,~] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,-1,copula,opts.knots_est,opts.parallel);

    if opts.verbose
    	disp('Estimating Mutual Information (this may take some time)...')
    end
    % Calculate MI:
    [MI,~,~,~] = NPC_kernelvineinfo(vine,copula,opts);
   
    if nargout > 1
        varargout{1} = copula;
    end
end

end

