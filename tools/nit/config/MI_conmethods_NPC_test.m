%Load the sample bivariate data:
load('NPC_X.mat')

%The parameters which can be selected in the copula fitting and information estimation are as follows:
opts.method = 'NPC';
opts.bw='LL1';                               %% LL1 or LL2 bandwidth methods
opts.knots_fit=30;                           %% number of bins used in estimating the copula 
opts.knots_est=opts.knots_fit;
opts.margins={'cont','cont'};                %% continous 'cont' or discrete 'discrete' marginals
opts.parallel=0;                             %% 0=non paralle, 1=parallel computing 
opts.alpha=0.05;                             %% alpha for monte-carlo sampling variance
opts.erreps=1e-3;                            %% variance threshhold for monte-carlo sampling for information estimation
opts.iter=2;                                 %% maximum number of monte-carlo iterations 
opts.cases=5000;                             %% number of samples in each iteration 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Computing I(X(:,1);X(:,2)):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = X';
[MI] = copulaInformation(X(1,:), X(2,:), opts);