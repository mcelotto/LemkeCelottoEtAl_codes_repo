% function info=NPC_information(X,opts)

clear all
close all


addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Selmaan_Switch/glmnet'));
addpath(genpath('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/TOOLS/ENTROPY_METHODS'));
path_list_cell = regexp(path,pathsep,'Split');
if ismember('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno',path_list_cell)
    rmpath(('/home/hs258/Codes_Folder/HOUMAN_GITLAB_IIT/Copula_Arno'))
end
addpath(('/home/hs258/Codes_Folder/Houman_Git/NPC_Info/NPC_Info'));



%%%%%%% load sample data
load('NPC_X.mat')


%%%%%% parameters
opts.bw='LL1';                          %% LL1 or LL2 bandwidth methods
opts.knots_fit=100;                     %% number of bins used in estimating the copula 
opts.knots_est=opts.knots_fit;
opts.type={'cont','cont'};              %% continous or discrete marginals
% opts.type={'discrete','discrete'};
opts.parallel=0;                        %% 0=non paralle, 1=parallel computing 
opts.alpha=0.05;                        %% alpha for monte-carlo sampling varinace
opts.erreps=1e-3;                       %% variance threshhold for monte-carlo sampling for information estimation
opts.iter=50;                           %% maximum number of monte-carlo iterations 
opts.cases=20000;                       %% number of samples in each iteration 
opts.plot=0;                            %% 0=no iteration plot (default), 1=with iteration plot



%%%% computing I(X(:,1);X(:,2))
range(1:2,1)=min(X(:,1:2))-1e-10;
range(1:2,2)=max(X(:,1:2))+1e-10;

[vine]=NPC_prep_copula(X(:,[1 2]),opts.type,range([1 2],:));

%%%%% fitting the bandwidths of copula
[density_X,~,copula,~,~] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,1,0,opts.knots_fit,opts.parallel);
%%%%% evaluating the copula
[~,~,copula,~,~] = NPC_Fit_vCopula(vine,X(1,:),opts.bw,-1,copula,opts.knots_est,opts.parallel);
%%%%% computing the information
[info,~,~,~] = NPC_kernelvineinfo(vine,copula,opts)





