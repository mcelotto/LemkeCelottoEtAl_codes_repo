function [MI] = ...
    svm_MI(MI_opts, MI_calc_mode, data, target_labels, test_set, tr_method, varargin)
%%% *function [MI] = svm_MI(MI_opts, MI_calc_mode, data, target_labels, test_idx, tr_method, varargin)*
%%%
%%% ### Description
%%% This function uses a SVM to calculate the mutual information using: a) the predicted stimuli labels for each trial by the SVM, b) the posterior probabilties obtained by the SVM.
%%%
%%% ### Inputs:
%%% - *MI_opts*: options to be passed to information.m
%%% - *MI_calc_mode*: operating mode for SVM calculation of MI. Should be either:
%%%           - *'Sp_S'* : calculates mutual information between stimulus and predicted stimulus by SVM.
%%%           - *'post_prob'* : uses stimulus posterior probabilities as given by the SVM.
%%% - *data*: *n_samples x n_features* predictor array.
%%% - *target_labels*: original applied stimulus for each trial.
%%% - *test_set*: fraction of trials that will be used for testing. Must be between 0 and 1. It can also be a vector with indices to select test trials from data.
%%% - *tr_method*: traning method for the SVM. Should be either:
%%%           - *'KFold'* : K-fold cross validation. This option should be followed by the number of folds as an additional argument.
%%%           - *'LeaveOneOut'* : equivalent to a K-fold cross validation where K is the number of samples in the database.
%%% - varargin{1} = *k* Folds if tr_methods is 'KFold'.
%%%
%%% ### Outputs:
%%% - MI: Mutual Information


% check inputs
assert(size(data,1) == size(target_labels,1),...
    "Number of labels should equal the number of entries in the dataset");
assert(strcmp(MI_calc_mode,'Sp_S') || strcmp(MI_calc_mode,'Post_prob'), ...
    "Calculation mode for MI should be either 'Sp_S' or 'Post_prob'. Specified MI_calc_mode option: '%s'", MI_calc_mode);
assert(strcmp(tr_method,'KFold') || strcmp(tr_method,'LeaveOneOut'),... 
    "Training method for MI should be either 'KFold' or 'LeaveOneOut'. Specified MI_calc_mode option: '%s'", tr_method);
assert(strcmp(tr_method,'KFold')  && nargin > 5,...
    "Training method 'KFold' requires the specification of the number of folds as additional argument");

% build SVM classifier on data
[predicted_labels, test_labels, posterior_probs, test_data] = ...
    buildML(data, target_labels, test_set, 'SVM', tr_method, varargin{1});

if strcmp(MI_calc_mode,'Sp_S')
    % calculate MI between test labels and predicted labels
    [R, nt] = buildx(test_labels, predicted_labels');
    MI_opts.nt = nt;
    MI = information(R, MI_opts, 'I');
else
    % calculate MI using posterior probabilities of the SVM
    n_data = length(test_data(:,1));
    X = unique(test_data,'rows');
    marg_prob_X = zeros(length(X),1);
    joint_prob = zeros(length(marg_prob_X), length(unique(target_labels)));
    % calculate p(X)
    for i = 1:length(X)
        l = ismember(test_data,X(i),'rows');
        marg_prob_X(i) = sum(l);
        posterior_probs_ave(i,:) = mean(posterior_probs(logical(l),:),1);
    end
    marg_prob_X = marg_prob_X/n_data;
    
    % calculate p(Y;X)
    for i = 1:length(joint_prob(:,1))
        for j = 1:length(joint_prob(1,:))
            joint_prob(i,j) = posterior_probs_ave(i,j)*marg_prob_X(i);
        end
    end
    MI = joint_prob_MI(joint_prob);
end
end

