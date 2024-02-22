function varargout = buildML(data,labels,testSet,algorithm,cv_type,varargin)
%%% *function varargout = buildML(data, labels, testSet, algorithm, cv_type, varargin)*
%%%
%%% ### Description
%%% buildML trains a Machine Learning model with the chosen cross-validation algorithm on a specific set of data points and corresponding labels.
%%%
%%% ### Inputs:
%%% - *data*: n_samples x n_features predictor array.
%%% - *labels*: n_samples x 1 array of labels to be learnt.
%%% - *testSet*: fraction of trials that will be used for testing. Must be between 0 and 1. It can also be a vector with indices to select test trials from data.
%%% - *algorithm*: string specifying the ML algorithm. Must be either *'SVM'* or *'GLM'*.
%%% - *cv_type*: string speciying the type of crossvalidation to be performed. Must be either *'KFold'* or *'LeaveOneOut'*.
%%% - *K* (only if using 'KFold'): number of sets for 'KFold' crossvalidation.
%%%
%%% ### Outputs:
%%% - *predictedLabels*: predicted labels for samples in the held out set.
%%% - *labelsTest*: true labels for samples in the held out set.
%%% - *dataTest*: data used for the held out set.
%%% - *posteriorProbs* (only with SVM): probabilities of each label given the trial.

if any(isnan(data))
    error("data contains NaNs. Aborting.")
end
if any(isnan(labels))
    error("labels contains NaNs. Aborting.")
end
if any(isnan(testSet))
    error("testSet contains NaNs. Aborting.")
end

assert(strcmp(cv_type,'KFold') || strcmp(cv_type,'LeaveOneOut'),...     %check if the CV method agrees with the expected ones
    'The chosen crossvalidation method is not acceptable')
assert(strcmp(algorithm,'SVM') || strcmp(algorithm,'GLM'),...           %check if the algorithm agrees with the expected ones
    'The chosen algorithm method is not acceptable')
assert(~(nargin < 6 && strcmp(cv_type,'KFold')),...                     %check if all the arguments are given together with KFold method
    ['With KFold crossvalidation you need to provide',...
    ' the amount of subsets (k)'])
if (nargin > 6 && strcmp(cv_type,'KFold')) || ...                       %warn the user if he/she is giving more arguments than the ones required
        (nargin > 5 && strcmp(cv_type,'LeaveOneOut'))
    warning('Extra arguments will be ignored')
end

SVM_flag = strcmp(algorithm,'SVM');                                     %which algorithm are we using
GLM_flag = strcmp(algorithm,'GLM');

data = data(~any(isnan(data), 2),:);                                    %cleaning data
labels=labels(~any(isnan(data), 2));                                    %cleaning data

nTrials = length(labels);                                               %total trials
if length(testSet) > 1
    nTrialsTest = length(testSet);
else
    nTrialsTest = floor(testSet*nTrials);                               %rials for test
end
nTrialsTrain = nTrials-nTrialsTest;                                     %trials for train

if nTrialsTrain==0 || nTrialsTest==0                                    %check if data is correctly splitted for train and test
    error(['The number of trials for test or train cannot be zero.'...
        ' Modify your test fraction'])
end

if length(testSet) > 1
    train_test_partition.test = testSet;
    train_test_partition.training = setdiff(1:size(data,1),testSet);
else
    train_test_partition = cvpartition(labels,'HoldOut',testSet);       %dividing data between train and test with a given fraction
end

dataTrain = data(train_test_partition.training',:);                     %labels for train
labelsTrain = labels(train_test_partition.training');                   %trials for train
dataTest = data(train_test_partition.test',:);                          %labels for test
labelsTest = labels(train_test_partition.test');                        %trials for test

if strcmp(cv_type,'KFold')
    K = varargin{1};                                                    %crossvalidation folds
    cvPartition = cvpartition(labelsTrain,'KFold',K);                   %split data for kFold crossvalidation
elseif strcmp(cv_type,'LeaveOneOut')
    cvPartition = cvpartition(labelsTrain,'LeaveOut');                  %split data for LeaveOneOut crossvalidation
end

setsCV = cvPartition.NumTestSets;                                       %crossvalidation sets

if SVM_flag                                                             %if SVM is chosen  
    for i = 1:setsCV                                                    %for each set of CV
        testIdx = find(cvPartition.test(i));                            %trials for test CV
        [parameters,accuracy(:,:,i)] = ...
            svm_training(dataTrain,labelsTrain,testIdx);                %get the c and gamma
    end
    
    mean_accuracy = mean(accuracy,3);                                   %get the average accuracy over CV sets
    [idx(:,1),idx(:,2)] = ...
        find(mean_accuracy == max(mean_accuracy(:)));                   %get the one that best generalizes
    idx = idx(1,:);
    c_opt = 2^parameters.c_vec(idx(1));                                 %get the value of the corresponding optimal c
    gamma_opt = 2^parameters.gamma_vec(idx(2));                         %get the value of the corresponding optimal gamma
    algorithm_opt = svmtrain(labelsTrain, dataTrain,...
        sprintf('-q -c %f -g %f -t 2 -b 1', c_opt, gamma_opt));         %training SVM with all trainingData and optimal params
    [predictedLabels, ~, posteriorProbs] = ...
        svmpredict(labelsTest, dataTest, algorithm_opt, '-q -b 1');     %predict using the testData
    varargout{1} = predictedLabels;                                     %assigning the outputs
    varargout{2} = labelsTest;
    varargout{3} = dataTest;
    varargout{4} = posteriorProbs;

elseif GLM_flag                                                         %if GLM chosen 
    foldId = zeros(size(labelsTrain));
    
    for i = 1:setsCV                                                    %for each set of CV
        foldId = foldId + i*cvPartition.test(i);
    end
    
    CVfit = cvglmnet(dataTrain,labelsTrain,...
        'multinomial',[],'mse',[],foldId);
    predictedLabels = cvglmnetPredict(CVfit,dataTest,...
        'lambda_min','class');
    varargout{1} = predictedLabels;                                     %assigning the outputs
    varargout{2} = labelsTest;
    varargout{3} = dataTest;
end
end

