function [parameters, accuracy]=...
    svm_training(data, target_labels, test_idxs)
%%% *function [parameters, accuracy] = svm_training(data, target_labels, test_idxs)*
%%% 
%%% ### Description
%%% Fits an rbf kernel svm to a series of data points and return the optimal fit parameters and accuracy on a held out set.
%%%
%%% ### Inputs: 
%%% - *data*: *n_samples x n_features* predictor array.
%%% - *target_labels*: *n_samples x 1* array of labels to be learnt.
%%% - *test_idxs*: positional indices of the samples to be held
%%% out of the SVM fitting procedure. 
%%%
%%% ### Outputs:
%%% - *parameters*: struct containing the values of c and gamma y separate
%%% vectors (c_vec and gamma_vec) respectively.
%%% - *accuracy*: matrix of *length(d_c) x length(d_gamma)* containing the accuracy associated to every possible pair of parameters.

    n_trials = size(data, 1);                                               %number of trials after cleaning
   
    train_idxs = setdiff(1:n_trials, test_idxs);                            %indices associated to the training trials

    c_exp = linspace(-5, 10, 16);                                           %c values to test
    gamma_exp = linspace(-15, 10, 26);                                      %gamma values to test
    
    d_c = numel(c_exp);                                                     %amount of values that c can take
    d_gamma = numel(gamma_exp);                                             %amount of values that gamma can take
    accuracy=zeros(d_c, d_gamma);                                           %initializing the array of accuracies corresponding to each pair of c and gamme

    trainData = data(train_idxs,:);                                         %assigning the training data
    trainLabel = double(target_labels(train_idxs,:));                       %assigning the training labels (and cast to double to keep libsvm happy)
    testData = data(test_idxs,:);                                           %assigning the test data
    testLabel = double(target_labels(test_idxs,:));                         %assigning the test labels (and cast to double to keep libsvm happy)

    testData=sparse(testData);                                              %converting trainData to sparse
    trainData=sparse(trainData);                                            %converting testData to sparse

    parameters.c_vec = c_exp;                                                 %saving the values of c in a struct for output
    parameters.gamma_vec = gamma_exp;                                         %saving the values of gamma in a struct for output
    
    for c_idx=1:d_c                                                         %for each value of c
        for gamma_idx=1:d_gamma                                             %for each value of gamma

            % Libsvm options:
            % -s svm_type : set type of SVM (default 0)
            % 	0 -- C-SVC
            % 	1 -- nu-SVC
            % 	2 -- one-class SVM
            % 	3 -- epsilon-SVR
            % 	4 -- nu-SVR
            % -t kernel_type : set type of kernel function (default 2)
            % 	0 -- linear: u'*v
            % 	1 -- polynomial: (gamma*u'*v + coef0)^degree
            % 	2 -- radial basis function: exp(-gamma*|u-v|^2)
            % 	3 -- sigmoid: tanh(gamma*u'*v + coef0)
            % -d degree : set degree in kernel function (default 3)
            % -g gamma : set gamma in kernel function (default 1/num_features)
            % -r coef0 : set coef0 in kernel function (default 0)
            % -c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)
            % -n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)
            % -p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)
            % -m cachesize : set cache memory size in MB (default 100)
            % -e epsilon : set tolerance of termination criterion (default 0.001)
            % -h shrinking: whether to use the shrinking heuristics, 0 or 1 (default 1)
            % -b probability_estimates: whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)
            % -wi weight: set the parameter C of class i to weight*C, for C-SVC (default 1)
            % -q : quiet mode (no outputs)
            
            c_aux=2^c_exp(c_idx);                                           %value of c to pass to svmtrain
            gamma_aux=2^gamma_exp(gamma_idx);                               %value of gamma to pass to svmtrain

            model = svmtrain(trainLabel, trainData, ...
                sprintf('-q -c %f -g %f -t 2 -b 1', c_aux, gamma_aux));     %training the model

            [~, accuracy_temp, ~] = ...
                svmpredict(testLabel, testData, model, '-q -b 1');          %run the model on the test data
            
            accuracy(c_idx,gamma_idx)=accuracy_temp(1);                     %saving the accuracy corresponding to each pair of c and gamma for output
        end
    end
end
