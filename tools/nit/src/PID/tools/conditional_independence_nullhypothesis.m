function [measure_null,pval] = conditional_independence_nullhypothesis(data, value, opts)
%%
% This function compute the 'conditional independence given the stimulus'
% null hypothesis for the intersection information (II) and the 
% feature-specific information transmission (FIT) measures.

% Inputs:
% data should be a N x nTrials array where:
% For II --> N = 3. row 1 = stimulus; row 2 = neural respose; row 3 = choices
% For FIT --> N = 4. row 1 = stimulus; row 2 = emitter paste; row 3 = receiver past; row 4 = receiver present
% value: information value on the original distribution

if size(data,1) == 3
    infotype = 'II';
    S = data(1,:)';
    X = data(2,:);
    Y = data(3,:)';
    Z = [];
elseif size(data,1) == 4
    infotype = 'FIT';
    S = data(1,:)';
    X = data(2,:);
    Y = data(3,:);
    Z = data(4,:);
end

n_permutations = opts.nh_perm;

nTrials = size(S,1);

measure_null = zeros(1, n_permutations);

Sval = unique(S);

for sidx = 1:n_permutations
    XSh = zeros(1,numel(X));

    for Ss = 1:numel(Sval)

        idx = (S == Sval(Ss));  % select trials where stim = Ss

        hX = X(idx);  % take X values on those trials

        ridx = randperm(sum(idx));  % generate random idxs

        XSh(1, idx) = hX(ridx);  % assign to XSh, on idxs where stim = Ss, reshuffled values of X at stim = Ss

    end
    
    if strcmp(infotype, 'II')
        % n_split_trials = nTrials because we compute the null hypothesis for the finite-sampling biased measure
        measure_null(sidx) = calculate_intersection_information(S, Y, XSh, nTrials, nTrials);
    elseif strcmp(infotype, 'FIT')
        measure_null(sidx) = calculate_feature_specific_information_transmission(S, XSh, Y, Z, nTrials, nTrials);
    end
end

pval = mean(measure_null > value);

end