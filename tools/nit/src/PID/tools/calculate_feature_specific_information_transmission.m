% Given 4 variables S, hX, hY and Y; this function computes the 
% feature-specific information transmission (FIT) defined as 
%
% FIT = min({hX}{Y}_S, {hX}{S}_Y)
%
% where the '_target' notation indicates the target of the information
% decomposition for the considered PID atom

function [FIT] = calculate_feature_specific_information_transmission(S,hX,hY,Y,n_trials,n_split_trials)
        % random sample n_split_trials from data in each of the
        % n_draws and calculate II on the reduced dataset
        ri = randperm(n_trials, n_split_trials);
        Sr = S(ri);
        hXr = hX(ri);
        hYr = hY(ri);
        Yr = Y(ri);
        [p_hXYhYS] = probabilityDist(hXr, Yr, hYr, Sr');
        FIT_S = computePiXY(p_hXYhYS);
        p_hXShYY = permute(p_hXYhYS, [1 4 3 2]);
        FIT_Y = computePiXY(p_hXShYY);
        FIT = min([FIT_S, FIT_Y]);
end