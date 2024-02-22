% Given 4 variables S, hX, hY and Y; this function computes the 
% feature-specific information transmission (FIT) defined as 
%
% FIT = min({hX}{Y}_S, {hX}{S}_Y)
%
% where the '_target' notation indicates the target of the information
% decomposition for the considered PID atom

function [FIT] = calculate_choice_specific_information_transmission(C,hX,hY,Y,n_trials,n_split_trials)
        % random sample n_split_trials from data in each of the
        % n_draws and calculate II on the reduced dataset
        ri = randperm(n_trials, n_split_trials);
        Cr = C(ri);
        hXr = hX(ri);
        hYr = hY(ri);
        Yr = Y(ri);
        [p_hXYhYC] = probabilityDist(hXr, Yr, hYr, Cr');
        FIT_C = computePiXY(p_hXYhYC);
        p_CYhYhX = permute(p_hXYhYC, [4 2 3 1]);
        FIT_hX = computePiXY(p_CYhYhX);
        FIT = min([FIT_C, FIT_hX]);
end