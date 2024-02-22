% Given 3 variables S, R and C; this function computes the 
% intersection information (II) from S to C trhough R defined as 
%
% II = min({S}{R}_C, {C}{R}_S)
%
% where the '_target' notation indicates the target of the information
% decomposition for the considered PID atom


function [II] = calculate_intersection_information(S,C,R,n_trials,n_split_trials)
        % random sample n_split_trials from data in each of the
        % n_draws and calculate II on the reduced dataset
        ri = randperm(n_trials, n_split_trials);
        Sr = S(ri);
        Cr = C(ri);
        Rr = R(ri);
        [prob_matrix, ~, nS, nR, nC, n_s_d] = build_p(Sr, Rr, Cr);
        python_ComputeII = py.ComputePID.ComputeII(prob_matrix,...
             nS, nR, nC,...
             n_s_d);
        II = double(python_ComputeII.calculate());
end