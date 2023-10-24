% Given 3 variables S, R and C; this function computes the 
% partial information decomposition of R and C about S


function PID = calculate_pid(S,C,R,n_trials,n_split_trials)
        % random sample n_split_trials from data in each of the
        % n_draws and calculate II on the reduced dataset
        ri = randperm(n_trials, n_split_trials);
        Sr = S(ri);
        Cr = C(ri);
        Rr = R(ri);
        [prob_matrix, ~, nS, nR, nC, n_s_d] = build_p(Sr', Rr', Cr');
        python_ComputePID = py.ComputePID.ComputePID(prob_matrix,...
             nS, nR, nC,...
             n_s_d);
        PID = cell2mat(cell(python_ComputePID.calculate()));
end