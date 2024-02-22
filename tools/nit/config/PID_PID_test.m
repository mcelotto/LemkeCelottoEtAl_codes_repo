%%

% generate random test dataset
trials = 200;

S = randi([0 3],1,trials);
R = zeros(1,trials);
R(S==0|S==1) = 1;
R2 = zeros(1,trials);
R2(S==1|S==3) = 1;
R2(S==0) = 2;


opts.bin_method_X1 = 'none';
opts.bin_method_X2 = 'none';
opts.bin_method_Y = 'none';
opts.n_bins_X1 = 4;
opts.n_bins_X2 = 4;
opts.n_bins_Y = 4;
opts.bias = 'qe';

[pid, pid_unbiased] = PID(S, R, R2, opts);

%%

% generate random test dataset
trials = 200;

S = randn(2,trials);
R = sum(S);
R2 = diff(S);


opts.bin_method_X1 = 'eqpop';
opts.bin_method_X2 = 'eqpop';
opts.bin_method_Y = 'eqpop';
opts.n_bins_X1 = 2;
opts.n_bins_X2 = 2;
opts.n_bins_Y = 2;
opts.bias = 'le';

[pid, pid_unbiased] = PID(S, R, R2, opts)