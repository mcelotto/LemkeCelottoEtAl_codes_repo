% In this test we load data from a real set of 17 neurons recorded during a 
% 2AFC task and compute II using compute_II_BROJA2

rng('default')

%load('test_data.mat')
S = [2*ones(1,50),ones(1,50);2*ones(1,50),ones(1,50)];
R = [randi(10,1,50), 5 + randi(10,1,50)];
C = double(R > 8);

opts.n_bins_stimulus = 2;
opts.n_bins_activity = 2;
opts.n_bins_choice = 2;
opts.bias = 'le';
opts.max_draws_per_split_number = 50;
opts.bin_method_stimulus = 'eqpop';
opts.bin_method_activity = 'eqpop';
opts.bin_method_choice = 'eqpop';
opts.nullhyp = 1;
opts.nh_perm = 40;
tic;
[II_biased, II_unbiased, II_null, pval] = II(S, R(1,:), C, opts);
t = toc;

disp(['The II test took ', num2str(t), ' seconds'])
disp(['II for the first neuron is ', num2str(II_biased(1))])


%% Plot
% figure()
% histogram(II_biased(2:end),20)
% xline(II_biased(1),'r','linewidth',1.5)
% xline(II_unbiased,'g','linewidth',1.5)
% lgd = legend('Null hyp distribution', 'II biased', 'II unbiased');
% lgd.FontSize = 14;
% title('FIT and its null hypothesis distribution','fontsize',16)
% xlabel('bits','fontsize',14)

% %% 2D II
% R_bins = 3;
% n1 = 2;
% n2 = 15;
% n3 = 13;
% R_2D = [R(n1,:)',R(n2,:)'];
% R_3D = [R(n1,:)',R(n2,:)',R(n3,:)'];
% II_2D = compute2D_II(S',R_2D,C',R_bins,1);
% II_per_unit(n1)+II_per_unit(n2)
% 
% II_ND2 = computeND_II(S',R_2D,C',R_bins,1);
% 
% II_ND3 = computeND_II(S',R_3D,C',R_bins,1);
