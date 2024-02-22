%% Content affects directionality of transfer
% Communication in presence of overlapping information timecourses
% Real transmission from X to Y only


clear all; %close all;

rng(1) % For reproducibility
save_results = 0; % set to 1 to save results file
results_path =  '\your_results_path'; % path to the directory to save results

% Simulation parameters
nTrials_per_stim = 500; % number of trials per stimulus value
simReps = 100; % repetitions of the simulation
nShuff = 10; % number of permutations (used for both FIT permutation tests)

% w_xy_noise = 0; % range of w_noise parameter
noise = 5; % standard deviation of gaussian noise in X_noise and Y

f_encoding_hX = [3 1.5 1];
f_encoding_Xt = [3 2.5 1];
f_encoding_hY = [1 3 1.5];
f_encoding_Yt = [3 1.5 1];

% Define information options
opts = [];
opts.verbose = false;
opts.method = "dr";
opts.bias = 'naive';
opts.btsp = 0;
opts.n_binsX = 5;
opts.n_binsY = 5; 
opts.n_binsC = 3; % Number of stimulus values
opts.bin_method_X = 'none';
opts.bin_method_Y = 'none';
opts.bin_method_C = 'none';

% Initialize structures
fit_xy = nan(1,simReps); di_xy = fit_xy; shared_info_xy = fit_xy; 
fit_yx = nan(1,simReps); di_yx = fit_xy; shared_info_yx = fit_xy; 
fitSh_xy.simple = nan(1,simReps,nShuff); diSh_xy.simple = fitSh_xy.simple; 
fitSh_yx.simple = nan(1,simReps,nShuff); diSh_yx.simple = fitSh_xy.simple; 
info.hX = fit_xy; info.Xt = fit_xy; info.hY = fit_xy; info.Yt = fit_xy;

%% Plot encoding functions
nTrials = nTrials_per_stim(1)*opts.n_binsC; % Compute number of trials
S = randi(opts.n_binsC,1,nTrials);

figure('Position',[4,311,1274,279])
subplot(1,4,1)
encoding_function(S, f_encoding_hX, 1, 1);
xlim([0.5 3.5]) 
xlabel('Max velocity')
xticklabels({'Slow','Mid','Fast'})
title('X_{past}')
subplot(1,4,2)
encoding_function(S, f_encoding_hY, 1, 1);
xlim([0.5 3.5])
xlabel('Max velocity') 
xticklabels({'Slow','Mid','Fast'})
title('Y_{past}')
subplot(1,4,3)
encoding_function(S, f_encoding_Xt, 1, 1);
xlim([0.5 3.5]) 
xlabel('Max velocity') 
xticklabels({'Slow','Mid','Fast'})
title('X_{pres}')
subplot(1,4,4)
encoding_function(S, f_encoding_Yt, 1, 1);
xlim([0.5 3.5]) 
xlabel('Max velocity') 
xticklabels({'Slow','Mid','Fast'})
title('Y_{pres}')

%% Run simulation

    
for trialIdx = 1:numel(nTrials_per_stim)
    disp(['Simulation for ',num2str(nTrials_per_stim(trialIdx)),' trials per stim'])
    for repIdx = 1:simReps
        %disp(['Repetition number ',num2str(repIdx)]);
        nTrials = nTrials_per_stim(trialIdx)*opts.n_binsC; % Compute number of trials

        % Draw the stimulus value for each trial
        S = randi(opts.n_binsC,1,nTrials);

        % add measurement noise
        hX = encoding_function(S,f_encoding_hX,1,0) + noise*randn(1,nTrials); 
        hY = encoding_function(S,f_encoding_hY,1,0) + noise*randn(1,nTrials); 
        Xt = encoding_function(S,f_encoding_Xt,1,0) + noise*randn(1,nTrials);
        Yt = hX;
        
        % Discretize neural activity
        bXpast = eqpop(hX, opts.n_binsX);
        bXt = eqpop(Xt, opts.n_binsX);

        bYt = eqpop(Yt, opts.n_binsY);
        bYpast = eqpop(hY, opts.n_binsY);

        [fit_xy(repIdx)]=...
            FIT_choice(S,bXpast,bYpast,bYt,opts);
        [fit_yx(repIdx)]=...
            FIT_choice(S,bYpast,bXpast,bXt,opts);
        
        tmpPID = PID(S,bXpast,bYt,opts);
        [shared_info_xy(repIdx)] = tmpPID.shared;
        
        tmpPID = PID(S,bYpast,bXt,opts);
        [shared_info_yx(repIdx)] = tmpPID.shared;
        
        info.hX(repIdx) = cell2mat(information(S,bXpast,opts,{'I'}));
        info.Xt(repIdx) = cell2mat(information(S,bXt,opts,{'I'}));
        info.hY(repIdx) = cell2mat(information(S,bYpast,opts,{'I'}));
        info.Yt(repIdx) = cell2mat(information(S,bYt,opts,{'I'}));
    end
end


if save_results
    fname = ['simulation_results.mat'];
    save([fname])
end

%% Plot 

x_vals = [1 2];

figure('Position',[282,268,409,281])
hold on
bar(x_vals(1),mean(fit_xy),'FaceColor','b');
bar(x_vals(2),mean(fit_yx),'FaceColor','r');
errorbar(x_vals,[mean(fit_xy) mean(fit_yx)], [std(fit_xy)/sqrt(simReps) std(fit_yx)/sqrt(simReps)],'ko');
set(gca, 'XTick',x_vals, 'XTickLabel', {'X -> Y', 'Y -> X'})

title(['FIT due to encoding format difference'])


figure('Position',[282,268,409,281])
hold on
bar(x_vals(1),mean(shared_info_xy),'FaceColor','b');
bar(x_vals(2),mean(shared_info_yx),'FaceColor','r');
errorbar(x_vals,[mean(shared_info_xy) mean(shared_info_yx)], [std(shared_info_xy)/sqrt(simReps) std(shared_info_yx)/sqrt(simReps)],'ko');
set(gca, 'XTick',x_vals, 'XTickLabel', {'X -> Y', 'Y -> X'})

title(['Shared info due to encoding format difference'])

cols = distinguishable_colors(4);

figure('Position',[282,268,409,281])
hold on
bar(1,mean(info.hX),'FaceColor',cols(1,:));
bar(2,mean(info.Xt),'FaceColor',cols(2,:));
bar(3,mean(info.hY),'FaceColor',cols(3,:));
bar(4,mean(info.Yt),'FaceColor',cols(4,:));
errorbar([1 2 3 4],[mean(info.hX) mean(info.Xt) mean(info.hY) mean(info.Yt)], [std(info.hX)/sqrt(simReps) std(info.Xt)/sqrt(simReps) std(info.hY)/sqrt(simReps) std(info.Yt)/sqrt(simReps)],'ko');
set(gca, 'XTick',[1 2 3 4], 'XTickLabel', {'X_{past}', 'X_{pres}','Y_{past}','Y_{pres}'})
title(['Info at different time points'])
