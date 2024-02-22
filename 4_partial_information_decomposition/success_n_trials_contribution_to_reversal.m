%% Run GLM within naive and skilled days to predict info directionality from success rate
% Select params to plot and compute pooled across-feature mean metrics
% Within learning stage correlation between reversal and success

function success_n_trials_contribution_to_reversal(data_path,neural_metric)

    clear all; close all
    
    load(data_path)
    all_chans_flow.M1_DLS_early = all_chans_flow.M1_DLS.early;
    all_chans_flow = rmfield(all_chans_flow,'M1_DLS');
    
    rng(0) % for reproducibility
    
    plot_corr_animals = 0;
    % delays windows corresponding to M1-to-DLS and DLS-to-M1 directions
    M1_DLS_lags = [27 51];
    DLS_M1_lags = [1 25];
    
    % Plot correlations across subjects params
    % neural_metric = 'info flow'; % either 'info flow' or 'fract pairs' or 'delays'
    sel_reach_feature = 'pooled';
    analysis_animals = [1 2 3 4 5 6 7 8];
    dpt = 0.1; scale = 1;
    
    % GLM params
    glm_reach_feat = 'pooled';
    % preds_labels = {'# train trials','succ. rate'};
    remove_out_window_pairs = 1; % if 1, remove all pairs not belonging to either M1->DLS or DLS->M1 delay window (i.e. zero lag)
    nResamplings = 10000;
    confidence_interv_prc = [2.5,97.5];
    cv_folds = 5;
    
    sigmoid = @(x) 1 ./ (1 + exp(-x));
    
    % Compute change lags
    for fidx = 1:numel(params.PIDreachFeatures)
        start_pair_early = 1;
        start_pair_late = 1;
    
        chan_pairs_lag.(params.PIDreachFeatures{fidx}) = [];
        chan_pairs_lag.(params.PIDreachFeatures{fidx}) = [];
        condition_info.(params.PIDreachFeatures{fidx}) = [];
        early_chans_animal.(params.PIDreachFeatures{fidx}) = [];
        late_chans_animal.(params.PIDreachFeatures{fidx}) = [];
        all_chans_flow.M1_DLS.(params.PIDreachFeatures{fidx})=[];
        all_chans_flow.DLS_M1.(params.PIDreachFeatures{fidx})=[];
        chans_delays.(params.PIDreachFeatures{fidx}) = [];
    
        change_lag.(params.PIDreachFeatures{fidx}).M1_DLS = [];
        change_lag.(params.PIDreachFeatures{fidx}).DLS_M1 = [];
        for animal = 1:numel(params.animals)
            fract_lag.(params.PIDreachFeatures{fidx}).early.M1_DLS(animal) = sum(early_peakDelay{fidx,animal}>=M1_DLS_lags(1) & early_peakDelay{fidx,animal}<=M1_DLS_lags(2))/length(early_peakDelay{fidx,animal});
            fract_lag.(params.PIDreachFeatures{fidx}).late.M1_DLS(animal)  = sum(late_peakDelay{fidx,animal}>=M1_DLS_lags(1) & late_peakDelay{fidx,animal}<=M1_DLS_lags(2))/length(late_peakDelay{fidx,animal});
            fract_lag.(params.PIDreachFeatures{fidx}).early.DLS_M1(animal)  = sum(early_peakDelay{fidx,animal}>=DLS_M1_lags(1) & early_peakDelay{fidx,animal}<=DLS_M1_lags(2))/length(early_peakDelay{fidx,animal});
            fract_lag.(params.PIDreachFeatures{fidx}).late.DLS_M1(animal)  = sum(late_peakDelay{fidx,animal}>=DLS_M1_lags(1) & late_peakDelay{fidx,animal}<=DLS_M1_lags(2))/length(late_peakDelay{fidx,animal});
    
    
            delays_animal.(params.PIDreachFeatures{fidx}).early(animal) = nanmean(-(early_peakDelay{fidx,animal}-26))*10;
            delays_animal.(params.PIDreachFeatures{fidx}).late(animal) = nanmean(-(late_peakDelay{fidx,animal}-26))*10;
            % Compute binarized channel pairs lags for GLM
            EL_chans_pairs = [early_peakDelay{fidx,animal},late_peakDelay{fidx,animal}];
            binarized_chans_lags = zeros(1,numel(EL_chans_pairs));
            binarized_chans_lags((EL_chans_pairs >=M1_DLS_lags(1) & EL_chans_pairs <=M1_DLS_lags(2))) = -1;
            binarized_chans_lags((EL_chans_pairs >=DLS_M1_lags(1) & EL_chans_pairs <=DLS_M1_lags(2))) = 1;
            chan_pairs_lag.(params.PIDreachFeatures{fidx}) = [chan_pairs_lag.(params.PIDreachFeatures{fidx}), binarized_chans_lags];
            condition_info.(params.PIDreachFeatures{fidx}) = [condition_info.(params.PIDreachFeatures{fidx}), zeros(1,numel(early_peakDelay{fidx,animal})),ones(1,numel(late_peakDelay{fidx,animal}))];
    
            early_chans_animal.(params.PIDreachFeatures{fidx})(animal) = numel(early_peakDelay{fidx,animal});
            late_chans_animal.(params.PIDreachFeatures{fidx})(animal) = numel(late_peakDelay{fidx,animal});
    
            early_pairs = start_pair_early:start_pair_early+early_chans_animal.(params.PIDreachFeatures{fidx})(animal)-1;
            late_pairs = start_pair_late:start_pair_late+late_chans_animal.(params.PIDreachFeatures{fidx})(animal)-1;
    
            % Info flow directionality magnitude
            all_chans_flow.M1_DLS.(params.PIDreachFeatures{fidx}) = [all_chans_flow.M1_DLS.(params.PIDreachFeatures{fidx}), all_chans_flow.M1_DLS_early.(params.PIDreachFeatures{fidx})(early_pairs)', all_chans_flow.M1_DLS_late.(params.PIDreachFeatures{fidx})(late_pairs)'];
            all_chans_flow.DLS_M1.(params.PIDreachFeatures{fidx}) = [all_chans_flow.DLS_M1.(params.PIDreachFeatures{fidx}), all_chans_flow.DLS_M1_early.(params.PIDreachFeatures{fidx})(early_pairs)', all_chans_flow.DLS_M1_late.(params.PIDreachFeatures{fidx})(late_pairs)'];
    
            start_pair_early = start_pair_early+early_chans_animal.(params.PIDreachFeatures{fidx})(animal);
            start_pair_late = start_pair_late+late_chans_animal.(params.PIDreachFeatures{fidx})(animal);
    
            % Delay per channel
            chans_delays.(params.PIDreachFeatures{fidx}) = [chans_delays.(params.PIDreachFeatures{fidx}), (-(early_peakDelay{fidx,animal}-26))*10, (-(late_peakDelay{fidx,animal}-26))*10];
        end         
        change_lag.(params.PIDreachFeatures{fidx}).M1_DLS = fract_lag.(params.PIDreachFeatures{fidx}).late.M1_DLS-fract_lag.(params.PIDreachFeatures{fidx}).early.M1_DLS;
        change_lag.(params.PIDreachFeatures{fidx}).DLS_M1 = fract_lag.(params.PIDreachFeatures{fidx}).late.DLS_M1-fract_lag.(params.PIDreachFeatures{fidx}).early.DLS_M1;
        change_lag.(params.PIDreachFeatures{fidx}).early = fract_lag.(params.PIDreachFeatures{fidx}).early.M1_DLS-fract_lag.(params.PIDreachFeatures{fidx}).early.DLS_M1;
        change_lag.(params.PIDreachFeatures{fidx}).late = fract_lag.(params.PIDreachFeatures{fidx}).late.M1_DLS-fract_lag.(params.PIDreachFeatures{fidx}).late.DLS_M1;
    end
    
    animalColors = distinguishable_colors(numel(params.animals));
    direct_flow_animals.pooled.early = [];
    direct_flow_animals.pooled.late = [];
    change_lag.pooled.early = [];
    change_lag.pooled.late = [];
    change_lag.pooled.M1_DLS = [];
    change_lag.pooled.DLS_M1 = [];
    delays_animal.pooled.early = [];
    delays_animal.pooled.late = [];
    
    direct_flow_animals.mean.early = [];
    direct_flow_animals.mean.late = [];
    change_lag.mean.early = [];
    change_lag.mean.late = [];
    change_lag.mean.M1_DLS = [];
    change_lag.mean.DLS_M1 = [];
    delays_animal.mean.early = [];
    delays_animal.mean.late = [];
    
    for fidx = 1:numel(params.PIDreachFeatures)
        direct_flow_animals.pooled.early = [direct_flow_animals.pooled.early,direct_flow_animals.(params.PIDreachFeatures{fidx}).early];
        direct_flow_animals.pooled.late = [direct_flow_animals.pooled.late,direct_flow_animals.(params.PIDreachFeatures{fidx}).late];
        change_lag.pooled.early = [change_lag.pooled.early,change_lag.(params.PIDreachFeatures{fidx}).early];
        change_lag.pooled.late = [change_lag.pooled.late,change_lag.(params.PIDreachFeatures{fidx}).late];
        change_lag.pooled.M1_DLS = [change_lag.pooled.M1_DLS,change_lag.(params.PIDreachFeatures{fidx}).M1_DLS];
        change_lag.pooled.DLS_M1 = [change_lag.pooled.DLS_M1,change_lag.(params.PIDreachFeatures{fidx}).DLS_M1];
    
        delays_animal.pooled.early = [delays_animal.pooled.early, delays_animal.(params.PIDreachFeatures{fidx}).early];
        delays_animal.pooled.late = [delays_animal.pooled.late, delays_animal.(params.PIDreachFeatures{fidx}).late];
    
        direct_flow_animals.mean.early = [direct_flow_animals.mean.early;direct_flow_animals.(params.PIDreachFeatures{fidx}).early];
        direct_flow_animals.mean.late = [direct_flow_animals.mean.late;direct_flow_animals.(params.PIDreachFeatures{fidx}).late];
        change_lag.mean.early = [change_lag.mean.early;change_lag.(params.PIDreachFeatures{fidx}).early];
        change_lag.mean.late = [change_lag.mean.late;change_lag.(params.PIDreachFeatures{fidx}).late];
        change_lag.mean.M1_DLS = [change_lag.mean.M1_DLS;change_lag.(params.PIDreachFeatures{fidx}).M1_DLS];
        change_lag.mean.DLS_M1 = [change_lag.mean.DLS_M1;change_lag.(params.PIDreachFeatures{fidx}).DLS_M1];
    
        delays_animal.mean.early = [delays_animal.mean.early; delays_animal.(params.PIDreachFeatures{fidx}).early];
        delays_animal.mean.late = [delays_animal.mean.late; delays_animal.(params.PIDreachFeatures{fidx}).late];
        
    end
    direct_flow_animals.mean.early = mean(direct_flow_animals.mean.early,1);
    direct_flow_animals.mean.late = mean(direct_flow_animals.mean.late,1);
    change_lag.mean.early = mean(change_lag.mean.early,1);
    change_lag.mean.late = mean(change_lag.mean.late,1);
    delays_animal.mean.early = mean(delays_animal.mean.early,1);
    delays_animal.mean.late = mean(delays_animal.mean.late,1);
    
    if strcmp(neural_metric,'info flow')
        direct_metric = direct_flow_animals;
        direct_metric_label = 'M1->DLS-DLS->M1 info [bits]';
        title_direct_metric = 'Info flow direction ';
        GLM_distribution = 'normal';
        GLM_link = 'identity';
    elseif strcmp(neural_metric,'fract pairs')
        direct_metric = change_lag;
        direct_metric_label = 'fraction M1->DLS-DLS->M1 [%]';
        title_direct_metric = 'Fract. directed pairs ';
        GLM_distribution = 'binomial';
        GLM_link = 'logit';
    elseif strcmp(neural_metric,'delays')
        direct_metric = delays_animal;
        direct_metric_label = 'delay [ms]';
        title_direct_metric = 'Delay ';
        GLM_distribution = 'normal';
        GLM_link = 'identity';
    end
    
    % delta_training_days = avg_training_day(:,2)-avg_training_day(:,1);
    train_trials(:,1) = mean_trials(:,1);
    for animal = 1:numel(params.animals)
       train_trials(animal,2) = sum(n_trials_per_day.maxVel{animal}) - mean_trials(animal,2);
    end
    
    tmp_training_trials = train_trials;
    delta_training_trials = train_trials(:,2)-train_trials(:,1);
    
    % repmat delta training
    % repmat delta taining and behav
    
    learn_stage = [zeros(1,numel(params.animals)),ones(1,numel(params.animals))];
    if strcmp(sel_reach_feature,'pooled')
        delta_training_trials = repmat(delta_training_trials,numel(params.PIDreachFeatures),1);
        mean_reachFeat = repmat(mean_reachFeat,1,numel(params.PIDreachFeatures),1);
        animalColors = repmat(animalColors,numel(params.PIDreachFeatures),1);
        learn_stage = repelem(learn_stage,1,2);
        training_days = [];
        analysis_animals = [analysis_animals analysis_animals+8];
        for fidx = 1:numel(params.PIDreachFeatures)
            training_days = [training_days;tmp_training_trials];
        end
    else
        training_days = tmp_training_trials;
    end
    
    % Build GLM predictors
    succ_idx = find(strcmp(params.reachFeatures,'success_rate'));
    for fidx = 1:numel(params.PIDreachFeatures)
        training_days_GLM.(params.PIDreachFeatures{fidx}) = [];
        animal_chan_pair_GLM.(params.PIDreachFeatures{fidx}) = [];
        success_rate_GLM.(params.PIDreachFeatures{fidx}) = [];
        for animal = 1:numel(params.animals)
            training_days_GLM.(params.PIDreachFeatures{fidx}) = [training_days_GLM.(params.PIDreachFeatures{fidx}) repelem(tmp_training_trials(animal,1),1,early_chans_animal.(params.PIDreachFeatures{fidx})(animal))];
            training_days_GLM.(params.PIDreachFeatures{fidx}) = [training_days_GLM.(params.PIDreachFeatures{fidx}) repelem(tmp_training_trials(animal,2),1,late_chans_animal.(params.PIDreachFeatures{fidx})(animal))];
    
            success_rate_GLM.(params.PIDreachFeatures{fidx}) = [success_rate_GLM.(params.PIDreachFeatures{fidx}) repelem(mean_reachFeat(succ_idx,animal,1),1,early_chans_animal.(params.PIDreachFeatures{fidx})(animal))];
            success_rate_GLM.(params.PIDreachFeatures{fidx}) = [success_rate_GLM.(params.PIDreachFeatures{fidx}) repelem(mean_reachFeat(succ_idx,animal,2),1,late_chans_animal.(params.PIDreachFeatures{fidx})(animal))];
    
            animal_chan_pair_GLM.(params.PIDreachFeatures{fidx}) = [animal_chan_pair_GLM.(params.PIDreachFeatures{fidx}) repelem(animal,1,early_chans_animal.(params.PIDreachFeatures{fidx})(animal))];
            animal_chan_pair_GLM.(params.PIDreachFeatures{fidx}) = [animal_chan_pair_GLM.(params.PIDreachFeatures{fidx}) repelem(animal,1,late_chans_animal.(params.PIDreachFeatures{fidx})(animal))];
        end
    end
        
    %% GLM across channels pairs
    
    all_chans_flow_diff.maxVel = (all_chans_flow.DLS_M1.maxVel-all_chans_flow.M1_DLS.maxVel);
    all_chans_flow_diff.distTrav = (all_chans_flow.DLS_M1.distTrav-all_chans_flow.M1_DLS.distTrav);
    
    chan_pairs_lag.pooled = [chan_pairs_lag.maxVel, chan_pairs_lag.distTrav];
    chans_delays.pooled = [chans_delays.maxVel, chans_delays.distTrav];
    all_chans_flow_diff.pooled = [all_chans_flow_diff.maxVel, all_chans_flow_diff.distTrav];
    condition_info.pooled = [condition_info.maxVel, condition_info.distTrav];
    success_rate_GLM.pooled = [success_rate_GLM.maxVel, success_rate_GLM.distTrav];
    training_days_GLM.pooled = [training_days_GLM.maxVel, training_days_GLM.distTrav];
    animal_chan_pair_GLM.pooled = [animal_chan_pair_GLM.maxVel, animal_chan_pair_GLM.distTrav];
    
    if remove_out_window_pairs
        if strcmp(neural_metric,'fract pairs') % for fract pairs we remove zero lag pairs
            rem_idxs = chan_pairs_lag.(glm_reach_feat)==0;
            chan_pairs_lag.(glm_reach_feat)(rem_idxs) = [];
            condition_info.(glm_reach_feat)(rem_idxs) = [];
            success_rate_GLM.(glm_reach_feat)(rem_idxs) = [];
            training_days_GLM.(glm_reach_feat)(rem_idxs) = [];
            animal_chan_pair_GLM.(glm_reach_feat)(rem_idxs) = [];
        end
        condition_info.(glm_reach_feat)(condition_info.(glm_reach_feat)==0) = -1; % Code naive condition as -1
        % Code M1->DLS as zero
        chan_pairs_lag.(glm_reach_feat)(chan_pairs_lag.(glm_reach_feat)==-1) = 0;
    else % Set M1->DLS lags to same value as out window
        chan_pairs_lag.(glm_reach_feat)(chan_pairs_lag.(glm_reach_feat)==-1) = 0;
    end
    
    predictors = [condition_info.(glm_reach_feat);training_days_GLM.(glm_reach_feat);success_rate_GLM.(glm_reach_feat)]';
    if strcmp(neural_metric,'info flow')
        y = all_chans_flow_diff.(glm_reach_feat)';
    elseif strcmp(neural_metric,'fract pairs')
        y = chan_pairs_lag.(glm_reach_feat)';
    elseif strcmp(neural_metric,'delays')
        y = chans_delays.(glm_reach_feat)';
    end
    
    
    ealy_late_label = {'naive','skilled'};
    for early_late = 1:2
        elLabel = ealy_late_label{early_late};
        betas.(elLabel) = nan(3,nResamplings+1,cv_folds);

        learn_stage_idx = condition_info.pooled==(early_late-1.5)*2;
        yStage = y(learn_stage_idx);
        predsStage = predictors(learn_stage_idx,:);
    
        % Z-score predictors
        % condition_info.(sel_reach_feat) = zscore(condition_info.(sel_reach_feat));
        predsStage(:,2) = zscore(predsStage(:,2));
        predsStage(:,3) = zscore(predsStage(:,3));
    
    
        for resIdx = 1:nResamplings+1
            resIdx
            if resIdx == 1 % original data
                indexes = 1:numel(yStage);
            else % random resampling with replacements
                indexes = randsample(1:numel(yStage),numel(yStage),true);
            end
        
            bootY = yStage(indexes);
            bootPreds = predsStage(indexes,:);
            cv = cvpartition(yStage,'KFold',cv_folds); % Should we have a more balanced cv?
            
            % Fit null model (no predictors to optimize with train/test
            [B, FitInfo_null] = lassoglm(ones(1,numel(bootY))', bootY, GLM_distribution, 'Link', GLM_link);
            B0 = B(1);
        
            for cvi = 1:length(cv.TrainSize)
                D_trn = bootPreds(cv.training(cvi),:);
                D_test = bootPreds(cv.test(cvi),:);
            
                bootY_trn = bootY(cv.training(cvi));
                y_meas = bootY(cv.test(cvi));
                    
                if strcmp(GLM_distribution,'binomial')
                    y_pred_null = sigmoid(linear_pred_null); % Poisson regression with log link
                    % y_pred_null(y_pred_null == 0) = eps;
                    dev_null = -2 * sum(y_meas .* log(y_pred_null) + (1 -y_meas) .* log(1-y_pred_null));
                    tmp_accu_null = mean(y_pred_null>0.5==y_meas);
                    accuracy.(elLabel)(1,resIdx,cvi) = tmp_accu_null;
                elseif strcmp(GLM_distribution,'normal')
                    y_pred_null = linear_pred_null; % Poisson regression with log link
                    % y_pred_null(y_pred_null == 0) = eps;
                    dev_null = sum((y_meas - y_pred_null).^2,'omitnan');
                    tmp_accu_null = mean(y_pred_null>0==y_meas>0);
                    accuracy.(elLabel)(1,resIdx,cvi) = tmp_accu_null;
                end
    
                b = [];
                preds_idx = [2 3];
                [B, FitInfo] = glmfit(D_trn(:,preds_idx), bootY_trn, GLM_distribution, 'Link', GLM_link);
                b = B;
            
                betas.(elLabel)(:,resIdx,cvi) = b;
            end
        end
    end
    
    %% Plot bias and naive vs skilled coefficients
    figure()
    hold on
    % Naive
    bootBetas = squeeze(mean(betas.naive(:,2:end,:),3));
    lowPrctile = squeeze(prctile(bootBetas,confidence_interv_prc(1),2));
    highPrctile = squeeze(prctile(bootBetas,confidence_interv_prc(2),2));
    originalBetas = squeeze(mean(betas.naive(:,1,:),3));
    pvals = min([sum(bootBetas<0,2),sum(bootBetas>0,2)]/nResamplings,[],2);
    pvals = pvals*2; % since it's a two tailed ttest
    
    bar([1 2 3],originalBetas)
    errorbar([1 2 3],originalBetas,originalBetas-lowPrctile,highPrctile-originalBetas,'ok')
    for pred = 1:3
        pvalues_plot(pvals(pred),pred,(originalBetas(pred)+highPrctile(pred)),dpt,scale,14,0,'k',0)
    end
    
    % Skilled
    bootBetas = squeeze(mean(betas.skilled(:,2:end,:),3));
    lowPrctile = squeeze(prctile(bootBetas,confidence_interv_prc(1),2));
    highPrctile = squeeze(prctile(bootBetas,confidence_interv_prc(2),2));
    originalBetas = squeeze(mean(betas.skilled(:,1,:),3));
    pvals = min([sum(bootBetas<0,2),sum(bootBetas>0,2)]/nResamplings,[],2);
    pvals = pvals*2; % since it's a two tailed ttest
    
    bar([5 6 7],originalBetas)
    errorbar([5 6 7],originalBetas,originalBetas-lowPrctile,highPrctile-originalBetas,'ok')
    for pred = 1:3
        pvalues_plot(pvals(pred),pred+4,(originalBetas(pred)+highPrctile(pred)),dpt,scale,14,0,'k',0)
    end
    
    ylabel('beta')
    title([neural_metric ' model coefficients'])
    
    set(gca, 'XTick',[1 2 3 5 6 7], 'XTickLabel',{'bias naive','# trials naive','succ naive','bias skilled','# trials skilled','success skilled'})
end
