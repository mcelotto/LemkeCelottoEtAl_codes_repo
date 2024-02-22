%% MASTER SCRIPT FOR ALL INFORMATION ANALYSES
%
% Main script for data anlyses in Lemke, Celotto, et al., 2023

%% INITIALIZE PARAMETERS
    
    clear all; clc; close all;
    
    %%% dataset parameters
        params.animals = {'T102','T107','T110','T200','T201','SLDD1','SLDD2','T398'};
        params.days = {[1:8],[2:6],[2:3 5:7],[2:10],[4:17],[3:14],[4:13],[3:5 7:9]};
        params.num_earlylate_days = {[4],[3],[3],[3],[4],[3],[4],[3]};
        params.lfpFeatures = {'reach_related_lfp_ref'};
        params.areas= {'M1','DLS'};
        params.reachFeatures = {'maxPosX','maxPosY','maxVelX','maxVelY','maxAccX','maxAccY','minAccX','minAccY','maxVel','maxAcc','minAcc','movDur','distTrav','PCAstab','success_rate'}; % added success_rate only for PID_sucess_failure

    %%% rearrange and binarize parameters
        bin_params.reachFeatureBin = 3;
        bin_params.LFPBin = 5;
        bin_params.binning_pool = 'pooling_times'; % either 'pooling_times' or 'pooling_times_and_trials'
        bin_params.binSize = 50; % bin size in ms to use for reach feature computation -- length of trajectories loaded divided by binSize must be an integer
        bin_params.prevBins = 4; % how many bins before pellet touch to include for NaN discarding AND to compute "max" or "min" reach features
        bin_params.nan_thresh = 50; % min number of consecutive NaNs (in ms) ater which we discard a trial for reach feature computation
        bin_params.nan_thresh_global = 100; % min number of consecutive NaNs (in ms) ater which we discard a trial for reach feature computation
        bin_params.num_PCs = 3; % number of PCs used for PCA stability metric        
        bin_params.oldBin = 10; % TIME bin size (in ms) specified in preprocessing script
        bin_params.newBin = 1; % number of TIME bins used for rebinning neural activity 
        bin_params.totBin = (bin_params.oldBin*bin_params.newBin); % final TIME bin size
        bin_params.doLowpass = 0;
        bin_params.doHighpass = 0;
        bin_params.frequencyCutoff = 0;

    %%% mutual information parameters
        MI_params.doShuff = 1;
        MI_params.nShuff = 100;
        MI_params.opts.bias = 'pt';
        MI_params.opts.method = 'dr';
        MI_params.opts.bin_methodX = 'none';
        MI_params.opts.bin_methodY = 'none';
        MI_params.doLFP = 1;
        MI_params.doSpikes = 1; 
        MI_params.windowSize = 5; % set to 1 for Individual Points Method
        MI_params.timeJump = 2; % set to 1 for Individual Points Method
        MI_params.instMIinterval = [-1, 0.5]; % time window where we compute instantaneous MI
    
    %%% partial information decomposition parameters
        PID_params.delay = 25;
        PID_params.windowSize = 10; % set to 1 for Individual Points Method
        PID_params.timeJump = 5; % set to 1 for Individual Points Method         
        PID_params.doShuff = 1;
        PID_params.nShuff = 10;
        PID_params.opts.method = 'dr';
        PID_params.opts.bias = 'naive';
        PID_params.opts.bin_method_X1 = 'none';
        PID_params.opts.bin_method_X2 = 'none';
        PID_params.opts.bin_method_Y = 'none';
        PID_params.MIclusterparams = [.05 .95];
        PID_params.MIsigtiming = [-300 0]; % in ms, relative to pellet touch
        PID_params.doLFP=1;
        PID_params.doSpikes=1;

    %%% FIT parameters
        FIT_params.MIclusterparams = [.05 .95];
        FIT_params.MIsigtiming = [-300 0]; % in ms, relative to pellet touch

    %%% TE parameters
        TE_params.doShuff = 1;
        TE_params.nShuff = 100;
        TE_params.opts.bias = 'naive';
        TE_params.opts.method = 'dr';
        TE_params.opts.bin_methodX = 'none';
        TE_params.opts.bin_methodY = 'none';
        TE_params.opts.n_binsX = 3;
        TE_params.opts.n_binsY = 3;
        TE_params.MIclusterparams = [.05 .95];
        TE_params.MIsigtiming = [-300 0]; % in ms, relative to pellet touch
        TE_params.doLFPearlylate = 1;
        TE_params.doLFPday = 0;
        TE_params.doSpikesday = 0; 
        TE_params.windowSize = 10; % set to 1 for Individual Points Method
        TE_params.timeJump = 5; % set to 1 for Individual Points Method
        TE_params.delay = 25; 
        TE_params.saveDI = 1;

    %%% set matlab path
        restoredefaultpath;   
        main_path = 'your_path_here';
        addpath(genpath([main_path,'\code\']));
        paths.neural_activity = [main_path,'\data\preprocessed_data'];
        paths.data_path = [main_path,'\data\animal_data'];
        paths.trajectories = [main_path,'\data\preprocessed_data\all_animals_trajectories_4001_to_5000_04-Jun-2021.mat'];
        paths.trajectories_forExample = [main_path,'\data\preprocessed_data\all_animals_trajectories_4001_to_5500_08-Jun-2021.mat'];

%% COMPUTE BEHAVIOR

    %%% compute reaching features
        [reach_features] = compute_reachFeatures(bin_params,params,paths);
    
    %%% load neural activity, and rearrange/binarize both reach features and neural activity 
        [LFP_early_late, LFP_early_late_noBin, LFP_early_late_noTrialMatch, ...
         reachFeatures_early_late, reachFeatures_early_late_noBin, reachFeatures_early_late_noTrialMatch, reachFeatures_early_late_noBinNoTrialMatch...
         spikes, ...
         reachFeatures, reachFeatures_noBin, ...
         tMin, tMax] = rearrange_and_binarize(params,reach_features,bin_params,paths);

    %%% load neural activity, and rearrange/binarize both reach features and neural activity 
        [LFP_early_late, LFP_early_late_noBin, LFP_early_late_noTrialMatch, ...
         reachFeatures_early_late, reachFeatures_early_late_noBin, reachFeatures_early_late_noTrialMatch, reachFeatures_early_late_noBinNoTrialMatch...
         spikes, ...
         reachFeatures, reachFeatures_noBin, ...
         tMin, tMax] = rearrange_and_binarize_alignMovOnset(params,reach_features,bin_params,paths);
        
%% PLOT BEHAVIOR

    %%% FIGURE 1 B-E (REACH FEATURES) & SUPPLEMENTAL FIGURE 1 A-E (INDIVIDUAL ANIMAL BEHAVIOR)
    
        %%% plot behavior
            save_params.save = 0;
            save_params.save_path = [main_path,'\figures\figure_1'];
            plot_behavior(reach_features,params,save_params);
    
        %%% trajectory example
            save_params.save = 0;
            save_params.save_path = [main_path,'\figures\figure_1'];    
            plot_trajectoryExample(params,bin_params,paths,save_params);
    
        %%% glme reach features -> success models
            save_params.save = 0;
            save_params.save_path = [main_path,'\figures\figure_1'];       
            glme_success_reachFeatures(reach_features,params,save_params);
    
        %%% plot behavior match
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_3'];         
            animals_to_run = 1:8;
            features_to_run = [9 12 13];    
            plot_beh_match_lateTOearly(reachFeatures_early_late_noBinNoTrialMatch,reachFeatures_early_late_noBin,params,animals_to_run,features_to_run,save_params)
    %         plot_beh_match_earlyTOlate(reachFeatures_early_late_noBinNoTrialMatch,reachFeatures_early_late_noBin,params,animals_to_run,features_to_run,save_params)

%% COMPUTE MI

    %%% same function used for pooled and individual points method - depending on windowSize and timeJump

        % with trial match
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\MI'];    
        filtered = 0;        
        trial_match = 1;
        [MIout] = compute_mutualInformation(LFP_early_late, spikes, reachFeatures_early_late, reachFeatures, params, MI_params, PID_params, tMin, bin_params.totBin,trial_match,filtered,save_params);

        % without trial matching
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\MI'];         
        filtered = 0;
        trial_match = 0;
        [MIout] = compute_mutualInformation(LFP_early_late_noTrialMatch, spikes, reachFeatures_early_late_noTrialMatch, reachFeatures, params, MI_params, PID_params, tMin, bin_params.totBin, trial_match,filtered,save_params);

        % filtered
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\MI\filtered'];        
        trial_match = 1;        
        filtered = 1;
        [MIout] = compute_mutualInformation(LFP_early_late, spikes, reachFeatures_early_late, reachFeatures, params, MI_params, PID_params, tMin, bin_params.totBin, trial_match,filtered,save_params);

        % behavior matched
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\MI\test2'];    
        features_to_compute = [1:14];        
        [MIout] = compute_mutualInformation_behaviorMatch(LFP_early_late_noTrialMatch, spikes, reachFeatures_early_late_noTrialMatch, reachFeatures_early_late_noBinNoTrialMatch, reachFeatures, reachFeatures_noBin, params, MI_params, PID_params, tMin, bin_params.totBin, save_params, features_to_compute,bin_params);

        % success fail
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\MI'];    
        features_to_compute = [1:14];        
        [MIout] = compute_mutualInformation_successFail(LFP_early_late_noTrialMatch, reachFeatures_early_late_noTrialMatch, spikes, reachFeatures, params, MI_params, PID_params, tMin, bin_params.totBin, save_params, features_to_compute,bin_params);

        % align movement onset
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\MI\'];    
        features_to_compute = [1:14];        
        filtered = 0;
        trial_match = 1;        
        [MIout] = compute_mutualInformation_alignMovOnset(LFP_early_late, spikes, reachFeatures_early_late, reachFeatures, params, MI_params, PID_params, tMin, bin_params.totBin, trial_match,filtered,save_params,features_to_compute);

        %%% Compute instant mutual information
        features_to_compute = {'vel'};
        [MIout] = compute_instMutualInformation(LFP, spikes, spikes_ave, reachFeatures, LFP_early_late, spikes_early_late, reachFeatures_early_late, features_to_compute, params, MI_params, PID_params, rf_params, tMin, bin_params.totBin);
        
        % Compute time-lagged MI (analogous to Transfer Entropy without conditioning on receiver's past)
        animals_to_run = 1:8;
        params.reachFeatures = {'maxVel','distTrav'};
        [lagMIout] = compute_MIlagged_script_Pooled(LFP, spikes, reachFeatures, LFP_early_late, reachFeatures_early_late, params, TE_params, tMin, bin_params.totBin,animals_to_run); % using DI_params since MIlagged is very similar to TE apart from not conditioning
%% PLOT MI
    
    %%% FIGURE 1 G-H (ALL FEATURE INFORMATION) & SUPPLEMENTAL FIGURE 1 F (NAIVE VS. SKILLED)

        %%% plot pooled LFP info for ALL FEATURES
            load([processed_data_path,'\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat']);
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];   
            plot_mutualInformation_LFP_comparison_newChan_allFeatures(MIout,LFP_early_late_noBin,params,clusterParams,save_params,newBinTimes); 
    
        %%% plot spiking info for ALL FEATURES
            load([processed_data_path,'\pooled\MI\MI_01-Jul-2021_LFPearlylate_0_LFPday_0_Spikesday_1_MatchTrials_1_POOLED.mat'])
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];   
            plot_mutualInformation_Spikes_comparison_newChan_allFeatures(MIout,LFP_early_late_noBin,params,clusterParams,save_params,newBinTimes); 

        %%% scatterplot success relevance from GLME model vs. amount of information in LFP
            load([processed_data_path,'\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat']);
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];    
            success_relevance_mutual_information_scatterplot(reach_features,MIout,LFP_early_late_noBin,params,clusterParams,save_params,newBinTimes); 

        %%% scatterplot success relevance from GLME model vs. amount of information in SPIKES
            load([processed_data_path,'\pooled\MI\MI_01-Jul-2021_LFPearlylate_0_LFPday_0_Spikesday_1_MatchTrials_1_POOLED.mat'])
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];    
            success_relevance_mutual_information_scatterplot_spikes(reach_features,MIout,LFP_early_late_noBin,params,clusterParams,save_params,newBinTimes); 
            
    %%% FIGURE 2 (ACTIVITY AND INFO DISSOCIATION)
    
        %%% plot LFP info vs. LFP activity comparisons
            load([processed_data_path,'\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat']);
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];        
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];               
            plot_LFP_info_activity_comparison(MIout,LFP_early_late_noBin,reachFeatures_early_late,params,clusterParams,save_params,newBinTimes); 
    
        %%% plot spiking info vs. spiking activity comparisons
            load([processed_data_path,'\pooled\MI\MI_01-Jul-2021_LFPearlylate_0_LFPday_0_Spikesday_1_MatchTrials_1_POOLED.mat'])
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];        
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];   
            plot_spikes_info_activity_comparison(MIout,spikes,reachFeatures,params,clusterParams,save_params,newBinTimes); 
    
    %%% FIGURE 3 (ACTIVITY AND INFO TIMING)

        %%% plot pooled LFP info
            load([processed_data_path,'\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat']);
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_2'];   
            features_to_plot = [9 13];
            plot_mutualInformation_LFP_comparison_newChan(MIout,LFP_early_late_noBin,params,clusterParams,features_to_plot,save_params,newBinTimes); 
    
        %%% plot pooled spikes info
            load([processed_data_path,'\pooled\MI\MI_01-Jul-2021_LFPearlylate_0_LFPday_0_Spikesday_1_MatchTrials_1_POOLED.mat'])
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_2'];   
            features_to_plot = [9 13];
            plot_mutualInformation_Spikes_comparison_newChan(MIout,spikes,params,clusterParams,features_to_plot,save_params,newBinTimes); 

    %%% SUPPLEMENTAL FIGURE 1 G (FREQ DECOMP)

        %%% frequency decomposition pooled LFP
            broadband_path = [processed_data_path,'\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat'];
            filtered_path = [processed_data_path,'\pooled\MI\filtered'];
            load([processed_data_path,'\realBinTimes.mat'])        
            clusterParams = [.05 .95];
            features_to_plot = [9 13];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_1'];   
            plot_mutualInformation_LFP_multFreq(params,clusterParams,features_to_plot,save_params,broadband_path,filtered_path,newBinTimes); 
        
    %%% SUPPLEMENTAL FIGURE 4 A (Instant velocity info)

            % path for MI channels selection in the submitted paper
            features_to_plot = [7];  
            load([processed_data_path,'\MIinst_12-Sep-2023_shuff10_LFPearlylate_1_LFPday_0_Spikesearlylate_0_Spikesday_0_SpikesdayAve_0_MatchTrials_1'])
            bin_params.TW_pellet = 25:55; bin_params.pellet_touch_time = 50; % time frames corresponding to time window we want to plot and to time of pellet touch
            plot_params.save = 0; 
            plot_params.save_path = [];
            [delays] = plot_instMutualInformation_LFP_comparison_v3(MIout,params,features_to_plot,bin_params,MI_params,plot_params);
            
    %%% FORMER SUPPLEMENTAL FIGURE 2 (MOVEMENT ONSET ALIGNED, not included in resubmitted paper due to lack of space)
                
        %%% plot pooled LFP info movement onset aligned
            load([processed_data_path,'\pooled\MI\MI_02-Sep-2023_LFPearlylate_1_Spikes_1_MatchTrials_1_movOnsetAligned.mat']);
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_2'];   
            features_to_plot = [9 13];
            plot_mutualInformation_LFP_comparison_newChan_movOnset(MIout,LFP_early_late_noBin,params,clusterParams,features_to_plot,save_params,newBinTimes); 
    
        %%% plot pooled spikes info movement onset aligned
            load([processed_data_path,'\pooled\MI\MI_02-Sep-2023_LFPearlylate_1_Spikes_1_MatchTrials_1_movOnsetAligned.mat']);
            load([processed_data_path,'\realBinTimes.mat'])
            clusterParams = [.05 .95];
            save_params.save = 1;
            save_params.save_path = [main_path,'\figures\figure_2'];   
            features_to_plot = [9 13];
            plot_mutualInformation_Spikes_comparison_newChan_movOnset(MIout,spikes,params,clusterParams,features_to_plot,save_params,newBinTimes); 
    
            
%% COMPUTE PID
    
    %%% compute PID pooled
        animals_to_run = 1:8;
        features_to_run = [9 13];
        MI_chanSelection_LFP = [processed_data_path,'\pooled\MI\MI_24-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_0_POOLED'];
        MI_chanSelection_spikes = [processed_data_path,'\pooled\MI\MI_24-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_0_POOLED'];
        [PIDout] = compute_PID_Pooled(MI_chanSelection_LFP,MI_chanSelection_spikes,LFP_early_late, reachFeatures_early_late, spikes, reachFeatures, params, PID_params, tMin, bin_params.totBin, animals_to_run,features_to_run)

    %%% compute PID beh match
        animals_to_run = [1:2];
        features_to_run = [9 12 13];
        features_to_match = [9 12 13];
        MI_chanSelection_LFP = [processed_data_path,'\pooled\MI\MI_24-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_0_POOLED'];
        MI_chanSelection_spikes = [processed_data_path,'\pooled\MI\MI_24-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_0_POOLED'];
        save_params.save = 1;
        save_params.save_path = [processed_data_path,'\pooled\PID\PID_behMatch\lateTOearly3\'];
        [PIDout] = compute_PID_Pooled_behaviorMatch_lateTOearly(MI_chanSelection_LFP,LFP_early_late_noTrialMatch, reachFeatures_early_late, reachFeatures_early_late_noBinNoTrialMatch, params, PID_params, tMin, bin_params, animals_to_run,features_to_run,features_to_match,save_params);

    %%% compute PID individual points
        animals_to_run = 1:8;
        features_to_run = {'maxVel','distTrav'};    
        [PIDout] = compute_PID_IndPoints(LFP, spikes, reachFeatures, LFP_early_late, reachFeatures_early_late, params, PID_params, tMin, bin_params.totBin, animals_to_run);
    
    %%% compute PID success/failure
        animals_to_run = 1:8;
        features_to_run = {'maxVel','distTrav'};    
        [PIDout] = compute_PID_Pooled_successFailure(LFP_early_late, reachFeatures_early_late, params, PID_params, tMin, bin_params.totBin, animals_to_run,features_to_run);

    %%% Compute PID instant velocity
        animals_to_run = 1:8;
        features_to_run = {'vel'};
        features_to_run_static = {'maxVel'}; % used for channels selection
        PIDtime = (26:225)-175;
        % Define start and end kinematic time points
        startKinTime = -70; % in tens of ms
        endKinTime = 5;
        timeJump = 5; % tiem jump for every step (5 = 50ms)
        winLength = 10; % window size used to pool time points (10=100ms)

        t1 = find(PIDtime==startKinTime);
        t2 = find(PIDtime==endKinTime);
        kin_times = [t1-winLength/2 : t1+(winLength-1)/2]; % time point to pool for the first kin time point (time center = -705ms -> [-750,-660]);
        
        % Initialize kinematic time points for each calculation (kinematic time centers spanning [-700, 50]ms)
        for i = 1:numel(startKinTime:timeJump:endKinTime) %
            kin_times = [kin_times; t1:t1+winLength-1];
            t1 = t1+timeJump; % for every step, jump of 50ms
        end
        for i = 1:size(kin_times,1)
            [PIDout] = compute_instPID_Pooled(LFP_early_late, spikes_early_late, reachFeatures, reachFeatures_early_late, params, PID_params, tMin, bin_params.totBin, animals_to_run,features_to_run,features_to_run_static,kin_times(i,:));
        end

    %%% SUPPLEMENTAL FIGURE 2 A,C
    %%% Overlapping encoding shared info directionality simulation
        shared_info_overlap_encoding_simulation;
%% PLOT PID

%%% FIGURE 4 (SHARED INFO)

    %%% plot PID combined features
        save_params.save = 1;
        save_params.save_path = [main_path,'\figures\figure_3'];
        clusterParams = [.05 .95];
        reach_bins = 16:26;    
        params.reachFeatures = {'maxVel','distTrav'};
        PIDpath = [processed_data_path,'\pooled\PID\November2023\maxVel_distTrav_100Sh_recRef_070621_traj\'];
        plot_PID_Pooled_Sh_newChan_combinedFeatures(PIDpath,clusterParams,params,reach_bins,save_params);

%%% SUPPLEMENTAL FIGURE 2 D-G (SHARED INFO SPIKING)

    %%% spikes

        plot_params.save = 1;
        save_params.save_path = [main_path,'\figures\figure_4'];
        clusterParams = [.05 .95];
        reach_bins = 22:28;  
        PIDpath = [processed_data_path,'\pooled\PID\spiking\receiver_ref_good\'];
        plot_PID_Spikes(PIDpath,params,reach_bins,plot_params,clusterParams)
        
%%% SUPPLEMENTAL FIGURE 3 (BEHAVIOR MATCHED PID)

    %%% reach feature matched
        save_params.save = 1;
        save_params.save_path = [main_path,'\figures\figure_4'];
        clusterParams = [.05 .95];
        reach_bins = 16:26;    
        params.reachFeatures = {'maxVel','distTrav'};
        matched_features = {'maxVel','distTrav','movDur'};
        PIDpath = [processed_data_path,'\pooled\PID\PID_behMatch\'];
        plot_PID_Pooled_Sh_newChan_combinedFeatures_matchedBeh(PIDpath,clusterParams,params,reach_bins,matched_features,save_params);

    %%% matched success rate 
        save_params.save = 1;
        save_params.save_path = [main_path,'\figures\figure_4'];
        clusterParams = [.05 .95];
        reach_bins = 16:26;    
        params.reachFeatures = {'maxVel','distTrav'};
        PIDpath = [processed_data_path,'\pooled\PID\maxVel_distTrav_succFail_matched_recRef_100Shuff\'];
        plot_PID_Pooled_Sh_newChan_combinedFeatures_successMatch(PIDpath,clusterParams,params,reach_bins,save_params);

    %%% within success/failure
        plot_params.save = 0;
        plot_params.save_path = 'F:\final_IIT\figures\figure_files\fig2';
        clusterParams = [.05 .95];
        reach_bins = 16:26;
        features_to_plot = [9 13];
        PIDpath = '/media/DATA/slemke/InfoFlowProject/processed_data/PID_successFailure/sept18_maxVel_distTrav_100Sh/';
        plot_PID_Pooled_successFailure_v7(PIDpath,clusterParams,features_to_plot,params,reach_bins,plot_params);
        
%%% SUPPLEMENTAL FIGURE 4 (INSTANT VELOCITY SHARED INFO)

        %%% SUPPLEMENTAL FIGURE 4 D,E
        % Plot instant PID results for two specific kinematic time points (-500ms and -100ms from pellet touch)
        PIDtime = (26:225)-175;
        kin_time_centers = [-50,-10]; % time centers of the two kinematic times to plot; in tens of ms
        winLength = 10; % window size used to pool time points (10=100ms)
        t1 = find(PIDtime==kin_time_centers(1));
        t2 = find(PIDtime==kin_time_centers(2));
        fixedKin = [t1-winLength/2 t1+(winLength-1)/2; t2-winLength/2 t2+(winLength-1)/2]; % time centers corresponding to -500ms and -100ms respectively
        plot_params.save = 0;
        plot_params.save_path = 'F:\final_IIT\figures\figure_files\fig2';
        plot_params.time_delay = 'maxmax'; % take peak shared into in the time-delay domain
        clusterParams = [.05 .95];
        reach_bins = 16:26;    
        features_to_plot = 7;
        for i = 1:size(fixedKin,1) 
            PIDpath = ['/media/DATA/slemke/InfoFlowProject/processed_data/instPID/100Shuff_recRef/fixedKin_',num2str(fixedKin(i,1)),'_',num2str(fixedKin(i,2)),'/'];
            centerTime = 0.5*(PIDtime(fixedKin(i,1)) + PIDtime(fixedKin(i,2)));
            plot_instPID_Pooled_Sh(PIDpath,clusterParams,params,features_to_plot,reach_bins,plot_params,centerTime);
        end
        
        %%% SUPPLEMENTAL FIGURE 4 B,C
        % Plot instant PID results across several kinematic points
        startKinTime = -70; % in tens of ms
        endKinTime = 5;
        timeJump = 5; % tiem jump for every step (5 = 50ms)
        winLength = 10; % window size used to pool time points (10=100ms)

        t1 = find(PIDtime==startKinTime);
        t2 = find(PIDtime==endKinTime);
        fixedKin = [t1-winLength/2 : t1+(winLength-1)/2]; % time point to pool for the first kin time point (time center = -705ms -> [-750,-660]);
        
        kinTimePoints = numel(startKinTime:timeJump:endKinTime); % number of kinematic time points
        fixedKin = repelem(fixedKin,kinTimePoints,1);
        fixedKin = fixedKin + (0:5:(kinTimePoints-1)*5)';
        plot_instPID_Pooled_Sh_multi_timePoints(clusterParams,params,features_to_plot,reach_bins,plot_params,fixedKin);
        

%%% SUPPLEMENTAL FIGURE 5 (ALL FEATURE INFO)
    
    %%% all feature PID

        save_params.save = 1;
        save_params.save_path = [main_path,'\figures\figure_4'];
        clusterParams = [.05 .95];
        reach_bins = 16:26;
        params.reachFeatures = {'maxPosX','maxPosY','maxVelX','maxVelY','maxAccX','maxAccY','minAccX','minAccY','maxVel','maxAcc','minAcc','movDur','distTrav','PCAstab','success_rate'}; % added success_rate only for PID_sucess_failure        
        plot_PID_Pooled_Sh_all_features(clusterParams,params,reach_bins,save_params);

%%% SUPPLEMENTAL FIGURE 3 S (success and # trials contribution to reversal)
        data_path = []
        neural_metric = 'info flow';
        success_n_trials_contribution_to_reversal(data_path,neural_metric)

%% PLOT CROSS-AREA INFO AND LFP COH

%%% SUPPLEMENTAL FIGURE 2 H-M

    %%% plot cross-area MI across all channels
    plot_params.save = 1;
    plot_params.save_path = [main_path,'\figures\figure_4'];
    clusterParams = [.05 .95];
    reach_bins = 16:26;
    infoPath = [processed_data_path,'\pooled\cross_area_MI\'];
    plot_crossAreaMI(infoPath,clusterParams,params,FIT_params,reach_bins,plot_params)

    %%% plot LFP coherence
    plot_params.save = 1;
    plot_params.save_path = [main_path,'\figures\figure_4'];
    coh_params.Fs= 100;
    coh_params.tapers=[3 5];
    coh_params.fpass = [1 30];
    coh_params.err = [2 0.05];
    coh_params.trialave = 1;
    coh_params.pad = 0;
    coh_params.movingwin = [.5 .05];
    plot_LFPcoherence(params,coh_params,plot_params)    

%% Compute Transfer Entropy (TE)
%%% compute TE pooled
        animals_to_run = 1:8;
        [TEout] = compute_TE_script_Pooled(LFP, spikes, reachFeatures, LFP_early_late, reachFeatures_early_late, params, TE_params, tMin, bin_params.totBin,animals_to_run);
      
%% PLOT TE

%%% FIGURE 5 E,F

    %%% plot TE
    plot_params.save = 1;
    plot_params.save_path = [main_path,'\figures\figure_4'];
    clusterParams = [.05 .95];
    reach_bins = 16:26;
    TEpath = [processed_data_path,'\pooled\DI\'];
    plot_TE_Pooled_Sh_newChan_combinedFeatures(TEpath,clusterParams,params,FIT_params,reach_bins,plot_params)    
    
%% PLOT FIT

%%% FIGURE 5 C,D

    %%% plot FIT combined features
    plot_params.save = 1;
    plot_params.save_path = [main_path,'\figures\figure_4'];
    clusterParams = [.05 .95];
    reach_bins = 16:26;
    params.reachFeatures = {'maxVel','distTrav'};
    FITpath = [processed_data_path,'\pooled\FIT\'];
    plot_FIT_Pooled_Sh_newChan_combinedFeatures(FITpath,clusterParams,params,reach_bins,plot_params)

%%% SUPPLEMENTAL FIGURE 7 (ALL FEATURE INFO)
    
    %%% all feature FIT
    save_params.save = 1;
    save_params.save_path = [main_path,'\figures\figure_4'];
    clusterParams = [.05 .95];
    reach_bins = 16:26;
    params.reachFeatures = {'maxPosX','maxPosY','maxVelX','maxVelY','maxAccX','maxAccY','minAccX','minAccY','maxVel','maxAcc','minAcc','movDur','distTrav','PCAstab','success_rate'}; % added success_rate only for PID_sucess_failure        
    plot_FIT_Pooled_Sh_all_features(clusterParams,params,reach_bins,save_params);
