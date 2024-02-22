function [PIDout] = compute_PID_Pooled_behaviorMatch(MI_chanSelection_LFP,LFP_early_late, reachFeatures_early_late, reachFeatures_early_late_noBinNoTrialMatch, params, PID_params, tMin, bin_params, animals_to_run,features_to_run)
   
%% calculate LFP PID by early/late
if PID_params.doLFP==1

    MI_info = load(MI_chanSelection_LFP);

    for animal = animals_to_run
        
        clearvars PIDout
        
        % determine channels with significant mutual information
        sig_channels = cell(length(params.lfpFeatures),length(params.reachFeatures),length(params.areas),2);
        for n_lfp = 1:length(params.lfpFeatures)
            for fidx = 1:length(params.reachFeatures)-1
                for aidx = 1:length(params.areas)
                    for early_late = 1:2
                        sig_chan = [];
                        for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                            infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                            infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,PID_params.MIclusterparams(1),PID_params.MIclusterparams(2));
                            infQuant(~tmp_sig) = 0;
                            [~,i] = max(infQuant);
                            sigTime_start = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(1));
                            sigTime_end = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(2));
                            sig_chan = [sig_chan i>=sigTime_start & i<=sigTime_end];
                        end
                        sig_channels{n_lfp,fidx,aidx,early_late} = sig_chan;
                    end
                end
            end
        end
        
        for fidx = features_to_run

            if fidx==9
                dist_thresh = .5;
            elseif fidx==13
                dist_thresh = 5;
            elseif fidx==14
                dist_thresh = .2;
            end

            %%% match late TO early
            all_late_vals = reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx});
            matched_trial_vals = cell(1,2);
            matched_trials = cell(1,2);
            for i = 1:length(reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx}))
                tmp_early_val = reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx})(i);
                [val, idx] = min(abs(all_late_vals-tmp_early_val));
                if val<dist_thresh
                    matched_trial_vals{1} = [matched_trial_vals{1} tmp_early_val];
                    matched_trial_vals{2} = [matched_trial_vals{2} all_late_vals(idx)];
                    matched_trials{1} = [matched_trials{1} i];
                    matched_trials{2} = [matched_trials{2} idx];
                end
                all_late_vals(idx) = NaN;
            end

            tmp_min = min([reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx}) reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx})]);
            tmp_max = max([reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx}) reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx})]);
            rangeX = [tmp_min:(tmp_max-tmp_min)/20:tmp_max];
            figure;
                subplot(1,2,1); hold on;
                    histogram(reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx}),rangeX);
                    histogram(reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx}),rangeX);
                    title('all')
                subplot(1,2,2); hold on;
                    histogram(matched_trial_vals{1},rangeX);
                    histogram(matched_trial_vals{2},rangeX);
                    title('matched');

            matched_trial_vals{1} = eqpop(matched_trial_vals{1},bin_params.reachFeatureBin);
            matched_trial_vals{2} = eqpop(matched_trial_vals{2},bin_params.reachFeatureBin);

            if PID_params.doShuff
                shuffIdxs = zeros(length(matched_trials{early_late}),PID_params.nShuff);
                for shIdx = 1:PID_params.nShuff
                    shuffIdxs(:,shIdx) = randperm(length(matched_trials{early_late}));
                end
            end

            for early_late = 1:2
                for n_lfp = 1:length(params.lfpFeatures)
                    
                    tmp_pairs = [...
                        repmat(find(sig_channels{n_lfp,fidx,1,early_late}),[1 sum(sig_channels{n_lfp,fidx,2,early_late})])' ...
                        reshape(repmat(find(sig_channels{n_lfp,fidx,2,early_late}),[sum(sig_channels{n_lfp,fidx,1,early_late}) 1]), ...
                        [sum(sig_channels{n_lfp,fidx,1,early_late})*sum(sig_channels{n_lfp,fidx,2,early_late}) 1])];
                    tmp_delays = -25:25;
                    tmp_timeStep = 1:floor((length(26:225)-(PID_params.windowSize-PID_params.timeJump))/PID_params.timeJump);
                    
                    tmp_shared = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_unique_m1 = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_unique_dls = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_complementary = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    
                    if PID_params.doShuff
                        tmp_shared_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                        tmp_unique_m1_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                        tmp_unique_dls_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                        tmp_complementary_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                    end
                    
                    parfor pairIdx = 1:size(tmp_pairs,1)
                        
                        m1_chan = tmp_pairs(pairIdx,1);
                        dls_chan = tmp_pairs(pairIdx,2);
                        disp(['Computing: ' params.animals{animal} ' | M1 channel ' num2str(m1_chan) ' & DLS channel ' num2str(dls_chan) ' | ' num2str(early_late) ' of 2 (early/late) | ' params.reachFeatures{fidx}])
                                
                        tmp_M1_LFP = squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1(m1_chan,:,:));
                        tmp_DLS_LFP = squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS(dls_chan,:,:));
                        
                        tmp_M1_LFP = tmp_M1_LFP(:,matched_trials{early_late});
                        tmp_DLS_LFP = tmp_DLS_LFP(:,matched_trials{early_late});
                        
                        pair_shared = zeros(length(tmp_timeStep),length(tmp_delays));
                        pair_unique_m1 = zeros(length(tmp_timeStep),length(tmp_delays));
                        pair_unique_dls = zeros(length(tmp_timeStep),length(tmp_delays));
                        pair_complementary = zeros(length(tmp_timeStep),length(tmp_delays));
                        
                        if PID_params.doShuff
                            pair_shared_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                            pair_unique_m1_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                            pair_unique_dls_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                            pair_complementary_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
                        end
                        
                        tmpBinTimes = [];
                        sTime = 26;
                        
                        for timeStep = 1:length(tmp_timeStep)
                            
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*bin_params.totBin)+tMin)-5001];
                            delayIdx = 1;
                            
                            for delay = tmp_delays
                                
                                eTime = sTime + PID_params.windowSize-1;
                                
                                X1 = tmp_M1_LFP(sTime:eTime,:);
                                X2 = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                
                                Y = matched_trial_vals{early_late};
                                Y = repmat(Y,size(X1,1),1);
                                
                                X1 = X1(:);
                                X2 = X2(:);
                                Y = Y(:);

                                I = PID(Y',X1',X2',PID_params.opts);
                                
                                pair_shared(timeStep,delayIdx) = I.shared;
                                pair_unique_m1(timeStep,delayIdx) = I.uniqueX1;
                                pair_unique_dls(timeStep,delayIdx) = I.uniqueX2;
                                pair_complementary(timeStep,delayIdx) = I.complementary;
                                
                                if PID_params.doShuff
                                    for shidx = 1:PID_params.nShuff
                                        
                                        X1 = tmp_M1_LFP(sTime:eTime,:);
                                        X2 = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                        
                                        Y = matched_trial_vals{early_late};
                                        Y = Y(shuffIdxs(:,shidx));
                                        Y = repmat(Y,size(X1,1),1);
                                        
                                        X1 = X1(:);
                                        X2 = X2(:);
                                        Y = Y(:);
                                        
                                        I = PID(Y',X1',X2',PID_params.opts);
                                        
                                        pair_shared_Sh(timeStep,delayIdx,shidx) = I.shared;
                                        pair_unique_m1_Sh(timeStep,delayIdx,shidx) = I.uniqueX1;
                                        pair_unique_dls_Sh(timeStep,delayIdx,shidx) = I.uniqueX2;
                                        pair_complementary_Sh(timeStep,delayIdx,shidx) = I.complementary;
                                        
                                    end
                                end
                                delayIdx = delayIdx + 1;
                            end
                            sTime = sTime + PID_params.timeJump;
                        end
                        
                        tmp_shared(pairIdx,:,:) = pair_shared;
                        tmp_unique_m1(pairIdx,:,:) = pair_unique_m1;
                        tmp_unique_dls(pairIdx,:,:) = pair_unique_dls;
                        tmp_complementary(pairIdx,:,:) = pair_complementary;
                        
                        if PID_params.doShuff
                            tmp_shared_Sh(pairIdx,:,:,:) = pair_shared_Sh;
                            tmp_unique_m1_Sh(pairIdx,:,:,:) = pair_unique_m1_Sh;
                            tmp_unique_dls_Sh(pairIdx,:,:,:) = pair_unique_dls_Sh;
                            tmp_complementary_Sh(pairIdx,:,:,:) = pair_complementary_Sh;
                        end
                        
                    end
       
                    PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared = tmp_shared;
                    PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1 = tmp_unique_m1;
                    PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dls = tmp_unique_dls;
                    PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementary = tmp_complementary;
                    
                    if PID_params.doShuff
                        PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh = tmp_shared_Sh;
                        PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1Sh = tmp_unique_m1_Sh;
                        PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dlsSh = tmp_unique_dls_Sh;
                        PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementarySh = tmp_complementary_Sh;
                    end
                    
                end
            end
        end
        
        % save results
        if PID_params.savePID
            save([PID_params.save_path 'PID_' params.animals{animal} '_' date ...
                '_LFPearlylate_' num2str(PID_params.doLFPearlylate) ...
                '_LFPday_' num2str(PID_params.doLFPday) ...
                '_Spikesday_' num2str(PID_params.doSpikesday) ...
                '_Pooled.mat'],'-v7.3','PIDout','params','PID_params')
        end

    end
end       
   

end