function [MIout] = compute_mutualInformation_behaviorMatch(LFP_early_late_noTrialMatch, spikes, reachFeatures_early_late, reachFeatures_early_late_noBinNoTrialMatch, reachFeatures, reachFeatures_noBin, params, MI_params, PID_params, tMin, binSize, save_params, features_to_compute, bin_params)

 parfor animal = 1:numel(params.animals)
   
    % calculate LFP MI by early/late
    if MI_params.doLFP==1
        for early_late = 1:2
            for n_lfp = 1:length(params.lfpFeatures)
                for aidx = 1:length(params.areas)                   
                    for fidx = features_to_compute

                        tmp_min = prctile([reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx}) reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx})],2);
                        tmp_max = prctile([reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx}) reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx})],98);
                        tmp_range = tmp_max-tmp_min;
                        dist_thresh = tmp_range/8;

%                         if fidx==9
%                             dist_thresh = .5;
%                         elseif fidx==13 
%                             dist_thresh = 5;
%                         elseif fidx==14
%                             dist_thresh = .2;
%                         end

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
                    
                         %%% match early TO late
%                         all_early_vals = reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx});
%                         matched_trial_vals = cell(1,2);
%                         matched_trials = cell(1,2);
%                         for i = 1:length(reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx}))
%                             tmp_late_val = reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx})(i);
%                             [val, idx] = min(abs(all_early_vals-tmp_late_val));
%                             if val<dist_thresh
%                                 matched_trial_vals{2} = [matched_trial_vals{2} tmp_late_val];
%                                 matched_trial_vals{1} = [matched_trial_vals{1} all_early_vals(idx)];
%                                 matched_trials{2} = [matched_trials{2} i];
%                                 matched_trials{1} = [matched_trials{1} idx];
%                             end
%                             all_early_vals(idx) = NaN;
%                         end

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

                        if MI_params.doShuff
                            shuffIdxs = zeros(length(matched_trials{early_late}),MI_params.nShuff);
                            for shIdx = 1:MI_params.nShuff
                                shuffIdxs(:,shIdx) = randperm(length(matched_trials{early_late}));
                            end
                        end

                        for chan = 1:size(LFP_early_late_noTrialMatch{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | ', params.lfpFeatures{n_lfp} ,' channel ', num2str(chan), ' | ' num2str(early_late) ' of 2 (early/late) | reach feature: ', params.reachFeatures{fidx}])

                            LFPsignal = squeeze(LFP_early_late_noTrialMatch{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:));
                            LFPsignal = LFPsignal(:,matched_trials{early_late});
                            
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(LFPsignal,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(LFPsignal,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = LFPsignal(sTime:eTime,:);
                                    Y = matched_trial_vals{early_late};
                                    Y = repmat(Y,size(X,1),1);
                                    
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,timeStepCount) = I1{1}(1);
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];

                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            
                                            Y = matched_trial_vals{early_late};
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);                                            
                                            
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,timeStepCount,shIdx) = I1{1}(1);
                                        end
                                    end
                                    
                                    sTime = sTime + MI_params.timeJump;
                                    timeStepCount = timeStepCount+1;
                                end
                            end
                            
                        end
                    end
                end
            end
        end
    end

    % calculate spikes MI by day
    if MI_params.doSpikes==1
        for day = 1:numel(spikes{animal})
            for aidx = 1:length(params.areas)
                for fidx = features_to_compute

                    if fidx==9
                        dist_thresh = .5;
                    elseif fidx==13
                        dist_thresh = 5;
                    elseif fidx==14
                        dist_thresh = .2;
                    end


                  %%% match late TO early                    
                    matched_trials = cell(1,numel(spikes{animal}));
                    matched_trials_idx = cell(1,numel(spikes{animal}));
                    for i = 1:params.num_earlylate_days{animal}
                        late_vals = reachFeatures_noBin{animal}{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i}.(params.reachFeatures{fidx});
                        for i2 = 1:length(reachFeatures_noBin{animal}{i}.(params.reachFeatures{fidx}))
                            tmp_early = reachFeatures_noBin{animal}{i}.(params.reachFeatures{fidx})(i2);
                            [val, idx] = min(abs(late_vals-tmp_early));
                            if val<dist_thresh
                                matched_trials{i} = [matched_trials{i} tmp_early];
                                matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} = [matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} late_vals(idx)];
                                matched_trials_idx{i} = [matched_trials_idx{i} i2];
                                matched_trials_idx{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} = [matched_trials_idx{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} idx];
                            end
                            late_vals(idx) = NaN;
                        end
                    end
                        
%                     %%% match early TO late
%                     matched_trials = cell(1,numel(spikes{animal}));
%                     matched_trials_idx = cell(1,numel(spikes{animal}));
%                     for i = 1:params.num_earlylate_days{animal}
%                         early_vals = reachFeatures_noBin{animal}{i}.(params.reachFeatures{fidx});
%                         for i2 = 1:length(reachFeatures_noBin{animal}{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i}.(params.reachFeatures{fidx}))
%                             tmp_late = reachFeatures_noBin{animal}{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i}.(params.reachFeatures{fidx})(i2);
%                             [val, idx] = min(abs(early_vals-tmp_late));
%                             if val<dist_thresh
%                                 matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} = [matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} tmp_late];
%                                 matched_trials{i} = [matched_trials{i} early_vals(idx)];
%                                 matched_trials_idx{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} = [matched_trials_idx{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i} i2];
%                                 matched_trials_idx{i} = [matched_trials_idx{i} idx];
%                             end
%                             early_vals(idx) = NaN;
%                         end
%                     end

                    disp(['ANIMAL ' num2str(animal) ' | ' num2str(cellfun(@length,matched_trials))])

                    if animal==2 && fidx==9
                        matched_trials{3} = matched_trials{3}(1:21);
                        matched_trials_idx{3} = matched_trials_idx{3}(1:21);
                    elseif animal==3 && fidx==9
                        matched_trials{3} = matched_trials{3}(1:37);                    
                        matched_trials_idx{3} = matched_trials_idx{3}(1:37);                    
                    end

                    if animal==2 && fidx==13
                        matched_trials{3} = matched_trials{3}(1:12);
                        matched_trials_idx{3} = matched_trials_idx{3}(1:12);
                    elseif animal==3 && fidx==13
                        matched_trials{3} = matched_trials{3}(1:12);
                        matched_trials_idx{3} = matched_trials_idx{3}(1:12);
                    end

                    if animal==2 && fidx==14
                        matched_trials{3} = matched_trials{3}(1:24);
                        matched_trials_idx{3} = matched_trials_idx{3}(1:24);
                    elseif animal==3 && fidx==14
                        matched_trials{3} = matched_trials{3}(1:37);                    
                        matched_trials_idx{3} = matched_trials_idx{3}(1:37);                    
                    end

                    tmp_min = [];
                    tmp_max = [];                    
                    for i = 1:params.num_earlylate_days{animal}
                        tmp_min = [tmp_min min(matched_trials{i}) min(matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i})];
                        tmp_max = [tmp_max max(matched_trials{i}) max(matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i})];
                    end
                    tmp_min = min(tmp_min);
                    tmp_max = min(tmp_max);

                    rangeX = [tmp_min:(tmp_max-tmp_min)/20:tmp_max];
                    figure;
                    for i = 1:params.num_earlylate_days{animal}
                        subplot(params.num_earlylate_days{animal},2,1+2*(i-1)); hold on;
                        histogram(reachFeatures_noBin{animal}{i}.(params.reachFeatures{fidx}),rangeX);
                        histogram(reachFeatures_noBin{animal}{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i}.(params.reachFeatures{fidx}),rangeX);
                        subplot(params.num_earlylate_days{animal},2,i*2); hold on;
                        histogram(matched_trials{i},rangeX);
                        histogram(matched_trials{length(reachFeatures_noBin{animal})-params.num_earlylate_days{animal}+i},rangeX);
                    end

                    if length(matched_trials{day})>10

                        for i = 1:length(matched_trials)
                            if ~isempty(matched_trials{day})
                                matched_trials{day} = eqpop(matched_trials{day},bin_params.reachFeatureBin);
                            end
                        end

                        if MI_params.doShuff
                            shuffIdxs = zeros(length(matched_trials{day}),MI_params.nShuff);
                            for shIdx = 1:MI_params.nShuff
                                shuffIdxs(:,shIdx) = randperm(length(matched_trials{day}));
                            end
                        end
                    
                        for unit = 1:size(spikes{animal}{day}.(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | unit ', num2str(unit), ' | day ' num2str(day) ' of ' num2str(numel(spikes{1,1})) ' | reach feature: ', params.reachFeatures{fidx}])
                            
                            spikeActivity = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:));
                            spikeActivity = spikeActivity(:,matched_trials_idx{day});
    
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(spikeActivity,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(spikeActivity,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = spikeActivity(sTime:eTime,:);
                                    Y = matched_trials{day};
                                    Y = repmat(Y,size(X,1),1);
                                    %%% SL 8.31.2022
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,timeStepCount) = I1{1}(1);
                                    
                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            
                                            Y = matched_trials{day};
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);
                                            
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,timeStepCount,shIdx) = I1{1}(1);
                                        end
                                    end
                                    
                                    sTime = sTime + MI_params.timeJump;
                                    timeStepCount = timeStepCount+1;
                                end
                            end
                            
                        end

                    end
                end
            end
        end
        
    end
    
end

% save results
if save_params.save
    save([save_params.save_path '\MI_behMatch_' date ...
        '_LFPearlylate_' num2str(MI_params.doLFP) ...
        '_Spikes_' num2str(MI_params.doSpikes) ...
        '.mat'],'-v7.3','MIout')
end

end
