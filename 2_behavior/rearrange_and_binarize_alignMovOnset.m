function [LFP_early_late, LFP_early_late_noBin, LFP_early_late_noTrialMatch, ...
          reachFeatures_early_late, reachFeatures_early_late_noBin, reachFeatures_early_late_noTrialMatch, reachFeatures_early_late_noBinNoTrialMatch, ...
          spikes, ...
          reachFeatures, reachFeatures_noBin, ...
          tMin, tMax] = rearrange_and_binarize_alignMovOnset(params,reach_features,bin_params,paths)

% intialize output
LFP = cell(1,numel(params.animals));
LFP_early_late = cell(1,numel(params.animals));
LFP_early_late_noBin = cell(1,numel(params.animals));
spikes = cell(1,numel(params.animals));
reachFeatures = cell(1,numel(params.animals));
reachFeatures_noBin = cell(1,numel(params.animals));
reachFeatures_early_late = cell(1,numel(params.animals));
reachFeatures_early_late_noBin = cell(1,numel(params.animals));

for animal = 1:numel(params.animals)

    % load animal specific neural/behavioral data
    animal_neuralActivity = load([paths.neural_activity '\neuralActivity_' params.animals{animal} '_10ms_04-Jun-2021.mat']);
    animal_neuralActivity_movOnset = animal_neuralActivity;

    for day = 1:length(animal_neuralActivity.LFP)
        for trial = 1:size(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1,3)
            if ~isnan(reach_features{animal}(day).movDur(trial))
                animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.M1(:,:,trial) = ...
                      [nan(size(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1,1),(reach_features{animal}(day).movDur(trial)/10)) ...
                      animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1(:,1:end-(reach_features{animal}(day).movDur(trial)/10),trial) ...
                      ];
                animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.DLS(:,:,trial) = ...
                      [nan(size(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.DLS,1),(reach_features{animal}(day).movDur(trial)/10)) ...
                      animal_neuralActivity.LFP{day}.reach_related_lfp_ref.DLS(:,1:end-(reach_features{animal}(day).movDur(trial)/10),trial) ...
                      ];                
            else
                animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.M1(:,:,trial) = nan(size(animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.M1(:,:,trial)));
                animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.DLS(:,:,trial) = nan(size(animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.DLS(:,:,trial)));
            end
        end
    end

    for day = 1:length(animal_neuralActivity.spikeTrains)
        for trial = 1:size(animal_neuralActivity.spikeTrains{day}.M1,3)
            if ~isnan(reach_features{animal}(day).movDur(trial))
                animal_neuralActivity_movOnset.spikeTrains{day}.M1(:,:,trial) = ...
                      [nan(size(animal_neuralActivity.spikeTrains{day}.M1,1),(reach_features{animal}(day).movDur(trial)/10)) ...
                      animal_neuralActivity.spikeTrains{day}.M1(:,1:end-(reach_features{animal}(day).movDur(trial)/10),trial) ...
                      ];
                animal_neuralActivity_movOnset.spikeTrains{day}.DLS(:,:,trial) = ...
                      [nan(size(animal_neuralActivity.spikeTrains{day}.DLS,1),(reach_features{animal}(day).movDur(trial)/10)) ...
                      animal_neuralActivity.spikeTrains{day}.DLS(:,1:end-(reach_features{animal}(day).movDur(trial)/10),trial) ...
                      ];                
            else
                animal_neuralActivity_movOnset.spikeTrains{day}.M1(:,:,trial) = nan(size(animal_neuralActivity_movOnset.spikeTrains{day}.M1(:,:,trial)));
                animal_neuralActivity_movOnset.spikeTrains{day}.DLS(:,:,trial) = nan(size(animal_neuralActivity_movOnset.spikeTrains{day}.DLS(:,:,trial)));
            end
        end
    end

    %%% PLOT LFP

    all_peth = [];
    mo_peth = [];
    for day = 1:length(animal_neuralActivity.LFP)
        for trial = 1:size(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1,3)
            all_peth = [all_peth; nanmean(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1(:,:,trial))];
            mo_peth = [mo_peth; nanmean(animal_neuralActivity_movOnset.LFP{day}.reach_related_lfp_ref.M1(:,:,trial))];
        end
    end

    for n = 1:size(all_peth,1)
        tmp_mean = mean(all_peth(n,~isnan(all_peth(n,:))));
        tmp_std = std(all_peth(n,~isnan(all_peth(n,:))));
        all_peth(n,:) = (all_peth(n,:)-tmp_mean)/tmp_std;
    end

    for n = 1:size(mo_peth,1)
        tmp_mean = mean(mo_peth(n,~isnan(mo_peth(n,:))));
        tmp_std = std(mo_peth(n,~isnan(mo_peth(n,:))));
        mo_peth(n,:) = (mo_peth(n,:)-tmp_mean)/tmp_std;
    end

    figure; 
        subplot(2,2,1);
            imagesc(all_peth)       
        subplot(2,2,3);
            plot(nanmean(all_peth))
        subplot(2,2,2);
            imagesc(mo_peth)       
        subplot(2,2,4);
            plot(nanmean(mo_peth))
        sgtitle(['animal ' num2str(animal)])

    %%% PLOT SPIKES

    all_peth = [];
    mo_peth = [];
    for day = 1:length(animal_neuralActivity.spikeTrains)
        for trial = 1:size(animal_neuralActivity.spikeTrains{day}.M1,3)
            all_peth = [all_peth; nanmean(animal_neuralActivity.spikeTrains{day}.M1(:,:,trial))];
            mo_peth = [mo_peth; nanmean(animal_neuralActivity_movOnset.spikeTrains{day}.M1(:,:,trial))];
        end
    end

    for n = 1:size(all_peth,1)
        tmp_mean = mean(all_peth(n,~isnan(all_peth(n,:))));
        tmp_std = std(all_peth(n,~isnan(all_peth(n,:))));
        all_peth(n,:) = (all_peth(n,:)-tmp_mean)/tmp_std;
    end

    for n = 1:size(mo_peth,1)
        tmp_mean = mean(mo_peth(n,~isnan(mo_peth(n,:))));
        tmp_std = std(mo_peth(n,~isnan(mo_peth(n,:))));
        mo_peth(n,:) = (mo_peth(n,:)-tmp_mean)/tmp_std;
    end

    figure; 
        subplot(2,2,1);
            imagesc(all_peth)       
        subplot(2,2,3);
            plot(nanmean(all_peth))
        subplot(2,2,2);
            imagesc(mo_peth)       
        subplot(2,2,4);
            plot(nanmean(mo_peth))
        sgtitle(['animal ' num2str(animal)])
    
    %%%

    animal_neuralActivity = animal_neuralActivity_movOnset;

    tMax = animal_neuralActivity.tMax;
    tMin = animal_neuralActivity.tMin;
    animal_reachFeatures = reach_features{animal};
   
    % remove trials based on consecurtive NaN threshold and apply temopral rebinning 
    for day = 1:numel(animal_reachFeatures)
        good_trial_idx = animal_reachFeatures(day).consNaNs<bin_params.nan_thresh;
        disp(['Animal # ' num2str(animal) ' | Trials = ' num2str(sum(good_trial_idx))]);
        for areaIdx = 1:length(params.areas)
            for lfpIdx = 1:length(params.lfpFeatures)
                LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = temporal_rebinning(animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,good_trial_idx),bin_params.newBin,'mean');
            end
            spikes{animal}{day}.(params.areas{areaIdx}) = temporal_rebinning(animal_neuralActivity.spikeTrains{day}.(params.areas{areaIdx})(:,:,good_trial_idx),bin_params.newBin,'mean');
        end
        for rfIdx=1:length(params.reachFeatures)
            reachFeatures{animal}{day}.(params.reachFeatures{rfIdx})= animal_reachFeatures(day).(params.reachFeatures{rfIdx})(good_trial_idx);
        end
    end
    
    % rearrange LFP into early and late days
    lastDay = numel(animal_reachFeatures);
    for areaIdx = 1:length(params.areas)
        for lfpIdx = 1:length(params.lfpFeatures)
            LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = [];
            LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = [];
            for day = 1:params.num_earlylate_days{animal}
                LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = cat(3,LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}));
                LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = cat(3,LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),LFP{animal}{lastDay-params.num_earlylate_days{animal}+day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}));
            end             
        end
    end
    
    % rearrange reach features into early and late days
    for rfIdx=1:length(params.reachFeatures)
        reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}) = [];
        reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}) = [];
        for day = 1:params.num_earlylate_days{animal}
            reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}) = cat(2,reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}),reachFeatures{animal}{day}.(params.reachFeatures{rfIdx}));
            reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}) = cat(2,reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}),reachFeatures{animal}{lastDay-params.num_earlylate_days{animal}+day}.(params.reachFeatures{rfIdx}));
        end
    end

    LFP_early_late_noBin{animal} = LFP_early_late{animal};
    reachFeatures_early_late_noBin{animal} = reachFeatures_early_late{animal};
    reachFeatures_early_late_noBinNoTrialMatch{animal} = reachFeatures_early_late{animal};    
    reachFeatures_noBin{animal} = reachFeatures{animal};

    % binarize spiking activity and reach features by day
    for day = 1:numel(animal_reachFeatures)
        for areaIdx = 1:length(params.areas)
            isnan_mat = isnan(spikes{animal}{day}.(params.areas{areaIdx}));
            tmp = double(spikes{animal}{day}.(params.areas{areaIdx})>0);
            tmp(isnan_mat) = NaN;
            spikes{animal}{day}.(params.areas{areaIdx}) = tmp;
        end
        for rfIdx=1:length(params.reachFeatures)
            if strcmp(params.reachFeatures{rfIdx},'success_rate')
                nanIdx = isnan(reachFeatures{animal}{day}.(params.reachFeatures{rfIdx}));
                reachFeatures{animal}{day}.(params.reachFeatures{rfIdx})(find(nanIdx))=NaN;
            else
                nanIdx = isnan(reachFeatures{animal}{day}.(params.reachFeatures{rfIdx}));
                reachFeatures{animal}{day}.(params.reachFeatures{rfIdx}) = eqpop(reachFeatures{animal}{day}.(params.reachFeatures{rfIdx}),bin_params.reachFeatureBin);
                reachFeatures{animal}{day}.(params.reachFeatures{rfIdx})(find(nanIdx))=NaN;
            end
        end
    end
    
    % binarize early and late LFP and reach features
    for early_late = 1:2
        for areaIdx = 1:length(params.areas)
            for lfpIdx = 1:length(params.lfpFeatures)
                for chan = 1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),1)
                    for trial = 1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3)
                        isnan_tmp = isnan(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial));
                        tmp = eqpop(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial), bin_params.LFPBin);
                        tmp(isnan_tmp) = NaN;
                        LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial) = tmp;
                    end
                end
            end
        end
        for rfIdx=1:length(params.reachFeatures)
            if strcmp(params.reachFeatures{rfIdx},'success_rate')
                nanIdx = isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{rfIdx}));
                reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{rfIdx})(find(nanIdx))=NaN;
            else
                nanIdx = isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{rfIdx}));
                reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{rfIdx}) = eqpop(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{rfIdx}),bin_params.reachFeatureBin);
                reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{rfIdx})(find(nanIdx))=NaN;
            end
        end
    end

    LFP_early_late_noTrialMatch{animal} = LFP_early_late{animal};
    reachFeatures_early_late_noTrialMatch{animal} = reachFeatures_early_late{animal};

    % match trials for early and late LFP activity and reach features
    for areaIdx = 1:length(params.areas)
        for lfpIdx = 1:length(params.lfpFeatures)
            min_trials = min([size(LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3) ...
                size(LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3)]);
            LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,1:min_trials);
            LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,end-(min_trials-1):end);
            LFP_early_late_noBin{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = LFP_early_late_noBin{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,1:min_trials);
            LFP_early_late_noBin{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = LFP_early_late_noBin{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,end-(min_trials-1):end);
        end
    end

    for rfIdx=1:length(params.reachFeatures)
        min_trials = min([numel(reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})) ...
            numel(reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}))]);
        reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})(1:min_trials);
        reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx})(end-(min_trials-1):end);
        reachFeatures_early_late_noBin{animal}{1}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late_noBin{animal}{1}.(params.reachFeatures{rfIdx})(1:min_trials);
        reachFeatures_early_late_noBin{animal}{2}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late_noBin{animal}{2}.(params.reachFeatures{rfIdx})(end-(min_trials-1):end);
        diffNanTrials = sum(~isnan(reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}))) - sum(~isnan(reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx})));
        if diffNanTrials > 0 % we have to remove trials from early
            nanId = ~isnan(reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}));
            remIdxs = find(nanId == 1);
            remIdxs = remIdxs(1:diffNanTrials);
            reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})(remIdxs) = nan;
            reachFeatures_early_late_noBin{animal}{1}.(params.reachFeatures{rfIdx})(remIdxs) = nan;
        elseif diffNanTrials < 0 % we have to remove trials from late
            nanId = ~isnan(reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}));
            remIdxs = find(nanId == 1);
            remIdxs = remIdxs(1:-diffNanTrials);
            reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx})(remIdxs) = nan;
            reachFeatures_early_late_noBin{animal}{2}.(params.reachFeatures{rfIdx})(remIdxs) = nan;
        end

        disp(['Animal # ' num2str(animal) ' | reach feature # ' num2str(rfIdx) ...
            ' | good trials early: ' num2str(sum(~isnan(reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})))) ...
            ' & late: ' num2str(sum(~isnan(reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}))))])

    end

                
end

end
   