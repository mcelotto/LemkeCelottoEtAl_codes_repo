function [LFP_early_late, LFP_early_late_noBin, LFP_early_late_noTrialMatch, ...
          reachFeatures_early_late, reachFeatures_early_late_noBin, reachFeatures_early_late_noTrialMatch, reachFeatures_early_late_noBinNoTrialMatch, ...
          spikes, ...
          reachFeatures, reachFeatures_noBin, ...
          tMin, tMax] = rearrange_and_binarize_v2(params,reach_features,bin_params,rf_params)

% intialize output
LFP = cell(1,numel(params.animals));
LFP_early_late = cell(1,numel(params.animals));
spikes = cell(1,numel(params.animals));
spikes_early_late = cell(1,numel(params.animals));
spikes_ave = cell(1,numel(params.animals));
reachFeatures = cell(1,numel(params.animals));
reachFeatures_early_late = cell(1,numel(params.animals));
reachFeatures_early_late_noBin = cell(1,numel(params.animals));

for animal = 1:numel(params.animals)

    % load animal specific neural/behavioral data
    animal_neuralActivity = load([bin_params.path_to_neural_activity_directory '/neuralActivity_' params.animals{animal} '_10ms_04-Jun-2021.mat']);
    
    tMax = animal_neuralActivity.tMax;
    tMin = animal_neuralActivity.tMin;
    animal_reachFeatures = reach_features{animal};
    
    %%% Lowpass/Highpass
    if params.doLowpass || params.doHighpass
        for day = 1:numel(animal_reachFeatures)
            disp(['Bandpassing animal # ', num2str(animal), ' day # ', num2str(day), ' | ' num2str(params.frequencyCutoff), 'Hz']);
            for areaIdx = 1:length(params.areas)
                for lfpIdx = 1:length(params.lfpFeatures)
                    for chan = 1:size(animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),1)
                        if params.doLowpass
                            animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,:) = lowpass(squeeze(animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,:)),params.frequencyCutoff,1000/bin_params.oldBin);
                        else
                            animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,:) = highpass(squeeze(animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,:)),params.frequencyCutoff,1000/bin_params.oldBin);
                        end
                    end
                end
            end
        end
    end
    
    % remove trials based on consecurtive NaN threshold and apply temopral rebinning 
    for day = 1:numel(animal_reachFeatures)
        good_trial_idx = animal_reachFeatures(day).consNaNs<rf_params.nan_thresh;
        disp(['Animal # ' num2str(animal) ' | Trials = ' num2str(sum(good_trial_idx))]);
        for areaIdx = 1:length(params.areas)
            for lfpIdx = 1:length(params.lfpFeatures)
                LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = temporal_rebinning(animal_neuralActivity.LFP{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,good_trial_idx),bin_params.newBin,'mean');
            end
            spikes{animal}{day}.(params.areas{areaIdx}) = temporal_rebinning(animal_neuralActivity.spikeTrains{day}.(params.areas{areaIdx})(:,:,good_trial_idx),bin_params.newBin,'mean');
            spikes_ave{animal}{day}.(params.areas{areaIdx}) = mean(temporal_rebinning(animal_neuralActivity.spikeTrains{day}.(params.areas{areaIdx})(:,:,good_trial_idx),bin_params.newBin,'mean'));
        end
        for rfIdx=1:length(params.reachFeatures)
            reachFeatures{animal}{day}.(params.reachFeatures{rfIdx})= animal_reachFeatures(day).(params.reachFeatures{rfIdx})(good_trial_idx);
        end
        for rfIdx=1:length(params.instReachFeatures)
            reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx})= animal_reachFeatures(day).instant.(params.instReachFeatures{rfIdx})(good_trial_idx,:)';
        end
    end
    
    % rearrange LFP/Spikes into early and late days
    lastDay = numel(animal_reachFeatures);
    for areaIdx = 1:length(params.areas)
        for lfpIdx = 1:length(params.lfpFeatures)
            LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = [];
            LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = [];
            spikes_early_late{animal}{1}.(params.areas{areaIdx}) = [];
            spikes_early_late{animal}{2}.(params.areas{areaIdx}) = [];
            for day = 1:params.num_earlylate_days{animal}
                LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = cat(3,LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}));
                LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = cat(3,LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),LFP{animal}{lastDay-params.num_earlylate_days{animal}+day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}));
                spikes_early_late{animal}{1}.(params.areas{areaIdx}) = cat(3,spikes_early_late{animal}{1}.(params.areas{areaIdx}),mean(spikes{animal}{day}.(params.areas{areaIdx}),1));
                spikes_early_late{animal}{2}.(params.areas{areaIdx}) = cat(3,spikes_early_late{animal}{2}.(params.areas{areaIdx}),mean(spikes{animal}{lastDay-params.num_earlylate_days{animal}+day}.(params.areas{areaIdx}),1));
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

    for rfIdx=1:length(params.instReachFeatures)
        reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx}) = [];
        reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx}) = [];
        for day = 1:params.num_earlylate_days{animal}
            reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx}) = cat(2,reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx}),reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx}));
            reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx}) = cat(2,reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx}),reachFeatures{animal}{lastDay-params.num_earlylate_days{animal}+day}.instant.(params.instReachFeatures{rfIdx}));
        end
    end
    
    % binarize spiking activity, LFP activity, and reach features by day
    for day = 1:numel(animal_reachFeatures)
        for areaIdx = 1:length(params.areas)
            spikes{animal}{day}.(params.areas{areaIdx}) = double(spikes{animal}{day}.(params.areas{areaIdx})>0);
            for trial = 1:size(spikes_ave{animal}{day}.(params.areas{areaIdx}),3)
                spikes_ave{animal}{day}.(params.areas{areaIdx})(1,:,trial) = binr(spikes_ave{animal}{day}.(params.areas{areaIdx})(1,:,trial),2,'eqspace');   
            end
            for lfpIdx = 1:length(params.lfpFeatures)
                for chan = 1:size(LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),1)
                    for trial = 1:size(LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3)
                        LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial) = eqpop(LFP{animal}{day}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial), bin_params.LFPBin);
                    end
                end
            end
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
        % Binarize instantaneous features
        for rfIdx=1:length(params.instReachFeatures)
            for t = 1:size(reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx}),1)
                nanIdx = isnan(reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx})(t,:));
                reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx})(t,:) = eqpop(reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx})(t,:),bin_params.reachFeatureBin);
                reachFeatures{animal}{day}.instant.(params.instReachFeatures{rfIdx})(t,find(nanIdx))=NaN;
            end
        end
    end
    
    % binarize early and late LFP/spiking activity and reach features
    reachFeatures_early_late_noBin{animal}{1} = reachFeatures_early_late{animal}{1};
    reachFeatures_early_late_noBin{animal}{2} = reachFeatures_early_late{animal}{2};
    reachFeatures_early_late_noBinNoTrialMatch{animal} = reachFeatures_early_late{animal}; 

    for early_late = 1:2
        for areaIdx = 1:length(params.areas)
            for trial = 1:size(spikes_early_late{animal}{early_late}.(params.areas{areaIdx}),3)
                spikes_early_late{animal}{early_late}.(params.areas{areaIdx})(:,:,trial) = binr(spikes_early_late{animal}{early_late}.(params.areas{areaIdx})(:,:,trial),2,'eqspace');
            end
            for lfpIdx = 1:length(params.lfpFeatures)
                for chan = 1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),1)
                    for trial = 1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3)
                        LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial) = eqpop(LFP_early_late{animal}{early_late}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(chan,:,trial), bin_params.LFPBin);
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
        % Binarize instantaneous features
        for rfIdx=1:length(params.instReachFeatures)
            for t = 1:size(reachFeatures_early_late{animal}{early_late}.instant.(params.instReachFeatures{rfIdx}),1)
                nanIdx = isnan(reachFeatures_early_late{animal}{early_late}.instant.(params.instReachFeatures{rfIdx})(t,:));
                reachFeatures_early_late{animal}{early_late}.instant.(params.instReachFeatures{rfIdx})(t,:) = eqpop(reachFeatures_early_late{animal}{early_late}.instant.(params.instReachFeatures{rfIdx})(t,:),bin_params.reachFeatureBin);
                reachFeatures_early_late{animal}{early_late}.instant.(params.instReachFeatures{rfIdx})(t,find(nanIdx))=NaN;
            end
        end
    end

    LFP_early_late_noTrialMatch{animal} = LFP_early_late{animal};
    reachFeatures_early_late_noTrialMatch{animal} = reachFeatures_early_late{animal};

    % match trials for early and late LFP activity and reach features    
        for areaIdx = 1:length(params.areas)
            min_trials = min([size(spikes_early_late{animal}{1}.(params.areas{areaIdx}),3) ...
                size(spikes_early_late{animal}{2}.(params.areas{areaIdx}),3)]);
            spikes_early_late{animal}{1}.(params.areas{areaIdx}) = spikes_early_late{animal}{1}.(params.areas{areaIdx})(:,:,1:min_trials);
            spikes_early_late{animal}{2}.(params.areas{areaIdx}) = spikes_early_late{animal}{2}.(params.areas{areaIdx})(:,:,end-(min_trials-1):end);
            for lfpIdx = 1:length(params.lfpFeatures) 
                min_trials = min([size(LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3) ...
                                  size(LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}),3)]);
                LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = LFP_early_late{animal}{1}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,1:min_trials);
                LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx}) = LFP_early_late{animal}{2}.(params.lfpFeatures{lfpIdx}).(params.areas{areaIdx})(:,:,end-(min_trials-1):end);
            end
        end
        for rfIdx=1:length(params.reachFeatures)
            min_trials = min([numel(reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})) ...
                              numel(reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}))]);         
            reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})(1:min_trials);
            reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx})(end-(min_trials-1):end);     
            %                 reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx})(1:min_trials);
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
        % Match instant features
        for rfIdx=1:length(params.instReachFeatures)
            min_trials = min([size(reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx}),2) ...
                              size(reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx}),2)]);         
            reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx}) = reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx})(:,1:min_trials);
            reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx}) = reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx})(:,end-(min_trials-1):end);     
            %                 reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}) = reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx})(1:min_trials);
            
            diffNanTrials = sum(~isnan(reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx})(1,:))) - sum(~isnan(reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx})(1,:)));
            if diffNanTrials > 0 % we have to remove trials from early
                nanId = ~isnan(reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx}));
                remIdxs = find(nanId == 1);
                remIdxs = remIdxs(1:diffNanTrials);
                reachFeatures_early_late{animal}{1}.instant.(params.instReachFeatures{rfIdx})(remIdxs) = nan;
            elseif diffNanTrials < 0 % we have to remove trials from late
                nanId = ~isnan(reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx}));
                remIdxs = find(nanId == 1);
                remIdxs = remIdxs(1:-diffNanTrials);
                reachFeatures_early_late{animal}{2}.instant.(params.instReachFeatures{rfIdx})(remIdxs) = nan;
            end
            
            disp(['Animal # ' num2str(animal) ' | reach feature # ' num2str(rfIdx) ...
                ' | good trials early: ' num2str(sum(~isnan(reachFeatures_early_late{animal}{1}.(params.reachFeatures{rfIdx})))) ...
                ' & late: ' num2str(sum(~isnan(reachFeatures_early_late{animal}{2}.(params.reachFeatures{rfIdx}))))])

    
    end
                
end

end
   