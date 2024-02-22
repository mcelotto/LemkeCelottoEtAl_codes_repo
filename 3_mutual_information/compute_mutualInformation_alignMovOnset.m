function [MIout] = compute_mutualInformation_alignMovOnset(LFP_early_late, spikes, reachFeatures_early_late, reachFeatures, params, MI_params, PID_params, tMin, binSize, trial_match,filtered,save_params,features_to_compute)

parfor animal = 1:numel(params.animals)
    
    % calculate LFP MI by early/late
    if MI_params.doLFP==1
        for early_late = 1:2
            for n_lfp = 1:length(params.lfpFeatures)
                for aidx = 1:length(params.areas)
                    
                    if MI_params.doShuff
                        shuffIdxs = zeros(size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),3),MI_params.nShuff);
                        for shIdx = 1:MI_params.nShuff
                            shuffIdxs(:,shIdx) = randperm(size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),3));
                        end
                    end
                    
                    for fidx = features_to_compute%1:length(params.reachFeatures)
                        for chan = 1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | ', params.lfpFeatures{n_lfp} ,' channel ', num2str(chan), ' | ' num2str(early_late) ' of 2 (early/late) | reach feature: ', params.reachFeatures{fidx}])
                            
                            LFPsignal = squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:));
                            nonNaNtrials = ~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})) & sum(isnan(LFPsignal(75:end,:)))==0;
                            LFPsignal = LFPsignal(:,nonNaNtrials);
                                  
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(LFPsignal,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(LFPsignal,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = LFPsignal(sTime:eTime,:);
                                    Y = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                    Y = repmat(Y,size(X,1),1);
                                    
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    if sum(isnan(X))>0
                                        MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,timeStepCount) = NaN;
                                        for shIdx = 1:MI_params.nShuff
                                            MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,timeStepCount,shIdx) = NaN;
                                        end
                                    else
                                        I1 = information(X',Y',MI_params.opts,{'I'});
                                        MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,timeStepCount) = I1{1}(1);
                                        if MI_params.doShuff
                                            for shIdx = 1:MI_params.nShuff
                                                Y = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                                tmpShuffle=shuffIdxs(:,shIdx);
                                                tmpShuffle=tmpShuffle(shuffIdxs(:,shIdx)<=sum(nonNaNtrials));
                                                Y = Y(tmpShuffle);
                                                Y = repmat(Y,MI_params.windowSize,1);
                                                Y = Y(:);                                            
                                                I1 = information(X',Y',MI_params.opts,{'I'});
                                                MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,timeStepCount,shIdx) = I1{1}(1);
                                            end
                                        end
                                    end
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];
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
                
                if MI_params.doShuff
                    shuffIdxs = zeros(size(spikes{animal}{day}.(params.areas{aidx}),3),MI_params.nShuff);
                    for shIdx = 1:MI_params.nShuff
                        shuffIdxs(:,shIdx) = randperm(size(spikes{animal}{day}.(params.areas{aidx}),3));
                    end
                end
                
                for fidx = features_to_compute
                    for unit = 1:size(spikes{animal}{day}.(params.areas{aidx}),1)
                        
                        disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | unit ', num2str(unit), ' | day ' num2str(day) ' of ' num2str(numel(spikes{1,1})) ' | reach feature: ', params.reachFeatures{fidx}])
                        
                        spikeActivity = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:));
                        nonNaNtrials = ~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})) & sum(isnan(spikeActivity(75:end,:)))==0;
                        spikeActivity = spikeActivity(:,nonNaNtrials);
                        
%                         tmpBinTimes = [];
                        newBinTimes = [];
%                         newBinTimes_totSpikes = [];
                        timeStepCount = 1;
                        sTime = 1;
                        for timeStep = 1:size(spikeActivity,1)/MI_params.timeJump
                            eTime = sTime + MI_params.windowSize - 1;
                            if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(spikeActivity,1)
                                sTime = sTime + MI_params.timeJump;
                            else
                                X = spikeActivity(sTime:eTime,:);

                                Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
                                Y = repmat(Y,size(X,1),1);
                                
                                X = X(:);
                                Y = Y(:);

                                if sum(isnan(X))>0
                                    MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,timeStepCount) = NaN;
                                    for shIdx = 1:MI_params.nShuff
                                        MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,timeStepCount,shIdx) = NaN;
                                    end
                                else
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,timeStepCount) = I1{1}(1);
                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            tmpShuffle=tmpShuffle(shuffIdxs(:,shIdx)<=sum(nonNaNtrials));
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,timeStepCount,shIdx) = I1{1}(1);
                                        end
                                    end
                                end

                                %%% SL 8.31.2022
                                newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];
%                                 newBinTimes_totSpikes = [newBinTimes_totSpikes; sum(spikeActivity(sTime:eTime,:),'all')];
                                %%%
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

    % save results
    if save_params.save
%         if filtered==0
            save([save_params.save_path 'MI_' date ...
                '_LFPearlylate_' num2str(MI_params.doLFP) ...
                '_Spikes_' num2str(MI_params.doSpikes) ...
                '_MatchTrials_' num2str(trial_match) ...
                '_movOnsetAligned.mat'],'-v7.3','MIout')
%         elseif filtered==1
%             save([MI_params.save_path 'MI_' date ...
%                 '_LFPearlylate_' num2str(MI_params.doLFPearlylate) ...
%                 '_LFPday_' num2str(MI_params.doLFPday) ...
%                 '_Spikesearlylate_' num2str(MI_params.doSpikesearlylate) ...
%                 '_Spikesday_' num2str(MI_params.doSpikesday) ...
%                 '_SpikesdayAve_' num2str(MI_params.doSpikesdayAve) ...
%                 '_MatchTrials_' num2str(trial_match) ...
%                 '_movOnsetAligned.mat'],'-v7.3','MIout','newBinTimes')
%         end
    end

end
