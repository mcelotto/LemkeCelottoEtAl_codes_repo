function [MIout] = compute_instMutualInformation(LFP, spikes, spikes_ave, reachFeatures, LFP_early_late, spikes_early_late, reachFeatures_early_late, features_to_compute, params, MI_params, PID_params, rf_params, tMin, binSize)

% select trajectories and neural activity time windows to compute MI 
rf_TW = [-1.49:0.01:0.5];
neural_TW = [-1.49:0.01:1];
sel_rf_TW(1) = find(rf_TW == MI_params.instMIinterval(1))+1;
sel_rf_TW(2) = find(rf_TW == MI_params.instMIinterval(2));
sel_neural_TW(1) = find(neural_TW == MI_params.instMIinterval(1))+1;
sel_neural_TW(2) = find(neural_TW == MI_params.instMIinterval(2));

rf_idxs = find(ismember(params.instReachFeatures, features_to_compute));

parfor animal = 1:numel(params.animals) % parallel loop over animals
    
    % calculate LFP MI by early/late
    if MI_params.doLFPearlylate==1
        for early_late = 1:2
            for n_lfp = 1:length(params.lfpFeatures)
                for aidx = 1:length(params.areas)
                    
                    if MI_params.doShuff
                        shuffIdxs = zeros(size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),3),MI_params.nShuff);
                        for shIdx = 1:MI_params.nShuff
                            shuffIdxs(:,shIdx) = randperm(size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),3));
                        end
                    end
                    
                    for fidx = rf_idxs
                        selReachFeature = reachFeatures_early_late{animal}{early_late}.instant.(params.instReachFeatures{fidx})(sel_rf_TW(1):sel_rf_TW(2),:);
                        
                        for chan = 1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | ', params.lfpFeatures{n_lfp} ,' channel ', num2str(chan), ' | ' num2str(early_late) ' of 2 (early/late) | reach feature: ', params.instReachFeatures{fidx}])
                            
                            LFPsignal = squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:));
                            nonNaNtrials = find(~isnan(sum(selReachFeature,1)));
                            LFPsignal = LFPsignal(:,nonNaNtrials);
                            LFPsignal = LFPsignal(sel_neural_TW(1):sel_neural_TW(2),:); % take signal in time window of interest
                            
                            tmpBinTimes = [];
                            timeStepCount1 = 1;
                            sTime1 = 1;
                            for timeStep1 = 1:size(LFPsignal,1)/MI_params.timeJump-2 % loop over neural activity time
                                eTime1 = sTime1 + MI_params.windowSize - 1;
                                X = LFPsignal(sTime1:eTime1,:);
                                X = X(:);
                                tmpBinTimes = [tmpBinTimes (((sTime1-1)*binSize)+tMin)-5001]; % not creating timeBins1 and timeBins2 since they would be the same
                                
                                timeStepCount2 = 1;
                                sTime2 = 1;
                                for timeStep2 = 1:size(selReachFeature,1)/MI_params.timeJump-2
                                    eTime2 = sTime2 + MI_params.windowSize - 1;
                                    
                                    if timeStep2 >= 1 % only compute info when feature time is after neural time
                                        Y = selReachFeature(sTime2:eTime2,nonNaNtrials);

    %                                     Y = repmat(Y,size(X,1),1);
                                        Y = Y(:);

                                        I1 = information(X',Y',MI_params.opts,{'I'});
                                        MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.instReachFeatures{fidx}).info(chan,timeStepCount1,timeStepCount2) = I1{1}(1);

                                        if MI_params.doShuff
                                            for shIdx = 1:MI_params.nShuff

                                                Y = selReachFeature(sTime2:eTime2,nonNaNtrials);
                                                tmpShuffle=shuffIdxs(:,shIdx);
                                                tmpShuffle=tmpShuffle(shuffIdxs(:,shIdx)<=length(nonNaNtrials));
                                                Y = Y(:,tmpShuffle);
    %                                             Y = repmat(Y,MI_params.windowSize,1);
                                                Y = Y(:);                                            

                                                I1 = information(X',Y',MI_params.opts,{'I'});
                                                MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.instReachFeatures{fidx}).infoSh(chan,timeStepCount1,timeStepCount2,shIdx) = I1{1}(1);
                                            end

                                        end
                                    end

                                    MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.instReachFeatures{fidx}).binTimes = tmpBinTimes;
                                    sTime2 = sTime2 + MI_params.timeJump;
                                    timeStepCount2 = timeStepCount2+1;
                                end
                                
                                sTime1 = sTime1 + MI_params.timeJump;
                                timeStepCount1 = timeStepCount1+1;
                            end
                            
                        end
                    end
                end
            end
        end
    end
    
    % calculate spikes MI by early/late
    if MI_params.doSpikesearlylate==1
        for early_late = 1:2
            for aidx = 1:length(params.areas)
                
                if MI_params.doShuff
                    shuffIdxs = zeros(size(spikes_early_late{animal}{early_late}.(params.areas{aidx}),3),MI_params.nShuff);
                    for shIdx = 1:MI_params.nShuff
                        shuffIdxs(:,shIdx) = randperm(size(spikes_early_late{animal}{early_late}.(params.areas{aidx}),3));
                    end
                end

                for fidx = rf_idxs
                    for unit = 1
                        
                        disp(['Computing... animal: ' params.animals{animal} ' | area: ' params.areas{aidx} ' | ' num2str(early_late) ' of 2 (early/late) | reach feature: ' params.instReachFeatures{fidx}])

                        Spikesignal = squeeze(spikes_early_late{animal}{early_late}.(params.areas{aidx})(unit,:,:));
                        nonNaNtrials = find(~isnan(sum(selReachFeature,1)));
                        Spikesignal = Spikesignal(sel_neural_TW(1):sel_neural_TW(2),nonNaNtrials);

                        tmpBinTimes = [];
                        timeStepCount1 = 1;
                        sTime1 = 1;
                        for timeStep1 = 1:size(Spikesignal,1)/MI_params.timeJump-2
                            eTime1 = sTime1 + MI_params.windowSize - 1;
                            
                            X = Spikesignal(sTime1:eTime1,:);
                            X = X(:);
                            tmpBinTimes = [tmpBinTimes (((sTime1-1)*binSize)+tMin)-5001];
                            
                            timeStepCount2 = 1;
                            sTime2 = 1;
                            for timeStep2 = 1:size(selReachFeature,1)/MI_params.timeJump-2
                                eTime2 = sTime2 + MI_params.windowSize - 1;
                                
                                if timeStep2 >= timeStep1 % only compute info when feature time is after neural time
                                    Y = selReachFeature(sTime2:eTime2,nonNaNtrials);
    %                                 Y = repmat(Y,size(X,1),1);
                                    Y = Y(:);

                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.spikes_early_late{early_late}.(params.areas{aidx}).(params.instReachFeatures{fidx}).info(unit,timeStepCount1,timeStepCount2) = I1{1}(1);

                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff

                                            Y = selReachFeature(sTime2:eTime2,nonNaNtrials);
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            tmpShuffle=tmpShuffle(shuffIdxs(:,shIdx)<=length(nonNaNtrials));
                                            Y = Y(:,tmpShuffle);
    %                                         Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);                                            

                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.spikes_early_late{early_late}.(params.areas{aidx}).(params.instReachFeatures{fidx}).infoSh(unit,timeStepCount1,timeStepCount2,shIdx) = I1{1}(1);
                                        end
                                    end
                                end
                                
                                MIout{animal}.spikes_early_late{early_late}.(params.areas{aidx}).(params.instReachFeatures{fidx}).binTimes = tmpBinTimes;
                                sTime2 = sTime2 + MI_params.timeJump;
                                timeStepCount2 = timeStepCount2+1;
                                
                            end
                            
                            sTime1 = sTime1 + MI_params.timeJump;
                            timeStepCount1 = timeStepCount1+1;
                        end

                        
                    end
                end
            end
        end
    end
    
    % calculate spikes MI by day
    if MI_params.doSpikesday==1
        for day = 1:numel(spikes{animal})
            for aidx = 1:length(params.areas)
                
                if MI_params.doShuff
                    shuffIdxs = zeros(size(spikes{animal}{day}.(params.areas{aidx}),3),MI_params.nShuff);
                    for shIdx = 1:MI_params.nShuff
                        shuffIdxs(:,shIdx) = randperm(size(spikes{animal}{day}.(params.areas{aidx}),3));
                    end
                end
                
                for fidx = rf_idxs
                    selReachFeature = reachFeatures{animal}{day}.instant.(params.instReachFeatures{fidx})(sel_rf_TW(1):sel_rf_TW(2),:);

                    for unit = 1:size(spikes{animal}{day}.(params.areas{aidx}),1)
                        
                        disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | unit ', num2str(unit), ' | day ' num2str(day) ' of ' num2str(numel(spikes{1,1})) ' | reach feature: ', params.instReachFeatures{fidx}])
                        
                        spikeActivity = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:));
                        nonNaNtrials = find(~isnan(sum(reachFeatures{animal}{day}.instant.(params.instReachFeatures{fidx}),1)));
                        spikeActivity = spikeActivity(sel_neural_TW(1):sel_neural_TW(2),nonNaNtrials);
                        
                        tmpBinTimes = [];
%                         newBinTimes = [];
%                         newBinTimes_totSpikes = [];
                        timeStepCount1 = 1;
                        sTime1 = 1;
                        for timeStep1 = 1:size(spikeActivity,1)/MI_params.timeJump-2 % loop over neural activity time
                            eTime1 = sTime1 + MI_params.windowSize - 1;
                            X = spikeActivity(sTime1:eTime1,:);
                            X = X(:);
                            tmpBinTimes = [tmpBinTimes (((sTime1-1)*binSize)+tMin)-5001]; % not creating timeBins1 and timeBins2 since they would be the same

                            timeStepCount2 = 1;
                            sTime2 = 1;
                            for timeStep2 = 1:size(selReachFeature,1)/MI_params.timeJump-2
                                eTime2 = sTime2 + MI_params.windowSize - 1;

                                if timeStep2 >= timeStep1 % only compute info when feature time is after neural time
                                    Y = selReachFeature(sTime2:eTime2,nonNaNtrials);

%                                     Y = repmat(Y,size(X,1),1);
                                    Y = Y(:);

                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.instReachFeatures{fidx}).info(unit,timeStepCount1,timeStepCount2) = I1{1}(1);

                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff

                                            Y = selReachFeature(sTime2:eTime2,nonNaNtrials);
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            tmpShuffle=tmpShuffle(shuffIdxs(:,shIdx)<=length(nonNaNtrials));
                                            Y = Y(:,tmpShuffle);
%                                             Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);                                            

                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.instReachFeatures{fidx}).infoSh(unit,timeStepCount1,timeStepCount2,shIdx) = I1{1}(1);
                                        end

                                    end
                                end

                                MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.instReachFeatures{fidx}).binTimes = tmpBinTimes;
                                sTime2 = sTime2 + MI_params.timeJump;
                                timeStepCount2 = timeStepCount2+1;
                            end

                            sTime1 = sTime1 + MI_params.timeJump;
                            timeStepCount1 = timeStepCount1+1;
                        end
                        
                    end
                end
            end
        end
        
    end
    
end

% save results
if MI_params.saveMI
    
    save([MI_params.save_path 'MIinst_all_indTpoints_' date ...
        '_shuff' num2str(MI_params.doLFPearlylate) ...
        '_LFPearlylate_' num2str(MI_params.doLFPearlylate) ...
        '_LFPday_' num2str(MI_params.doLFPday) ...
        '_Spikesearlylate_' num2str(MI_params.doSpikesearlylate) ...
        '_Spikesday_' num2str(MI_params.doSpikesday) ...
        '_SpikesdayAve_' num2str(MI_params.doSpikesdayAve) ...
        '_MatchTrials_' num2str(rf_params.match_trials) ...
        '.mat'],'-v7.3','MIout')

%     save([MI_params.save_path 'MI_' date ...
%         '_LFPearlylate_' num2str(MI_params.doLFPearlylate) ...
%         '_LFPday_' num2str(MI_params.doLFPday) ...
%         '_Spikesearlylate_' num2str(MI_params.doSpikesearlylate) ...
%         '_Spikesday_' num2str(MI_params.doSpikesday) ...
%         '_SpikesdayAve_' num2str(MI_params.doSpikesdayAve) ...
%         '_MatchTrials_' num2str(rf_params.match_trials) ...
%         '_25Hzhighpass.mat'],'-v7.3','MIout')
        
end

end
