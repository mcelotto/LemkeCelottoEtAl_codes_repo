function [MIout] = compute_mutualInformation_successFail(LFP_early_late_noTrialMatch, reachFeatures_early_late_noTrialMatch, spikes, reachFeatures, params, MI_params, PID_params, tMin, binSize, save_params, features_to_compute, bin_params)

parfor animal = 1:numel(params.animals)
  
    % calculate LFP MI by early/late
    if MI_params.doLFP==1
        for early_late = 1:2
            for n_lfp = 1:length(params.lfpFeatures)
                for aidx = 1:length(params.areas)                   
                    for fidx = features_to_compute

                        early_success = reachFeatures_early_late_noTrialMatch{animal}{1}.success_rate;
                        late_success = reachFeatures_early_late_noTrialMatch{animal}{2}.success_rate;
                        earlySuccessIdx = find(early_success==1);
                        lateSuccessIdx = find(late_success==1);                        
                        early_nonNaNtrials = find(~isnan(reachFeatures_early_late_noTrialMatch{animal}{1}.(params.reachFeatures{fidx})));
                        late_nonNaNtrials = find(~isnan(reachFeatures_early_late_noTrialMatch{animal}{2}.(params.reachFeatures{fidx})));
    
                        earlySuccessIdx = intersect(earlySuccessIdx,early_nonNaNtrials);
                        lateSuccessIdx = intersect(lateSuccessIdx,late_nonNaNtrials);
                        tmp_min = min([length(earlySuccessIdx) length(lateSuccessIdx)]);
                        earlySuccessIdx = earlySuccessIdx(1:tmp_min);
                        lateSuccessIdx = lateSuccessIdx(end-tmp_min+1:end);
                        success_trials = {earlySuccessIdx,lateSuccessIdx};

                        early_success = reachFeatures_early_late_noTrialMatch{animal}{1}.success_rate;
                        late_success = reachFeatures_early_late_noTrialMatch{animal}{2}.success_rate;
                        earlyFailIdx = find(early_success==0);
                        lateFailIdx = find(late_success==0);                        
                        early_nonNaNtrials = find(~isnan(reachFeatures_early_late_noTrialMatch{animal}{1}.(params.reachFeatures{fidx})));
                        late_nonNaNtrials = find(~isnan(reachFeatures_early_late_noTrialMatch{animal}{2}.(params.reachFeatures{fidx})));
    
                        earlyFailIdx = intersect(earlyFailIdx,early_nonNaNtrials);
                        lateFailIdx = intersect(lateFailIdx,late_nonNaNtrials);
                        tmp_min = min([length(earlyFailIdx) length(lateFailIdx)]);
                        earlyFailIdx = earlyFailIdx(1:tmp_min);
                        lateFailIdx = lateFailIdx(end-tmp_min+1:end);
                        fail_trials = {earlyFailIdx,lateFailIdx};

                        %%% SUCCESS

                        if MI_params.doShuff
                            shuffIdxs = zeros(length(earlySuccessIdx),MI_params.nShuff);
                            for shIdx = 1:MI_params.nShuff
                                shuffIdxs(:,shIdx) = randperm(length(earlySuccessIdx));
                            end
                        end

                        for chan = 1:size(LFP_early_late_noTrialMatch{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | ', params.lfpFeatures{n_lfp} ,' channel ', num2str(chan), ' | ' num2str(early_late) ' of 2 (early/late) | reach feature: ', params.reachFeatures{fidx}])

                            LFPsignal = squeeze(LFP_early_late_noTrialMatch{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:));
                            LFPsignal = LFPsignal(:,success_trials{early_late});
                            
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(LFPsignal,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(LFPsignal,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = LFPsignal(sTime:eTime,:);
                                    Y = reachFeatures_early_late_noTrialMatch{animal}{early_late}.(params.reachFeatures{fidx})(success_trials{early_late});
                                    Y = repmat(Y,size(X,1),1);
                                    
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess(chan,timeStepCount) = I1{1}(1);
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];

                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            
                                            Y = reachFeatures_early_late_noTrialMatch{animal}{early_late}.(params.reachFeatures{fidx})(success_trials{early_late});
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);                                            
                                            
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoShSuccess(chan,timeStepCount,shIdx) = I1{1}(1);
                                        end
                                    end
                                    sTime = sTime + MI_params.timeJump;
                                    timeStepCount = timeStepCount+1;
                                end
                            end
                        end

                        %%% FAIL

                        if MI_params.doShuff
                            shuffIdxs = zeros(length(earlyFailIdx),MI_params.nShuff);
                            for shIdx = 1:MI_params.nShuff
                                shuffIdxs(:,shIdx) = randperm(length(earlyFailIdx));
                            end
                        end

                        for chan = 1:size(LFP_early_late_noTrialMatch{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | ', params.lfpFeatures{n_lfp} ,' channel ', num2str(chan), ' | ' num2str(early_late) ' of 2 (early/late) | reach feature: ', params.reachFeatures{fidx}])

                            LFPsignal = squeeze(LFP_early_late_noTrialMatch{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:));
                            LFPsignal = LFPsignal(:,fail_trials{early_late});
                            
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(LFPsignal,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(LFPsignal,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = LFPsignal(sTime:eTime,:);
                                    Y = reachFeatures_early_late_noTrialMatch{animal}{early_late}.(params.reachFeatures{fidx})(fail_trials{early_late});
                                    Y = repmat(Y,size(X,1),1);
                                    
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail(chan,timeStepCount) = I1{1}(1);
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];

                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            
                                            Y = reachFeatures_early_late_noTrialMatch{animal}{early_late}.(params.reachFeatures{fidx})(fail_trials{early_late});
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);                                            
                                            
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoShFail(chan,timeStepCount,shIdx) = I1{1}(1);
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

                    day_success = reachFeatures{animal}{day}.success_rate;
                    day_success_idx = find(day_success==1);
                    day_fail_idx = find(day_success==0);
                    day_nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
                    day_success_idx = intersect(day_success_idx,day_nonNaNtrials);
                    day_fail_idx = intersect(day_fail_idx,day_nonNaNtrials);

                    
                    %%% SUCCESS

                    if length(day_success_idx)>10
    
                        if MI_params.doShuff
                            shuffIdxs = zeros(length(day_success_idx),MI_params.nShuff);
                            for shIdx = 1:MI_params.nShuff
                                shuffIdxs(:,shIdx) = randperm(length(day_success_idx));
                            end
                        end
                    
                        for unit = 1:size(spikes{animal}{day}.(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | unit ', num2str(unit), ' | day ' num2str(day) ' of ' num2str(numel(spikes{1,1})) ' | reach feature: ', params.reachFeatures{fidx}])
                            
                            spikeActivity = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:));
                            spikeActivity = spikeActivity(:,day_success_idx);
    
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(spikeActivity,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(spikeActivity,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = spikeActivity(sTime:eTime,:);
                                    Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(day_success_idx);
                                    Y = repmat(Y,size(X,1),1);
                                    %%% SL 8.31.2022
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess(unit,timeStepCount) = I1{1}(1);
                                    
                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            
                                            Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(day_success_idx);
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);
                                            
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoShSuccess(unit,timeStepCount,shIdx) = I1{1}(1);
                                        end
                                    end
                                    
                                    sTime = sTime + MI_params.timeJump;
                                    timeStepCount = timeStepCount+1;
                                end
                            end
                        end

                    end

                    %%% FAIL

                    if length(day_fail_idx)>10

                        if MI_params.doShuff
                            shuffIdxs = zeros(length(day_fail_idx),MI_params.nShuff);
                            for shIdx = 1:MI_params.nShuff
                                shuffIdxs(:,shIdx) = randperm(length(day_fail_idx));
                            end
                        end
                    
                        for unit = 1:size(spikes{animal}{day}.(params.areas{aidx}),1)
                            
                            disp(['Computing... animal: ', params.animals{animal}, ' | area: ', params.areas{aidx}, ' | unit ', num2str(unit), ' | day ' num2str(day) ' of ' num2str(numel(spikes{1,1})) ' | reach feature: ', params.reachFeatures{fidx}])
                            
                            spikeActivity = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:));
                            spikeActivity = spikeActivity(:,day_fail_idx);
    
                            newBinTimes = [];
                            timeStepCount = 1;
                            sTime = 1;
                            for timeStep = 1:size(spikeActivity,1)/MI_params.timeJump
                                eTime = sTime + MI_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(spikeActivity,1)
                                    sTime = sTime + MI_params.timeJump;
                                else
                                    X = spikeActivity(sTime:eTime,:);
                                    Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(day_fail_idx);
                                    Y = repmat(Y,size(X,1),1);
                                    %%% SL 8.31.2022
                                    newBinTimes = [newBinTimes; tMin+(sTime*10)-5000 tMin+(eTime*10)-5000];
                                    X = X(:);
                                    Y = Y(:);
                                    
                                    I1 = information(X',Y',MI_params.opts,{'I'});
                                    MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail(unit,timeStepCount) = I1{1}(1);
                                    
                                    if MI_params.doShuff
                                        for shIdx = 1:MI_params.nShuff
                                            
                                            Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(day_fail_idx);
                                            tmpShuffle=shuffIdxs(:,shIdx);
                                            Y = Y(tmpShuffle);
                                            Y = repmat(Y,MI_params.windowSize,1);
                                            Y = Y(:);
                                            
                                            I1 = information(X',Y',MI_params.opts,{'I'});
                                            MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoShFail(unit,timeStepCount,shIdx) = I1{1}(1);
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
    save([save_params.save_path '\MI_successFail_' date ...
        '_LFPearlylate_' num2str(MI_params.doLFP) ...
        '_Spikes_' num2str(MI_params.doSpikes) ...
        '.mat'],'-v7.3','MIout')
end

end
