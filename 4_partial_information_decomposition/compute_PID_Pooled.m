function [PIDout] = compute_PID_Pooled(LFP, spikes, spikes_ave, reachFeatures, LFP_early_late, spikes_early_late, reachFeatures_early_late, params, PID_params, tMin, binSize, animals_to_run)

% load MI results
MI_info = load(PID_params.path_to_MI);
    
%% calculate LFP PID by early/late
if PID_params.doLFPearlylate==1
    for animal = animals_to_run
        
        clearvars PIDout
        
        % determine channels with significant mutual information
        sig_channels = cell(length(params.lfpFeatures),length(params.reachFeatures),length(params.areas),2);
        for n_lfp = 1:length(params.lfpFeatures)
            for fidx = 1:length(params.reachFeatures)
                for aidx = 1:length(params.areas)
                    for early_late = 1:2
                        sig_chan = [];
                        for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                            infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                            infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                            tmp_sig = clusterStat(infQuant,infQuantSh,PID_params.MIclusterparams(1),PID_params.MIclusterparams(2));
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
        
        for fidx = 1:length(params.reachFeatures)
            for early_late = 1:2
                for n_lfp = 1:length(params.lfpFeatures)
                    
                    nonNaNtrials = find(~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})));
                    
                    if PID_params.doShuff
                        shuffIdxs = zeros(length(nonNaNtrials),PID_params.nShuff);
                        for shidx = 1:PID_params.nShuff
                            shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
                        end
                    end
                    
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
                        
                        nonNaNtrials = find(~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})));
                        tmp_M1_LFP = tmp_M1_LFP(:,nonNaNtrials);
                        tmp_DLS_LFP = tmp_DLS_LFP(:,nonNaNtrials);
                        
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
                            
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
                            delayIdx = 1;
                            
                            for delay = tmp_delays
                                
                                eTime = sTime + PID_params.windowSize-1;
                                
                                if delay <= 0 % M1 is the receiver (DLS prior to M1)
                                    X1 = tmp_M1_LFP(sTime:eTime,:);
                                    X2 = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                elseif delay > 0 % DLS is the receiver (M1 prior to DLS)
                                    X1 = tmp_DLS_LFP(sTime:eTime,:);
                                    X2 = tmp_M1_LFP(sTime-delay:eTime-delay,:);
                                end

                                Y = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                Y = repmat(Y,size(X1,1),1);
                                
                                X1 = X1(:);
                                X2 = X2(:);
                                Y = Y(:);

                                
                                pair_shared(timeStep,delayIdx) = I.shared;
                                pair_unique_m1(timeStep,delayIdx) = I.uniqueX1;
                                pair_unique_dls(timeStep,delayIdx) = I.uniqueX2;
                                pair_complementary(timeStep,delayIdx) = I.complementary;
                                
                                if PID_params.doShuff
                                    for shidx = 1:PID_params.nShuff
                                        
                                        if delay <= 0 % M1 is the receiver (DLS prior to M1)
                                            X1 = tmp_M1_LFP(sTime:eTime,:);
                                            X2 = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                        elseif delay > 0 % DLS is the receiver (M1 prior to DLS)
                                            X1 = tmp_DLS_LFP(sTime:eTime,:);
                                            X2 = tmp_M1_LFP(sTime-delay:eTime-delay,:);
                                        end

                                        Y = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
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

%% calculate Spiking PID by early/late
% if PID_params.doSpikesearlylate==1
%     for animal = animals_to_run
% 
%         clearvars PIDout
% 
%         for fidx = 1:length(params.reachFeatures)
%             for early_late = 1:2
%                 for n_lfp = 1:length(params.lfpFeatures)
% 
%                     nonNaNtrials = find(~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})));
% 
%                     if PID_params.doShuff
%                         shuffIdxs = zeros(length(nonNaNtrials),PID_params.nShuff);
%                         for shidx = 1:PID_params.nShuff
%                             shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
%                         end
%                     end
% 
%                     tmp_delays = -25:25;
%                     tmp_timeStep = 1:floor((length(26:225)-(PID_params.windowSize-PID_params.timeJump))/PID_params.timeJump);
% 
%                     tmp_shared = zeros(1,length(tmp_timeStep),length(tmp_delays));
%                     tmp_unique_m1 = zeros(1,length(tmp_timeStep),length(tmp_delays));
%                     tmp_unique_dls = zeros(1,length(tmp_timeStep),length(tmp_delays));
%                     tmp_complementary = zeros(1,length(tmp_timeStep),length(tmp_delays));
% 
%                     if PID_params.doShuff
%                         tmp_shared_Sh = zeros(1,length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                         tmp_unique_m1_Sh = zeros(1,length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                         tmp_unique_dls_Sh = zeros(1,length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                         tmp_complementary_Sh = zeros(1,length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                     end
% 
%                     for pairIdx = 1
% 
%                         disp(['Computing: ' params.animals{animal} ' | ' num2str(early_late) ' of 2 (early/late) | ' params.reachFeatures{fidx}])
% 
%                         tmp_M1_spikes = squeeze(spikes_early_late{animal}{early_late}.M1(1,:,:));
%                         tmp_DLS_spikes = squeeze(spikes_early_late{animal}{early_late}.DLS(1,:,:));
% 
%                         nonNaNtrials = find(~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})));
%                         tmp_M1_spikes = tmp_M1_spikes(:,nonNaNtrials);
%                         tmp_DLS_spikes = tmp_DLS_spikes(:,nonNaNtrials);
% 
%                         pair_shared = zeros(length(tmp_timeStep),length(tmp_delays));
%                         pair_unique_m1 = zeros(length(tmp_timeStep),length(tmp_delays));
%                         pair_unique_dls = zeros(length(tmp_timeStep),length(tmp_delays));
%                         pair_complementary = zeros(length(tmp_timeStep),length(tmp_delays));
% 
%                         if PID_params.doShuff
%                             pair_shared_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                             pair_unique_m1_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                             pair_unique_dls_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                             pair_complementary_Sh = zeros(length(tmp_timeStep),length(tmp_delays),PID_params.nShuff);
%                         end
% 
%                         tmpBinTimes = [];
%                         sTime = 26;
% 
%                         for timeStep = 1:length(tmp_timeStep)
% 
%                             tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
%                             delayIdx = 1;
% 
%                             for delay = tmp_delays
% 
%                                 eTime = sTime + PID_params.windowSize-1;
% 
%                                 X1 = tmp_M1_spikes(sTime:eTime,:);
%                                 X2 = tmp_DLS_spikes(sTime+delay:eTime+delay,:);
% 
%                                 Y = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
%                                 Y = repmat(Y,size(X1,1),1);
% 
%                                 X1 = X1(:);
%                                 X2 = X2(:);
%                                 Y = Y(:);
% 
%                                 I = PID(Y',X1',X2',PID_params.opts);
% 
%                                 pair_shared(timeStep,delayIdx) = I.shared;
%                                 pair_unique_m1(timeStep,delayIdx) = I.uniqueX1;
%                                 pair_unique_dls(timeStep,delayIdx) = I.uniqueX2;
%                                 pair_complementary(timeStep,delayIdx) = I.complementary;
% 
%                                 if PID_params.doShuff
%                                     for shidx = 1:PID_params.nShuff
% 
%                                         X1 = tmp_M1_LFP(sTime:eTime,:);
%                                         X2 = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
% 
%                                         Y = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
%                                         Y = Y(shuffIdxs(:,shidx));
%                                         Y = repmat(Y,size(X1,1),1);
% 
%                                         X1 = X1(:);
%                                         X2 = X2(:);
%                                         Y = Y(:);
% 
%                                         I = PID(Y',X1',X2',PID_params.opts);
% 
%                                         pair_shared_Sh(timeStep,delayIdx,shidx) = I.shared;
%                                         pair_unique_m1_Sh(timeStep,delayIdx,shidx) = I.uniqueX1;
%                                         pair_unique_dls_Sh(timeStep,delayIdx,shidx) = I.uniqueX2;
%                                         pair_complementary_Sh(timeStep,delayIdx,shidx) = I.complementary;
% 
%                                     end
%                                 end
%                                 delayIdx = delayIdx + 1;
%                             end
%                             sTime = sTime + PID_params.timeJump;
%                         end
% 
%                         tmp_shared(pairIdx,:,:) = pair_shared;
%                         tmp_unique_m1(pairIdx,:,:) = pair_unique_m1;
%                         tmp_unique_dls(pairIdx,:,:) = pair_unique_dls;
%                         tmp_complementary(pairIdx,:,:) = pair_complementary;
% 
%                         if PID_params.doShuff
%                             tmp_shared_Sh(pairIdx,:,:,:) = pair_shared_Sh;
%                             tmp_unique_m1_Sh(pairIdx,:,:,:) = pair_unique_m1_Sh;
%                             tmp_unique_dls_Sh(pairIdx,:,:,:) = pair_unique_dls_Sh;
%                             tmp_complementary_Sh(pairIdx,:,:,:) = pair_complementary_Sh;
%                         end
% 
%                     end
% 
%                     PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared = tmp_shared;
%                     PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1 = tmp_unique_m1;
%                     PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dls = tmp_unique_dls;
%                     PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementary = tmp_complementary;
% 
%                     if PID_params.doShuff
%                         PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh = tmp_shared_Sh;
%                         PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1Sh = tmp_unique_m1_Sh;
%                         PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dlsSh = tmp_unique_dls_Sh;
%                         PIDout.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementarySh = tmp_complementary_Sh;
%                     end
% 
%                 end
%             end
%         end
% 
%         % save results
%         if PID_params.savePID
%             save([PID_params.save_path 'PID_' params.animals{animal} '_' date ...
%                 '_LFPearlylate_' num2str(PID_params.doLFPearlylate) ...
%                 '_LFPday_' num2str(PID_params.doLFPday) ...
%                 '_Spikesday_' num2str(PID_params.doSpikesday) ...
%                 '_Pooled.mat'],'-v7.3','PIDout','params','PID_params')
%         end
% 
%     end
% end                                  

%% calculate Spikes PID by day Ave

% if PID_params.doSpikesdayAve==1
%     for animal = animals_to_run
% 
%         tmpPIDout = cell(length(params.reachFeatures),length(MI_info.MIout{animal}.spikesAveDay));
% 
%         for fidx = 1:length(params.reachFeatures)        
%             for day = 1:length(MI_info.MIout{animal}.spikesAveDay)
% 
%                 nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
%                 if PID_params.doShuff
%                     shuffIdxs = zeros(length(nonNaNtrials),PID_params.nShuff);
%                     for shidx = 1:PID_params.nShuff
%                         shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
%                     end
%                 end
% 
%                 pair_count = 1;
%                 for m1_unit = 1
%                     for dls_unit = 1
% 
%                         disp(['Computing... ' params.animals{animal} ' | day ' num2str(day) ' | ' params.reachFeatures{fidx}])
% 
%                         tmp_M1_spiking = squeeze(spikes_ave{animal}{day}.M1(m1_unit,:,:));
%                         tmp_DLS_spiking = squeeze(spikes_ave{animal}{day}.DLS(dls_unit,:,:));
% 
%                         nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
%                         tmp_M1_spiking = tmp_M1_spiking(:,nonNaNtrials);
%                         tmp_DLS_spiking = tmp_DLS_spiking(:,nonNaNtrials);
% 
%                         tmpBinTimes = [];
%                         sTime = 1;
%                         timestep_count = 1;
%                         for timeStep = 1:size(tmp_M1_spiking,1)/PID_params.timeJump
% 
%                             eTime = sTime + PID_params.windowSize - 1;
%                             if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(tmp_M1_spiking,1)
%                                 sTime = sTime + PID_params.timeJump;
%                             else
%                                 delay_idx = 1;
%                                 for delay = -25:25
%                                     eTime = sTime + PID_params.windowSize - 1;
% 
%                                     X1 = tmp_M1_spiking(sTime:eTime,:);
%                                     X2 = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
% 
%                                     Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
%                                     Y = repmat(Y,size(X1,1),1);
% 
%                                     X1 = X1(:);
%                                     X2 = X2(:);
%                                     Y = Y(:);
% 
%                                     I = PID(Y',X1',X2',PID_params.opts);
% 
%                                     tmpPIDout{fidx}{day}.shared(pair_count,timestep_count,delay_idx) = I.shared;
%                                     tmpPIDout{fidx}{day}.unique_m1(pair_count,timestep_count,delay_idx) = I.uniqueX1;
%                                     tmpPIDout{fidx}{day}.unique_dls(pair_count,timestep_count,delay_idx) = I.uniqueX2;
%                                     tmpPIDout{fidx}{day}.complementary(pair_count,timestep_count,delay_idx) = I.complementary;
% 
%                                     if PID_params.doShuff
% 
%                                         tmp_shared = zeros(1,PID_params.nShuff);
%                                         tmp_unique_m1 = zeros(1,PID_params.nShuff);
%                                         tmp_unique_dls = zeros(1,PID_params.nShuff);
%                                         tmp_complementary = zeros(1,PID_params.nShuff);
% 
%                                         parfor shidx = 1:PID_params.nShuff
% 
%                                             X1 = tmp_M1_spiking(sTime:eTime,:);
%                                             X2 = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
% 
%                                             Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
%                                             Y = Y(shuffIdxs(:,shidx));
%                                             Y = repmat(Y,size(X1,1),1);
% 
%                                             X1 = X1(:);
%                                             X2 = X2(:);
%                                             Y = Y(:);
% 
%                                             I = PID(Y',X1',X2',PID_params.opts);
% 
%                                             tmp_shared(shidx) = I.shared;
%                                             tmp_unique_m1(shidx) = I.uniqueX1;
%                                             tmp_unique_dls(shidx) = I.uniqueX2;
%                                             tmp_complementary(shidx) = I.complementary;
%                                         end
% 
%                                         tmpPIDout{fidx}{day}.sharedSh(pair_count,timestep_count,delay_idx,:) = tmp_shared;
%                                         tmpPIDout{fidx}{day}.unique_m1Sh(pair_count,timestep_count,delay_idx,:) = tmp_unique_m1;
%                                         tmpPIDout{fidx}{day}.unique_dlsSh(pair_count,timestep_count,delay_idx,:) = tmp_unique_dls;
%                                         tmpPIDout{fidx}{day}.complementarySh(pair_count,timestep_count,delay_idx,:) = tmp_complementary;
% 
%                                     end
% 
%                                     delay_idx = delay_idx + 1;
%                                 end
%                                 timestep_count = timestep_count + 1;
%                                 sTime = sTime + PID_params.timeJump;
%                             end
%                         end
%                     end
%                     pair_count = pair_count + 1;
%                 end
%             end
%         end
% 
%         for fidx = 1:length(params.reachFeatures)
%             for day = 1:length(MI_info.MIout{animal}.spikesAveDay)
%                 if ~isempty(tmpPIDout{fidx})
%                     if day<=length(tmpPIDout{fidx})
%                         if ~isempty(tmpPIDout{fidx}{day})
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).shared = tmpPIDout{fidx}{day}.shared;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).unique_m1 = tmpPIDout{fidx}{day}.unique_m1;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).unique_dls = tmpPIDout{fidx}{day}.unique_dls;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).complementary = tmpPIDout{fidx}{day}.complementary;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).sharedSh = tmpPIDout{fidx}{day}.sharedSh;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).unique_m1Sh = tmpPIDout{fidx}{day}.unique_m1Sh;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).unique_dlsSh = tmpPIDout{fidx}{day}.unique_dlsSh;
%                             PIDout{animal}.spikesAveDay{day}.(params.reachFeatures{fidx}).complementarySh = tmpPIDout{fidx}{day}.complementarySh;
%                         end
%                     end
%                 end
%             end
%         end
% 
%         % save results
%         if PID_params.savePID
%             save([PID_params.save_path 'PID_' params.animals{animal} '_' date ...
%                 '_LFPearlylate_' num2str(PID_params.doLFPearlylate) ...
%                 '_LFPday_' num2str(PID_params.doLFPday) ...
%                 '_Spikesday_' num2str(PID_params.doSpikesday) ...
%                 '_SpikesdayAve_' num2str(PID_params.doSpikesdayAve) ...
%                 '_Pooled.mat'],'-v7.3','PIDout','params','PID_params')
%         end
% 
%     end
% end

%% calculate LFP PID by day
% if PID_params.doLFPday==1
%     for animal = 1:numel(params.animals)
% 
%         % determine channels with significant mutual information
%         sig_channels = cell(length(params.lfpFeatures),length(params.reachFeatures),length(params.areas),length(MI_info.MIout{animal}.LFPday));
%         for n_lfp = 1:length(params.lfpFeatures)
%             for fidx = 1:length(params.reachFeatures)
%                 for aidx = 1:length(params.areas)
%                     for day = 1:length(MI_info.MIout{animal}.LFPday)
%                         sig_chan = [];
%                         for chan = 1:size(MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                             infQuant = MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
%                             infQuantSh = MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
%                             tmp_sig = clusterStat(infQuant,infQuantSh,PID_params.MIclusterparams(1),PID_params.MIclusterparams(2));
%                             infQuant(~tmp_sig) = 0;
%                             [~,i] = max(infQuant);
%                             sigTime_start = find(MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(1));
%                             sigTime_end = find(MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(2));
%                             sig_chan = [sig_chan i>=sigTime_start & i<=sigTime_end];
%                         end
%                         sig_channels{n_lfp,fidx,aidx,day} = sig_chan;
%                     end
%                 end
%             end
%         end
%         
%         tmpPIDout = cell(length(params.reachFeatures));
%         for fidx = 1:length(params.reachFeatures)        
%             for day = 1:length(MI_info.MIout{animal}.LFPday)
%                 for n_lfp = 1:length(params.lfpFeatures)
%                     
%                     if PID_params.doShuff
%                         shuffIdxs = zeros(size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).(params.areas{1}),3),PID_params.nShuff);
%                         for shidx = 1:PID_params.nShuff
%                             shuffIdxs(:,shidx) = randperm(size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).(params.areas{1}),3));
%                         end
%                     end
%                     
%                     pair_count = 1;
%                     for m1_chan = 1:size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).M1,1)
%                         for dls_chan = 1:size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).DLS,1)
%                             if sig_channels{n_lfp,fidx,1,day}(m1_chan) && sig_channels{n_lfp,fidx,2,day}(dls_chan)
%                                 disp(['Computing: ' params.animals{animal} ' | ' params.lfpFeatures{n_lfp} ' | M1 channel ' num2str(m1_chan) ' & DLS channel ' num2str(dls_chan) ' | day ' num2str(day) ' of ' num2str(length(MI_info.MIout{animal}.LFPday)) ' | ' params.reachFeatures{fidx}])
%                                 tmp_M1_spiking = squeeze(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).M1(m1_chan,:,:));
%                                 tmp_DLS_spiking = squeeze(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).DLS(dls_chan,:,:));
%                                 
%                                 tmpBinTimes = [];
%                                 sTime = 1;
%                                 timestep_count = 1;
%                                 for timeStep = 1:size(tmp_M1_spiking,1)/PID_params.timeJump
%                                     
%                                     eTime = sTime + PID_params.windowSize - 1;
%                                     if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(tmp_M1_spiking,1)
%                                         sTime = sTime + PID_params.timeJump;
%                                     else
%                                         delay_idx = 1;
%                                         for delay = -25:25
%                                             eTime = sTime + PID_params.windowSize - 1;
%                                             X1 = tmp_M1_spiking(sTime:eTime,:);
%                                             X2 = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
%                                             
%                                             Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
%                                             Y = repmat(Y,size(X1,1),1);
%                                             
%                                             X1 = X1(:);
%                                             X2 = X2(:);
%                                             Y = Y(:);
%                                             
%                                             I = PID(Y',X1',X2',PID_params.opts);
%                                             
%                                             tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).shared(pair_count,timestep_count,delay_idx) = I.shared;
%                                             tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_m1(pair_count,timestep_count,delay_idx) = I.uniqueX1;
%                                             tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_dls(pair_count,timestep_count,delay_idx) = I.uniqueX2;
%                                             tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).complementary(pair_count,timestep_count,delay_idx) = I.complementary;
%                                             
%                                             if PID_params.doShuff
%                                                 
%                                                 tmp_shared = zeros(1,PID_params.nShuff);
%                                                 tmp_unique_m1 = zeros(1,PID_params.nShuff);
%                                                 tmp_unique_dls = zeros(1,PID_params.nShuff);
%                                                 tmp_complementary = zeros(1,PID_params.nShuff);
%                                                 
%                                                 for shidx = 1:PID_params.nShuff
%                                                     Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(shuffIdxs(:,shidx));
%                                                     Y = repmat(Y,PID_params.windowSize,1);
%                                                     Y = Y(:);
%                                                     
%                                                     I = PID(Y',X1',X2',PID_params.opts);
%                                                     
%                                                     tmp_shared(shidx) = I.shared;
%                                                     tmp_unique_m1(shidx) = I.uniqueX1;
%                                                     tmp_unique_dls(shidx) = I.uniqueX2;
%                                                     tmp_complementary(shidx) = I.complementary;
%                                                 end
%                                                 
%                                                 tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).sharedSh(pair_count,timestep_count,delay_idx,:) = tmp_shared;
%                                                 tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_m1Sh(pair_count,timestep_count,delay_idx,:) = tmp_unique_m1;
%                                                 tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_dlsSh(pair_count,timestep_count,delay_idx,:) = tmp_unique_dls;
%                                                 tmpPIDout{fidx}{day}.(params.lfpFeatures{n_lfp}).complementarySh(pair_count,timestep_count,delay_idx,:) = tmp_complementary;
%                                                 
%                                             end
%                                             
%                                             delay_idx = delay_idx + 1;
%                                         end
%                                         timestep_count = timestep_count + 1;
%                                         sTime = sTime + PID_params.timeJump;
%                                     end
%                                 end
%                             end
%                             pair_count = pair_count + 1;
%                         end
%                     end
%                     
%                 end
%             end
%         end
%                                
%         for fidx = 1:length(params.reachFeatures)
%             for day = 1:length(MI_info.MIout{animal}.LFPday)
%                 for n_lfp = 1:length(params.lfpFeatures)
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).shared;
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1 = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_m1;
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dls = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_dls;
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementary = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).complementary;            
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).sharedSh;
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1Sh = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_m1Sh;
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dlsSh = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_dlsSh;
%                     PIDout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementarySh = tmpPIDout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).complementarySh;                
%                 end
%             end
%         end
%       
%     end
% end

%% calculate Spikes PID by day
if PID_params.doSpikesday==1
    for animal = animals_to_run
        
        % determine channels with significant mutual information
        sig_units = cell(length(params.reachFeatures),length(params.areas),length(MI_info.MIout{animal}.spikesDay));           
        for fidx = 1:length(params.reachFeatures)
            for aidx = 1:length(params.areas)
                for day = 1:length(MI_info.MIout{animal}.spikesDay)
                    sig_unit = [];
                    for unit = 1:size(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                        infQuantSh = MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                        tmp_sig = clusterStat(infQuant,infQuantSh,PID_params.MIclusterparams(1),PID_params.MIclusterparams(2));
                        infQuant(~tmp_sig) = 0;
                        [~,i] = max(infQuant);
                        sigTime_start = find(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(1));
                        sigTime_end = find(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(2));
                        sig_unit = [sig_unit i>=sigTime_start & i<=sigTime_end];
                    end
                    sig_units{fidx,aidx,day} = sig_unit;
                end
            end
        end

        tmpPIDout = cell(length(params.reachFeatures),length(MI_info.MIout{animal}.spikesDay));
        
        for fidx = 1:length(params.reachFeatures)        
            for day = 1:length(MI_info.MIout{animal}.spikesDay)

                nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
                if PID_params.doShuff
                    shuffIdxs = zeros(length(nonNaNtrials),PID_params.nShuff);
                    for shidx = 1:PID_params.nShuff
                        shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
                    end
                end

                pair_count = 1;
                for m1_unit = 1:size(MI_info.MIout{animal}.spikesDay{day}.M1.(params.reachFeatures{fidx}).info,1)
                    for dls_unit = 1:size(MI_info.MIout{animal}.spikesDay{day}.DLS.(params.reachFeatures{fidx}).info,1)
                         if sig_units{fidx,1,day}(m1_unit) && sig_units{fidx,2,day}(dls_unit)
                            disp(['Computing... ' params.animals{animal} ' | M1 unit ' num2str(m1_unit) ' & DLS unit ' num2str(dls_unit) ' | day ' num2str(day) ' | ' params.reachFeatures{fidx}])

                            tmp_M1_spiking = squeeze(spikes{animal}{day}.M1(m1_unit,:,:));
                            tmp_DLS_spiking = squeeze(spikes{animal}{day}.DLS(dls_unit,:,:));
                            
                            nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
                            tmp_M1_spiking = tmp_M1_spiking(:,nonNaNtrials);
                            tmp_DLS_spiking = tmp_DLS_spiking(:,nonNaNtrials);
                            
                            tmpBinTimes = [];
                            sTime = 1;
                            timestep_count = 1;
                            for timeStep = 1:size(tmp_M1_spiking,1)/PID_params.timeJump

                                eTime = sTime + PID_params.windowSize - 1;
                                if sTime-PID_params.delay<1 || eTime+PID_params.delay>size(tmp_M1_spiking,1)
                                    sTime = sTime + PID_params.timeJump;
                                else
                                    delay_idx = 1;
                                    for delay = -25:25
                                        eTime = sTime + PID_params.windowSize - 1;

                                        if delay <= 0 % M1 is the receiver (DLS prior to M1)
                                            X1 = tmp_M1_spiking(sTime:eTime,:);
                                            X2 = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
                                        elseif delay > 0 % DLS is the receiver (M1 prior to DLS)
                                            X1 = tmp_DLS_spiking(sTime:eTime,:);
                                            X2 = tmp_M1_spiking(sTime-delay:eTime-delay,:);
                                        end

                                        Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
                                        Y = repmat(Y,size(X1,1),1);

                                        X1 = X1(:);
                                        X2 = X2(:);
                                        Y = Y(:);

                                        I = PID(Y',X1',X2',PID_params.opts);

                                        tmpPIDout{fidx}{day}.shared(pair_count,timestep_count,delay_idx) = I.shared;
                                        tmpPIDout{fidx}{day}.unique_m1(pair_count,timestep_count,delay_idx) = I.uniqueX1;
                                        tmpPIDout{fidx}{day}.unique_dls(pair_count,timestep_count,delay_idx) = I.uniqueX2;
                                        tmpPIDout{fidx}{day}.complementary(pair_count,timestep_count,delay_idx) = I.complementary;

                                        if PID_params.doShuff

                                            tmp_shared = zeros(1,PID_params.nShuff);
                                            tmp_unique_m1 = zeros(1,PID_params.nShuff);
                                            tmp_unique_dls = zeros(1,PID_params.nShuff);
                                            tmp_complementary = zeros(1,PID_params.nShuff);

                                            parfor shidx = 1:PID_params.nShuff

                                                if delay <= 0 % M1 is the receiver (DLS prior to M1)
                                                    X1 = tmp_M1_LFP(sTime:eTime,:);
                                                    X2 = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                                elseif delay > 0 % DLS is the receiver (M1 prior to DLS)
                                                    X1 = tmp_DLS_LFP(sTime:eTime,:);
                                                    X2 = tmp_M1_LFP(sTime-delay:eTime-delay,:);
                                                end
                                                
                                                Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
                                                Y = Y(shuffIdxs(:,shidx));
                                                Y = repmat(Y,size(X1,1),1);
                                                
                                                X1 = X1(:);
                                                X2 = X2(:);
                                                Y = Y(:);
                                                
                                                I = PID(Y',X1',X2',PID_params.opts);
                                                    
                                                tmp_shared(shidx) = I.shared;
                                                tmp_unique_m1(shidx) = I.uniqueX1;
                                                tmp_unique_dls(shidx) = I.uniqueX2;
                                                tmp_complementary(shidx) = I.complementary;
                                            end

                                            tmpPIDout{fidx}{day}.sharedSh(pair_count,timestep_count,delay_idx,:) = tmp_shared;
                                            tmpPIDout{fidx}{day}.unique_m1Sh(pair_count,timestep_count,delay_idx,:) = tmp_unique_m1;
                                            tmpPIDout{fidx}{day}.unique_dlsSh(pair_count,timestep_count,delay_idx,:) = tmp_unique_dls;
                                            tmpPIDout{fidx}{day}.complementarySh(pair_count,timestep_count,delay_idx,:) = tmp_complementary;

                                        end

                                        delay_idx = delay_idx + 1;
                                    end
                                    timestep_count = timestep_count + 1;
                                    sTime = sTime + PID_params.timeJump;
                                end
                            end
                        end
                        pair_count = pair_count + 1;
                    end
                end
                
            end
        end
        
        for fidx = 1:length(params.reachFeatures)
            for day = 1:length(MI_info.MIout{animal}.spikesDay)
                if ~isempty(tmpPIDout{fidx})
                    if day<=length(tmpPIDout{fidx})
                        if ~isempty(tmpPIDout{fidx}{day})
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared = tmpPIDout{fidx}{day}.shared;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).unique_m1 = tmpPIDout{fidx}{day}.unique_m1;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).unique_dls = tmpPIDout{fidx}{day}.unique_dls;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).complementary = tmpPIDout{fidx}{day}.complementary;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).sharedSh = tmpPIDout{fidx}{day}.sharedSh;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).unique_m1Sh = tmpPIDout{fidx}{day}.unique_m1Sh;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).unique_dlsSh = tmpPIDout{fidx}{day}.unique_dlsSh;
                            PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).complementarySh = tmpPIDout{fidx}{day}.complementarySh;
                        end
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
                '_Pooled_neg600to0_maxVel.mat'],'-v7.3','PIDout','params','PID_params')
        end

    end
end       

end