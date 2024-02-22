function [FITout] = compute_FIT(LFP, spikes, reachFeatures, LFP_early_late, reachFeatures_early_late, params, FIT_params, tMin, binSize,animals_to_run)

% load MI results
MI_info = load(FIT_params.path_to_MI);
    
% calculate LFP FIT by early/late
if FIT_params.doLFPearlylate==1
    for animal = animals_to_run
        
        clearvars FITout

        % determine channels with significant mutual information
        sig_channels = cell(length(params.lfpFeatures),length(params.reachFeatures),length(params.areas),2);
        for n_lfp = 1:length(params.lfpFeatures)
            for fidx = 1:length(params.reachFeatures)
                for aidx = 1:length(params.areas)
                    for early_late = 1:2
                        sig_chan = [];
                        for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                            
                            infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                            infQuant = temporal_rebinning(infQuant,5,'movmean',2);
                            
                            infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                            infQuantSh = temporal_rebinning(infQuantSh,5,'movmean',2);
                            
                            tmp_sig = clusterStat(infQuant,infQuantSh,FIT_params.MIclusterparams(1),FIT_params.MIclusterparams(2));
                            infQuant(~tmp_sig) = 0;
                            [~,i] = max(infQuant);
                            
                            rebinned_binTimes = temporal_rebinning(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes,5,'movmean',2);
                            sigTime_start = find(rebinned_binTimes==FIT_params.MIsigtiming(1));
                            sigTime_end = find(rebinned_binTimes==FIT_params.MIsigtiming(2));
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
                    
                    if FIT_params.doShuff
                        shuffIdxs = zeros(length(nonNaNtrials),FIT_params.nShuff);
                        for shidx = 1:FIT_params.nShuff
                            shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
                        end
                    end
                    
                    tmp_pairs = [...
                        repmat(find(sig_channels{n_lfp,fidx,1,early_late}),[1 sum(sig_channels{n_lfp,fidx,2,early_late})])' ...
                        repmat(find(sig_channels{n_lfp,fidx,2,early_late}),[1 sum(sig_channels{n_lfp,fidx,1,early_late})])'];
                    tmp_delays = -25:-1;
                    tmp_timeStep = 26:225;
                    
                    tmp_FIT_negd_M1receiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
%                     tmp_FIT_neg1_M1receiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_FIT_negd_DLSreceiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
%                     tmp_FIT_neg1_DLSreceiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    
                    if FIT_params.doShuff
                        tmp_FIT_negd_M1receiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
%                         tmp_FIT_neg1_M1receiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                        tmp_FIT_negd_DLSreceiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
%                         tmp_FIT_neg1_DLSreceiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
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
                        
                        pair_FIT_negd_M1receiver = zeros(length(tmp_timeStep),length(tmp_delays));
%                         pair_FIT_neg1_M1receiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        pair_FIT_negd_DLSreceiver = zeros(length(tmp_timeStep),length(tmp_delays));
%                         pair_FIT_neg1_DLSreceiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        
                        if FIT_params.doShuff
                            pair_FIT_negd_M1receiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
%                             pair_FIT_neg1_M1receiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                            pair_FIT_negd_DLSreceiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
%                             pair_FIT_neg1_DLSreceiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                        end
                        
                        tmpBinTimes = [];
                        
                        for timeStep = 1:length(tmp_timeStep)
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
                            for delay = 1:length(tmp_delays)
                                
                                %%% M1 receiver
                                
                                C = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                
                                % -d delay
                                Y = tmp_M1_LFP(tmp_timeStep(timeStep),:);
                                hY = tmp_M1_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                hX = tmp_DLS_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                C = repmat(C,size(Y,1),1);
                                C = C(:)';
                                Y = Y(:)';
                                hY = hY(:)';
                                hX = hX(:)';
                                
                                tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                pair_FIT_negd_M1receiver(timeStep,delay) = tmpFIT;
                                
                                % -1 delay
%                                 hY = tmp_M1_LFP(tmp_timeStep(timeStep)-1,:);
%                                 hY = hY(:)';
%                                 
%                                 tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
%                                 pair_FIT_neg1_M1receiver(timeStep,delay) = tmpFIT;
                                
                                %%% DLS receiver
                                
                                % -d delay
                                Y = tmp_DLS_LFP(tmp_timeStep(timeStep),:);
                                hY = tmp_DLS_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                hX = tmp_M1_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                Y = Y(:)';
                                hY = hY(:)';
                                hX = hX(:)';
                                
                                tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                pair_FIT_negd_DLSreceiver(timeStep,delay) = tmpFIT;
                                
                                % -1 delay
%                                 hY = tmp_DLS_LFP(tmp_timeStep(timeStep)-1,:);
%                                 hY = hY(:)';
%                                 
%                                 tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
%                                 pair_FIT_neg1_DLSreceiver(timeStep,delay) = tmpFIT;
                                
                                %%% FIT shuffle
                                
                                if FIT_params.doShuff
                                    for shidx = 1:FIT_params.nShuff
                                        
                                        %%% M1 receiver
                                        
                                        C = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                        C = C(shuffIdxs(:,shidx));
                                        
                                        % -d delay
                                        Y = tmp_M1_LFP(tmp_timeStep(timeStep),:);
                                        hY = tmp_M1_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                        hX = tmp_DLS_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                        C = repmat(C,size(Y,1),1);
                                        C = C(:)';
                                        Y = Y(:)';
                                        hY = hY(:)';
                                        hX = hX(:)';
                                        
                                        tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                        pair_FIT_negd_M1receiver_Sh(timeStep,delay,shidx) = tmpFIT;
                                        
                                        % -1 delay
%                                         hY = tmp_M1_LFP(tmp_timeStep(timeStep)-1,:);
%                                         hY = hY(:)';
%                                         
%                                         tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
%                                         pair_FIT_neg1_M1receiver_Sh(timeStep,delay,shidx) = tmpFIT;
                                        
                                        %%% DLS receiver
                                        
                                        % -d delay
                                        Y = tmp_DLS_LFP(tmp_timeStep(timeStep),:);
                                        hY = tmp_DLS_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                        hX = tmp_M1_LFP(tmp_timeStep(timeStep)+tmp_delays(delay),:);
                                        Y = Y(:)';
                                        hY = hY(:)';
                                        hX = hX(:)';
                                        
                                        tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                        pair_FIT_negd_DLSreceiver_Sh(timeStep,delay,shidx) = tmpFIT;
                                        
                                        % -1 delay
%                                         hY = tmp_DLS_LFP(tmp_timeStep(timeStep)-1,:);
%                                         hY = hY(:)';
%                                         
%                                         tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
%                                         pair_FIT_neg1_DLSreceiver_Sh(timeStep,delay,shidx) = tmpFIT;
                                        
                                    end
                                end
                            end
                        end
                                
                        tmp_FIT_negd_M1receiver(pairIdx,:,:) = pair_FIT_negd_M1receiver;
%                         tmp_FIT_neg1_M1receiver(pairIdx,:,:) = pair_FIT_neg1_M1receiver;
                        tmp_FIT_negd_DLSreceiver(pairIdx,:,:) = pair_FIT_negd_DLSreceiver;
%                         tmp_FIT_neg1_DLSreceiver(pairIdx,:,:) = pair_FIT_neg1_DLSreceiver;
                        
                        if FIT_params.doShuff
                            tmp_FIT_negd_M1receiver_Sh(pairIdx,:,:,:) = pair_FIT_negd_M1receiver_Sh;
%                             tmp_FIT_neg1_M1receiver_Sh(pairIdx,:,:,:) = pair_FIT_neg1_M1receiver_Sh;
                            tmp_FIT_negd_DLSreceiver_Sh(pairIdx,:,:,:) = pair_FIT_negd_DLSreceiver_Sh;
%                             tmp_FIT_neg1_DLSreceiver_Sh(pairIdx,:,:,:) = pair_FIT_neg1_DLSreceiver_Sh;
                        end
                        
                    end
                    
                    FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver = tmp_FIT_negd_M1receiver;
%                     FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_neg1_M1receiver = tmp_FIT_neg1_M1receiver;
                    FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver = tmp_FIT_negd_DLSreceiver;
%                     FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_neg1_DLSreceiver = tmp_FIT_neg1_DLSreceiver;
   
                    if FIT_params.doShuff
                        FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiverSh = tmp_FIT_negd_M1receiver_Sh;
%                         FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_neg1_M1receiverSh = tmp_FIT_neg1_M1receiver_Sh;
                        FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh = tmp_FIT_negd_DLSreceiver_Sh;
%                         FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_neg1_DLSreceiverSh = tmp_FIT_neg1_DLSreceiver_Sh;
                    end
                    
                end
            end
        end
        
        %% save results

        if FIT_params.saveFIT
            save([FIT_params.save_path 'FIT_' params.animals{animal} '_' date ...
                '_LFPearlylate_' num2str(FIT_params.doLFPearlylate) ...
                '_LFPday_' num2str(FIT_params.doLFPday) ...
                '_Spikesday_' num2str(FIT_params.doSpikesday) ...
                '.mat'],'-v7.3','FITout','params','FIT_params')
        end
        
    end
end

%% calculate LFP FIT by day

% if FIT_params.doLFPday==1
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
%                             tmp_sig = clusterStat(infQuant,infQuantSh,FIT_params.MIclusterparams(1),FIT_params.MIclusterparams(2));
%                             infQuant(~tmp_sig) = 0;
%                             [~,i] = max(infQuant);
%                             sigTime_start = find(MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(1));
%                             sigTime_end = find(MI_info.MIout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(2));
%                             sig_chan = [sig_chan i>=sigTime_start & i<=sigTime_end];
%                         end
%                         sig_channels{n_lfp,fidx,aidx,day} = sig_chan;
%                     end
%                 end
%             end
%         end
%         
%         tmpFITout = cell(length(params.reachFeatures));
%         parfor fidx = 1:length(params.reachFeatures)        
%             for day = 1:length(MI_info.MIout{animal}.LFPday)
%                 for n_lfp = 1:length(params.lfpFeatures)
%                     
%                     if FIT_params.doShuff
%                         shuffIdxs = zeros(size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).(params.areas{1}),3),FIT_params.nShuff);
%                         for shidx = 1:FIT_params.nShuff
%                             shuffIdxs(:,shidx) = randperm(size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).(params.areas{1}),3));
%                         end
%                     end
%                     
%                     pair_count = 1;
%                     for m1_chan = 1:size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).M1,1)
%                         for dls_chan = 1:size(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).DLS,1)
%                             if sig_channels{n_lfp,fidx,1,day}(m1_chan) && sig_channels{n_lfp,fidx,2,day}(dls_chan)
%                                 disp(['Computing: ' params.animals{animal} ' | ' params.lfpFeatures{n_lfp} ' | M1 channel ' num2str(m1_chan) ' & DLS channel ' num2str(dls_chan) ' | day ' num2str(day) ' of ' num2str(length(MI_info.MIout{animal}.LFPday)) ' | ' params.reachFeatures{fidx}])
%                                 tmp_M1_lfp = squeeze(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).M1(m1_chan,:,:));
%                                 tmp_DLS_lfp = squeeze(LFP{animal}{day}.(params.lfpFeatures{n_lfp}).DLS(dls_chan,:,:));
%                                 
%                                 tmpBinTimes = [];
%                                 sTime = 1;
%                                 timestep_count = 1;
%                                 for timeStep = 1:size(tmp_M1_lfp,1)/FIT_params.timeJump
%                                     
%                                     eTime = sTime + FIT_params.windowSize - 1;
%                                     if sTime-FIT_params.delay<1 || eTime+FIT_params.delay>size(tmp_M1_lfp,1)
%                                         sTime = sTime + FIT_params.timeJump;
%                                     else
%                                         delay_idx = 1;
%                                         for delay = -25:25
%                                             eTime = sTime + FIT_params.windowSize - 1;
%                                             X1 = tmp_M1_lfp(sTime:eTime,:);
%                                             X2 = tmp_DLS_lfp(sTime+delay:eTime+delay,:);
%                                             
%                                             Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
%                                             Y = repmat(Y,size(X1,1),1);
%                                             
%                                             X1 = X1(:);
%                                             X2 = X2(:);
%                                             Y = Y(:);
%                                             
%                                             I = FIT(Y',X1',X2',FIT_params.opts);
%                                             
%                                             tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).shared(pair_count,timestep_count,delay_idx) = I.shared;
%                                             tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_m1(pair_count,timestep_count,delay_idx) = I.uniqueX1;
%                                             tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_dls(pair_count,timestep_count,delay_idx) = I.uniqueX2;
%                                             tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).complementary(pair_count,timestep_count,delay_idx) = I.complementary;
%                                             
%                                             if FIT_params.doShuff
%                                                 
%                                                 tmp_shared = zeros(1,FIT_params.nShuff);
%                                                 tmp_unique_m1 = zeros(1,FIT_params.nShuff);
%                                                 tmp_unique_dls = zeros(1,FIT_params.nShuff);
%                                                 tmp_complementary = zeros(1,FIT_params.nShuff);
%                                                 
%                                                 for shidx = 1:FIT_params.nShuff
%                                                     Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(shuffIdxs(:,shidx));
%                                                     Y = repmat(Y,FIT_params.windowSize,1);
%                                                     Y = Y(:);
%                                                     
%                                                     I = FIT(Y',X1',X2',FIT_params.opts);
%                                                     
%                                                     tmp_shared(shidx) = I.shared;
%                                                     tmp_unique_m1(shidx) = I.uniqueX1;
%                                                     tmp_unique_dls(shidx) = I.uniqueX2;
%                                                     tmp_complementary(shidx) = I.complementary;
%                                                 end
%                                                 
%                                                 tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).sharedSh(pair_count,timestep_count,delay_idx,:) = tmp_shared;
%                                                 tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_m1Sh(pair_count,timestep_count,delay_idx,:) = tmp_unique_m1;
%                                                 tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).unique_dlsSh(pair_count,timestep_count,delay_idx,:) = tmp_unique_dls;
%                                                 tmpFITout{fidx}{day}.(params.lfpFeatures{n_lfp}).complementarySh(pair_count,timestep_count,delay_idx,:) = tmp_complementary;
%                                                 
%                                             end
%                                             
%                                             delay_idx = delay_idx + 1;
%                                         end
%                                         timestep_count = timestep_count + 1;
%                                         sTime = sTime + FIT_params.timeJump;
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
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).shared;
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1 = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_m1;
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dls = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_dls;
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementary = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).complementary;            
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).sharedSh;
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_m1Sh = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_m1Sh;
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).unique_dlsSh = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).unique_dlsSh;
%                     FITout{animal}.LFPday{day}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).complementarySh = tmpFITout{fidx}{early_late}.(params.lfpFeatures{n_lfp}).complementarySh;                
%                 end
%             end
%         end
%       
%     end
% end

%% calculate Spikes FIT by day

% if FIT_params.doSpikesday==1
%     for animal = 1:numel(params.animals)
%         
%         % determine channels with significant mutual information
%         sig_units = cell(length(params.reachFeatures),length(params.areas),length(MI_info.MIout{animal}.spikesDay));
%         for fidx = 1:length(params.reachFeatures)
%             for aidx = 1:length(params.areas)
%                 for day = 1:length(MI_info.MIout{animal}.spikesDay)
%                     sig_unit = [];
%                     for unit = 1:size(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                         infQuant = MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
%                         infQuantSh = MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
%                         tmp_sig = clusterStat(infQuant,infQuantSh,FIT_params.MIclusterparams(1),FIT_params.MIclusterparams(2));
%                         infQuant(~tmp_sig) = 0;
%                         [~,i] = max(infQuant);
%                         sigTime_start = find(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(1));
%                         sigTime_end = find(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(2));
%                         sig_unit = [sig_unit i>=sigTime_start & i<=sigTime_end];
%                     end
%                     sig_units{fidx,aidx,day} = sig_unit;
%                 end
%             end
%         end
%         
%         tmpFITout = cell(1,length(params.reachFeatures));
%         
%         parfor fidx = 1:length(params.reachFeatures)
%             for day = 1:length(MI_info.MIout{animal}.spikesDay)
%                 if FIT_params.doShuff
%                     shuffIdxs = zeros(size(spikes{animal}{day}.(params.areas{1}),3),FIT_params.nShuff);
%                     for shidx = 1:FIT_params.nShuff
%                         shuffIdxs(:,shidx) = randperm(size(spikes{animal}{day}.(params.areas{1}),3));
%                     end
%                 end
%                 pair_count = 1;
%                 for m1_unit = 1:size(MI_info.MIout{animal}.spikesDay{day}.M1.(params.reachFeatures{fidx}).info,1)
%                     for dls_unit = 1:size(MI_info.MIout{animal}.spikesDay{day}.DLS.(params.reachFeatures{fidx}).info,1)
%                         if sig_units{fidx,1,day}(m1_unit) && sig_units{fidx,2,day}(dls_unit)
%                             
%                             disp(['Computing... ' params.animals{animal} ' | M1 unit ' num2str(m1_unit) ' & DLS unit ' num2str(dls_unit) ' | day ' num2str(day) ' | ' params.reachFeatures{fidx}])
%                             
%                             tmp_M1_spiking = squeeze(spikes{animal}{day}.M1(m1_unit,:,:));
%                             tmp_DLS_spiking = squeeze(spikes{animal}{day}.DLS(dls_unit,:,:));
%                             
%                             nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
%                             tmp_M1_spiking = tmp_M1_spiking(:,nonNaNtrials);
%                             tmp_DLS_spiking = tmp_DLS_spiking(:,nonNaNtrials);
%                             
%                             tmpBinTimes = [];
%                             sTime = 1;
%                             timestep_count = 1;
%                             for timeStep = 1:size(tmp_M1_spiking,1)/FIT_params.timeJump
%                                 
%                                 eTime = sTime + FIT_params.windowSize - 1;
%                                 if sTime-FIT_params.delay<1 || eTime+FIT_params.delay>size(tmp_M1_spiking,1)
%                                     sTime = sTime + FIT_params.timeJump;
%                                 else
%                                     delay_idx = 1;
%                                     for delay = -25:25
%                                         eTime = sTime + FIT_params.windowSize - 1;
%                                         
%                                         Y = tmp_M1_spiking(sTime:eTime,:);
%                                         hY = tmp_M1_spiking(sTime+delay:eTime+delay,:);
%                                         hX = tmp_M1_spiking(sTime+delay:eTime+delay,:);
%                                         
%                                         C = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
%                                         C = repmat(C,size(Y,1),1);
%                                         
%                                         Y = Y(:)';
%                                         hY = hY(:)';
%                                         hX = hX(:)';
%                                         C = C(:)';
%                                         
%                                         FIT_biased = FIT_choice(C,hX,hY,Y,FIT_params.opts);
%                                         
%                                         tmpFITout{fidx}{day}.FIT(pair_count,timestep_count,delay_idx) = FIT_biased;
%                                         tmpBinTimes = [tmpBinTimes (((sTime-1)*binSize)+tMin)-5001];
%                                         %
%                                         %                                         if FIT_params.doShuff
%                                         %
%                                         %                                             tmp_shared = zeros(1,FIT_params.nShuff);
%                                         %                                             tmp_unique_m1 = zeros(1,FIT_params.nShuff);
%                                         %                                             tmp_unique_dls = zeros(1,FIT_params.nShuff);
%                                         %                                             tmp_complementary = zeros(1,FIT_params.nShuff);
%                                         %
%                                         %                                             for shidx = 1:FIT_params.nShuff
%                                         %                                                 Y = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(shuffIdxs(:,shidx));
%                                         %                                                 Y = repmat(Y,FIT_params.windowSize,1);
%                                         %                                                 Y = Y(:);
%                                         %
%                                         %                                                 I = FIT(Y',X1',X2',FIT_params.opts);
%                                         %
%                                         %                                                 tmp_shared(shidx) = I.shared;
%                                         %                                                 tmp_unique_m1(shidx) = I.uniqueX1;
%                                         %                                                 tmp_unique_dls(shidx) = I.uniqueX2;
%                                         %                                                 tmp_complementary(shidx) = I.complementary;
%                                         %                                             end
%                                         %
%                                         %                                             tmpFITout{fidx}{day}.sharedSh(pair_count,timestep_count,delay_idx,:) = tmp_shared;
%                                         %                                             tmpFITout{fidx}{day}.unique_m1Sh(pair_count,timestep_count,delay_idx,:) = tmp_unique_m1;
%                                         %                                             tmpFITout{fidx}{day}.unique_dlsSh(pair_count,timestep_count,delay_idx,:) = tmp_unique_dls;
%                                         %                                             tmpFITout{fidx}{day}.complementarySh(pair_count,timestep_count,delay_idx,:) = tmp_complementary;
%                                         %
%                                         %                                         end
%                                         
%                                         delay_idx = delay_idx + 1;
%                                     end
%                                     timestep_count = timestep_count + 1;
%                                     sTime = sTime + FIT_params.timeJump;
%                                 end
%                             end
%                             
%                         end
%                         pair_count = pair_count + 1;
%                     end
%                 end
%             end
%         end
%         
%         for fidx = 1:length(params.reachFeatures)
%             if ~isempty(tmpFITout{fidx})
%                 for day = 1:length(MI_info.MIout{animal}.spikesDay)
%                     if day<=length(tmpFITout{fidx})
%                         if ~isempty(tmpFITout{fidx}{day})
%                             FITout{animal}.spikesDay{day}.(params.reachFeatures{fidx}).FIT = tmpFITout{fidx}{day}.FIT;
%                         end
%                     end
%                 end
%             end
%         end
%         
%     end
% end

%% save results

% if FIT_params.saveFIT
%     save([FIT_params.save_path 'FIT_' date ...
%         '_LFPearlylate_' num2str(FIT_params.doLFPearlylate) ...
%         '_LFPday_' num2str(FIT_params.doLFPday) ...
%         '_Spikesday_' num2str(FIT_params.doSpikesday) ...
%         '.mat'],'-v7.3','FITout','params','FIT_params')
% end

end

