function [FITout] = compute_FIT_Pooled(LFP, spikes, spikes_ave, reachFeatures, LFP_early_late, spikes_early_late, reachFeatures_early_late, params, FIT_params, tMin, binSize,animals_to_run)

% load MI results
MI_info = load(FIT_params.path_to_MI);
    
%% calculate LFP FIT by early/late
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
                            infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                            tmp_sig = clusterStat(infQuant,infQuantSh,FIT_params.MIclusterparams(1),FIT_params.MIclusterparams(2));
                            infQuant(~tmp_sig) = 0;
                            [~,i] = max(infQuant);
                            sigTime_start = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(1));
                            sigTime_end = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(2));
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
                    tmp_timeStep = 1:floor((length(26:225)-(FIT_params.windowSize-FIT_params.timeJump))/FIT_params.timeJump);
                    
                    tmp_FIT_negd_M1receiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_FIT_negd_DLSreceiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    if FIT_params.doShuff
                        tmp_FIT_negd_M1receiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                        tmp_FIT_negd_DLSreceiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
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
                        pair_FIT_negd_DLSreceiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        
                        if FIT_params.doShuff
                            pair_FIT_negd_M1receiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                            pair_FIT_negd_DLSreceiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                        end
                        
                        tmpBinTimes = [];
                        sTime = 26;

                        for timeStep = 1:length(tmp_timeStep)
                            
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
                            delayIdx = 1;

                            for delay = tmp_delays
                                
                                eTime = sTime + FIT_params.windowSize-1;

                                %%% M1 receiver
                                
                                C = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                
                                % -d delay
                                Y = tmp_M1_LFP(sTime:eTime,:);
                                hY = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                hX = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                C = repmat(C,size(Y,1),1);
                                C = C(:)';
                                Y = Y(:)';
                                hY = hY(:)';
                                hX = hX(:)';
                                
                                tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                pair_FIT_negd_M1receiver(timeStep,delayIdx) = tmpFIT;

                                %%% DLS receiver
                                
                                % -d delay
                                Y = tmp_DLS_LFP(sTime:eTime,:);
                                hY = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                hX = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                Y = Y(:)';
                                hY = hY(:)';
                                hX = hX(:)';
                                
                                tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                pair_FIT_negd_DLSreceiver(timeStep,delayIdx) = tmpFIT;

                                %%% FIT shuffle
                                
                                if FIT_params.doShuff
                                    for shidx = 1:FIT_params.nShuff
                                        
                                        %%% M1 receiver
                                        
                                        C = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})(nonNaNtrials);
                                        C = C(shuffIdxs(:,shidx));
                                        
                                        % -d delay
                                        Y = tmp_M1_LFP(sTime:eTime,:);
                                        hY = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                        hX = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                        C = repmat(C,size(Y,1),1);
                                        C = C(:)';
                                        Y = Y(:)';
                                        hY = hY(:)';
                                        hX = hX(:)';
                                        
                                        tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                        pair_FIT_negd_M1receiver_Sh(timeStep,delayIdx,shidx) = tmpFIT;
                                        
                                        %%% DLS receiver
                                        
                                        % -d delay
                                        Y = tmp_DLS_LFP(sTime:eTime,:);
                                        hY = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                        hX = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                        Y = Y(:)';
                                        hY = hY(:)';
                                        hX = hX(:)';
                                        
                                        tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                        pair_FIT_negd_DLSreceiver_Sh(timeStep,delayIdx,shidx) = tmpFIT;
                                        
                                    end
                                end
                                delayIdx = delayIdx + 1;
                            end
                            sTime = sTime + FIT_params.timeJump;
                        end
                                
                        tmp_FIT_negd_M1receiver(pairIdx,:,:) = pair_FIT_negd_M1receiver;
                        tmp_FIT_negd_DLSreceiver(pairIdx,:,:) = pair_FIT_negd_DLSreceiver;
                        
                        if FIT_params.doShuff
                            tmp_FIT_negd_M1receiver_Sh(pairIdx,:,:,:) = pair_FIT_negd_M1receiver_Sh;
                            tmp_FIT_negd_DLSreceiver_Sh(pairIdx,:,:,:) = pair_FIT_negd_DLSreceiver_Sh;
                        end
                        
                    end
                    
                    FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver = tmp_FIT_negd_M1receiver;
                    FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver = tmp_FIT_negd_DLSreceiver;
   
                    if FIT_params.doShuff
                        FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiverSh = tmp_FIT_negd_M1receiver_Sh;
                        FITout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh = tmp_FIT_negd_DLSreceiver_Sh;
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
                '_Pooled.mat'],'-v7.3','FITout','params','FIT_params')
        end
        
    end
end

%% calculate Spikes FIT by day
if FIT_params.doSpikesday==1
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
                        tmp_sig = clusterStat(infQuant,infQuantSh,FIT_params.MIclusterparams(1),FIT_params.MIclusterparams(2));
                        infQuant(~tmp_sig) = 0;
                        [~,i] = max(infQuant);
                        sigTime_start = find(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(1));
                        sigTime_end = find(MI_info.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==FIT_params.MIsigtiming(2));
                        sig_unit = [sig_unit i>=sigTime_start & i<=sigTime_end];
                    end
                    sig_units{fidx,aidx,day} = sig_unit;
                end
            end
        end
        
        for fidx = 1:length(params.reachFeatures)    
        
            tmpFITout = cell(1,length(MI_info.MIout{animal}.spikesDay));

            for day = 1:length(MI_info.MIout{animal}.spikesDay)

                nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
                if FIT_params.doShuff
                    shuffIdxs = zeros(length(nonNaNtrials),FIT_params.nShuff);
                    for shidx = 1:FIT_params.nShuff
                        shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
                    end
                end

                tmp_pairs = [...
                    repmat(find(sig_units{fidx,1,day}),[1 sum(sig_units{fidx,2,day})])' ...
                    repmat(find(sig_units{fidx,2,day}),[1 sum(sig_units{fidx,1,day})])'];
                tmp_delays = -25:-1;
                tmp_timeStep = 1:floor((length(26:225)-(FIT_params.windowSize-FIT_params.timeJump))/FIT_params.timeJump);
                    
                tmp_FIT_negd_M1receiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                tmp_FIT_negd_DLSreceiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                if FIT_params.doShuff
                    tmp_FIT_negd_M1receiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                    tmp_FIT_negd_DLSreceiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                end
                
                parfor pairIdx = 1:size(tmp_pairs,1)

                    m1_unit = tmp_pairs(pairIdx,1);
                    dls_unit = tmp_pairs(pairIdx,2);
                    disp(['Computing: ' params.animals{animal} ' | M1 unit ' num2str(m1_unit) ' & DLS unit ' num2str(dls_unit) ' | ' num2str(day) ' day | ' params.reachFeatures{fidx}])

                    tmp_M1_spiking = squeeze(spikes{animal}{day}.M1(m1_unit,:,:));
                    tmp_DLS_spiking = squeeze(spikes{animal}{day}.DLS(dls_unit,:,:));
                    
                    nonNaNtrials = find(~isnan(reachFeatures{animal}{day}.(params.reachFeatures{fidx})));
                    tmp_M1_spiking = tmp_M1_spiking(:,nonNaNtrials);
                    tmp_DLS_spiking = tmp_DLS_spiking(:,nonNaNtrials);

                    pair_FIT_negd_M1receiver = zeros(length(tmp_timeStep),length(tmp_delays));
                    pair_FIT_negd_DLSreceiver = zeros(length(tmp_timeStep),length(tmp_delays));

                    if FIT_params.doShuff
                        pair_FIT_negd_M1receiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                        pair_FIT_negd_DLSreceiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),FIT_params.nShuff);
                    end

                    tmpBinTimes = [];
                    sTime = 26;

                    for timeStep = 1:length(tmp_timeStep)
                            
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
                            delayIdx = 1;

                            for delay = tmp_delays
                                
                                eTime = sTime + FIT_params.windowSize-1;

                                %%% M1 receiver
                                
                                C = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
                                
                                % -d delay
                                Y = tmp_M1_spiking(sTime:eTime,:);
                                hY = tmp_M1_spiking(sTime+delay:eTime+delay,:);
                                hX = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
                                C = repmat(C,size(Y,1),1);
                                C = C(:)';
                                Y = Y(:)';
                                hY = hY(:)';
                                hX = hX(:)';
                                
                                tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                pair_FIT_negd_M1receiver(timeStep,delayIdx) = tmpFIT;

                                %%% DLS receiver
                                
                                % -d delay
                                Y = tmp_DLS_spiking(sTime:eTime,:);
                                hY = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
                                hX = tmp_M1_spiking(sTime+delay:eTime+delay,:);
                                Y = Y(:)';
                                hY = hY(:)';
                                hX = hX(:)';
                                
                                tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                pair_FIT_negd_DLSreceiver(timeStep,delayIdx) = tmpFIT;

                                %%% FIT shuffle
                                
                                if FIT_params.doShuff
                                    for shidx = 1:FIT_params.nShuff
                                        
                                        %%% M1 receiver
                                        
                                        C = reachFeatures{animal}{day}.(params.reachFeatures{fidx})(nonNaNtrials);
                                        C = C(shuffIdxs(:,shidx));
                                        
                                        % -d delay
                                        Y = tmp_M1_spiking(sTime:eTime,:);
                                        hY = tmp_M1_spiking(sTime+delay:eTime+delay,:);
                                        hX = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
                                        C = repmat(C,size(Y,1),1);
                                        C = C(:)';
                                        Y = Y(:)';
                                        hY = hY(:)';
                                        hX = hX(:)';
                                        
                                        tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                        pair_FIT_negd_M1receiver_Sh(timeStep,delayIdx,shidx) = tmpFIT;
                                        
                                        %%% DLS receiver
                                        
                                        % -d delay
                                        Y = tmp_DLS_spiking(sTime:eTime,:);
                                        hY = tmp_DLS_spiking(sTime+delay:eTime+delay,:);
                                        hX = tmp_M1_spiking(sTime+delay:eTime+delay,:);
                                        Y = Y(:)';
                                        hY = hY(:)';
                                        hX = hX(:)';
                                        
                                        tmpFIT = FIT_choice(C,hX,hY,Y,FIT_params.opts);
                                        pair_FIT_negd_DLSreceiver_Sh(timeStep,delayIdx,shidx) = tmpFIT;
                                        
                                    end
                                end
                                delayIdx = delayIdx + 1;
                            end
                            sTime = sTime + FIT_params.timeJump;
                        end
                                
                        tmp_FIT_negd_M1receiver(pairIdx,:,:) = pair_FIT_negd_M1receiver;
                        tmp_FIT_negd_DLSreceiver(pairIdx,:,:) = pair_FIT_negd_DLSreceiver;
                        
                        if FIT_params.doShuff
                            tmp_FIT_negd_M1receiver_Sh(pairIdx,:,:,:) = pair_FIT_negd_M1receiver_Sh;
                            tmp_FIT_negd_DLSreceiver_Sh(pairIdx,:,:,:) = pair_FIT_negd_DLSreceiver_Sh;
                        end
                        
                    end
                    
                    FITout{animal}.spikesDay{day}.(params.reachFeatures{fidx}).FIT_negd_M1receiver = tmp_FIT_negd_M1receiver;
                    FITout{animal}.spikesDay{day}.(params.reachFeatures{fidx}).FIT_negd_DLSreceiver = tmp_FIT_negd_DLSreceiver;
   
                    if FIT_params.doShuff
                        FITout{animal}.spikesDay{day}.(params.reachFeatures{fidx}).FIT_negd_M1receiverSh = tmp_FIT_negd_M1receiver_Sh;
                        FITout{animal}.spikesDay{day}.(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh = tmp_FIT_negd_DLSreceiver_Sh;
                    end
                    
            end
            
            % save results
            if FIT_params.saveFIT
                save([FIT_params.save_path 'FIT_' params.animals{animal} '_' date ...
                    '_LFPearlylate_' num2str(FIT_params.doLFPearlylate) ...
                    '_LFPday_' num2str(FIT_params.doLFPday) ...
                    '_Spikesday_' num2str(FIT_params.doSpikesday) ...
                    '_' params.reachFeatures{fidx} '.mat'],'-v7.3','FITout','params','FIT_params')
            end
            
        end
    end
end
  
end

