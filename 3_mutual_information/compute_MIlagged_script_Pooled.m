function [lagMIout] = compute_MIlagged_script_Pooled(LFP, spikes, reachFeatures, LFP_early_late, reachFeatures_early_late, params, DI_params, tMin, binSize,animals_to_run)
% We use same parameters as for DI 
% load MI results    
MI_info = load(DI_params.path_to_MI);

% calculate LFP DI by early/late
if DI_params.doLFPearlylate==1
    for animal = animals_to_run  
        clearvars DIout

        for fidx = 1:length(params.reachFeatures)
            for early_late = 1:2
                for n_lfp = 1:length(params.lfpFeatures)

                    nonNaNtrials = find(~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})));
                    
                    if DI_params.doShuff
                        shuffIdxs = zeros(length(nonNaNtrials),DI_params.nShuff);
                        for shidx = 1:DI_params.nShuff
                            shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
                        end
                    end

                    tmp_pairs = [...
                        reshape(repmat(1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1,1), ...
                        [size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS,1) 1]), ...
                        [1 size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1,1)*size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS,1)])' ...
                        repmat(1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS,1),[1 size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1,1)])'];
                    tmp_delays = -25:-1;
                    tmp_timeStep = 1:floor((length(26:225)-(DI_params.windowSize-DI_params.timeJump))/DI_params.timeJump);
                    
                    tmp_lagMI_negd_M1receiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_lagMI_negd_DLSreceiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    if DI_params.doShuff
                        tmp_lagMI_negd_M1receiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),DI_params.nShuff);
                        tmp_lagMI_negd_DLSreceiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),DI_params.nShuff);
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
                        
                        pair_lagMI_negd_M1receiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        pair_lagMI_negd_DLSreceiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        
                        if DI_params.doShuff
                            pair_lagMI_negd_M1receiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),DI_params.nShuff);
                            pair_lagMI_negd_DLSreceiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),DI_params.nShuff);
                        end
                        
                        tmpBinTimes = [];
                        sTime = 26;

                        for timeStep = 1:length(tmp_timeStep)
                            
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
                            delayIdx = 1;

                            for delay = tmp_delays
                                
                                eTime = sTime + DI_params.windowSize-1;

                                %%% M1 receiver
                                
                                    % -d delay
                                    Y = tmp_M1_LFP(sTime:eTime,:);
                                    hY = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                    hX = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                    Y = Y(:)';
                                    hY = hY(:)';
                                    hX = hX(:)';
                                
                                    lagMItmp = cell2mat(information(hX,Y,DI_params.opts,{'I'}));
                                    pair_lagMI_negd_M1receiver(timeStep,delayIdx) = lagMItmp;
                                
                                %%% DLS receiver
                                
                                    % -d delay
                                    Y = tmp_DLS_LFP(sTime:eTime,:);
                                    hY = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                    hX = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                    Y = Y(:)';
                                    hY = hY(:)';
                                    hX = hX(:)';
                                    
                                    lagMItmp = cell2mat(information(hX,Y,DI_params.opts,{'I'}));
                                    pair_lagMI_negd_DLSreceiver(timeStep,delayIdx) = lagMItmp;

                                if DI_params.doShuff
                                    for shidx = 1:DI_params.nShuff
                                        
                                        %%% M1 receiver
                                        
                                            % -d delay
                                            hX = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                            hX = hX(:,shuffIdxs(:,shidx));
                                            Y = tmp_M1_LFP(sTime:eTime,:);
                                            hY = tmp_M1_LFP(sTime+delay:eTime+delay,:);

                                            Y = Y(:)';
                                            hY = hY(:)';
                                            hX = hX(:)';

                                            lagMItmp = cell2mat(information(hX,Y,DI_params.opts,{'I'}));
                                            pair_lagMI_negd_M1receiver_Sh(timeStep,delayIdx,shidx) = lagMItmp;

                                        %%% DLS receiver
                                        
                                            % -d delay
                                            hX = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                            hX = hX(:,shuffIdxs(:,shidx));
                                            Y = tmp_DLS_LFP(sTime:eTime,:);
                                            hY = tmp_DLS_LFP(sTime+delay:eTime+delay,:);

                                            Y = Y(:)';
                                            hY = hY(:)';
                                            hX = hX(:)';

                                            lagMItmp = cell2mat(information(hX,Y,DI_params.opts,{'I'}));
                                            pair_lagMI_negd_DLSreceiver_Sh(timeStep,delayIdx,shidx) = lagMItmp;
                                        
                                    end
                                end
                                delayIdx = delayIdx + 1;
                            end
                            sTime = sTime + DI_params.timeJump;
                        end
                                
                        tmp_lagMI_negd_M1receiver(pairIdx,:,:) = pair_lagMI_negd_M1receiver;
                        tmp_lagMI_negd_DLSreceiver(pairIdx,:,:) = pair_lagMI_negd_DLSreceiver;
                        
                        if DI_params.doShuff
                            tmp_lagMI_negd_M1receiver_Sh(pairIdx,:,:,:) = pair_lagMI_negd_M1receiver_Sh;
                            tmp_lagMI_negd_DLSreceiver_Sh(pairIdx,:,:,:) = pair_lagMI_negd_DLSreceiver_Sh;
                        end
                        
                    end      
                                    
                    lagMIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).DI_negd_M1receiver = tmp_lagMI_negd_M1receiver;
                    lagMIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).DI_negd_DLSreceiver = tmp_lagMI_negd_DLSreceiver;
                
                    if DI_params.doShuff
                        lagMIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).DI_negd_M1receiverSh = tmp_lagMI_negd_M1receiver_Sh;
                        lagMIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).DI_negd_DLSreceiverSh = tmp_lagMI_negd_DLSreceiver_Sh;
                    end
                end
            end
        end
        
        % save results
        if DI_params.saveDI
            save([DI_params.save_path 'MIlag_' params.animals{animal} '_' date ...
                '_LFPearlylate_' num2str(DI_params.doLFPearlylate) ...
                '_LFPday_' num2str(DI_params.doLFPday) ...
                '_Spikesday_' num2str(DI_params.doSpikesday) ...
                '_Pooled.mat'],'-v7.3','lagMIout','DI_params','params')
        end
        
    end
end

end
