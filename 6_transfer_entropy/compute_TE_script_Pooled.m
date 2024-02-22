function [TEout] = compute_TE_script_Pooled(LFP, spikes, reachFeatures, LFP_early_late, reachFeatures_early_late, params, TE_params, tMin, binSize,animals_to_run)

% load MI results    
MI_info = load(TE_params.path_to_MI);

% calculate LFP DI by early/late
if TE_params.doLFPearlylate==1
    for animal = animals_to_run  
        clearvars DIout

        for fidx = 1:length(params.reachFeatures)
            for early_late = 1:2
                for n_lfp = 1:length(params.lfpFeatures)

                    nonNaNtrials = find(~isnan(reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx})));
                    
                    if TE_params.doShuff
                        shuffIdxs = zeros(length(nonNaNtrials),TE_params.nShuff);
                        for shidx = 1:TE_params.nShuff
                            shuffIdxs(:,shidx) = randperm(length(nonNaNtrials));
                        end
                    end

                    tmp_pairs = [...
                        reshape(repmat(1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1,1), ...
                        [size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS,1) 1]), ...
                        [1 size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1,1)*size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS,1)])' ...
                        repmat(1:size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).DLS,1),[1 size(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).M1,1)])'];
                    tmp_delays = -25:-1;
                    tmp_timeStep = 1:floor((length(26:225)-(TE_params.windowSize-TE_params.timeJump))/TE_params.timeJump);
                    
                    tmp_TE_negd_M1receiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    tmp_TE_negd_DLSreceiver = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays));
                    if TE_params.doShuff
                        tmp_TE_negd_M1receiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),TE_params.nShuff);
                        tmp_TE_negd_DLSreceiver_Sh = zeros(size(tmp_pairs,1),length(tmp_timeStep),length(tmp_delays),TE_params.nShuff);
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
                        
                        pair_TE_negd_M1receiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        pair_TE_negd_DLSreceiver = zeros(length(tmp_timeStep),length(tmp_delays));
                        
                        if TE_params.doShuff
                            pair_TE_negd_M1receiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),TE_params.nShuff);
                            pair_TE_negd_DLSreceiver_Sh = zeros(length(tmp_timeStep),length(tmp_delays),TE_params.nShuff);
                        end
                        
                        tmpBinTimes = [];
                        sTime = 26;

                        for timeStep = 1:length(tmp_timeStep)
                            
                            tmpBinTimes = [tmpBinTimes (((tmp_timeStep(timeStep)-1)*binSize)+tMin)-5001];
                            delayIdx = 1;

                            for delay = tmp_delays
                                
                                eTime = sTime + TE_params.windowSize-1;

                                %%% M1 receiver
                                
                                    % -d delay
                                    Y = tmp_M1_LFP(sTime:eTime,:);
                                    hY = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                    hX = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                    Y = Y(:)';
                                    hY = hY(:)';
                                    hX = hX(:)';
                                
                                    TEtmp = compute_TE(hX,hY,Y,TE_params.opts);
                                    pair_TE_negd_M1receiver(timeStep,delayIdx) = TEtmp;
                                
                                %%% DLS receiver
                                
                                    % -d delay
                                    Y = tmp_DLS_LFP(sTime:eTime,:);
                                    hY = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                    hX = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                    Y = Y(:)';
                                    hY = hY(:)';
                                    hX = hX(:)';
                                    
                                    TEtmp = compute_TE(hX,hY,Y,TE_params.opts);
                                    pair_TE_negd_DLSreceiver(timeStep,delayIdx) = TEtmp;

                                if TE_params.doShuff
                                    for shidx = 1:TE_params.nShuff
                                        
                                        %%% M1 receiver
                                        
                                            % -d delay
                                            hX = tmp_DLS_LFP(sTime+delay:eTime+delay,:);
                                            hX = hX(:,shuffIdxs(:,shidx));
                                            Y = tmp_M1_LFP(sTime:eTime,:);
                                            hY = tmp_M1_LFP(sTime+delay:eTime+delay,:);

                                            Y = Y(:)';
                                            hY = hY(:)';
                                            hX = hX(:)';

                                            TEtmp = compute_TE(hX,hY,Y,TE_params.opts);
                                            pair_TE_negd_M1receiver_Sh(timeStep,delayIdx,shidx) = TEtmp;

                                        %%% DLS receiver
                                        
                                            % -d delay
                                            hX = tmp_M1_LFP(sTime+delay:eTime+delay,:);
                                            hX = hX(:,shuffIdxs(:,shidx));
                                            Y = tmp_DLS_LFP(sTime:eTime,:);
                                            hY = tmp_DLS_LFP(sTime+delay:eTime+delay,:);

                                            Y = Y(:)';
                                            hY = hY(:)';
                                            hX = hX(:)';

                                            TEtmp = compute_TE(hX,hY,Y,TE_params.opts);
                                            pair_TE_negd_DLSreceiver_Sh(timeStep,delayIdx,shidx) = TEtmp;
                                        
                                    end
                                end
                                delayIdx = delayIdx + 1;
                            end
                            sTime = sTime + TE_params.timeJump;
                        end
                                
                        tmp_TE_negd_M1receiver(pairIdx,:,:) = pair_TE_negd_M1receiver;
                        tmp_TE_negd_DLSreceiver(pairIdx,:,:) = pair_TE_negd_DLSreceiver;
                        
                        if TE_params.doShuff
                            tmp_TE_negd_M1receiver_Sh(pairIdx,:,:,:) = pair_TE_negd_M1receiver_Sh;
                            tmp_TE_negd_DLSreceiver_Sh(pairIdx,:,:,:) = pair_TE_negd_DLSreceiver_Sh;
                        end
                        
                    end      
                                    
                    TEout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).TE_negd_M1receiver = tmp_TE_negd_M1receiver;
                    TEout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).TE_negd_DLSreceiver = tmp_TE_negd_DLSreceiver;
                
                    if TE_params.doShuff
                        TEout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).TE_negd_M1receiverSh = tmp_TE_negd_M1receiver_Sh;
                        TEout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).TE_negd_DLSreceiverSh = tmp_TE_negd_DLSreceiver_Sh;
                    end
                end
            end
        end
        
        % save results
        if TE_params.saveTE
            save([TE_params.save_path 'TE_' params.animals{animal} '_' date ...
                '_LFPearlylate_' num2str(TE_params.doLFPearlylate) ...
                '_LFPday_' num2str(TE_params.doLFPday) ...
                '_Spikesday_' num2str(TE_params.doSpikesday) ...
                '_Pooled.mat'],'-v7.3','TEout','TE_params','params')
        end
        
    end
end

end

function [info] = compute_TE(hX, hY, Y, opts)

%%% Script to compute Directed Information (DI) from the emitter X to the receiver Y

%%% - *hX*: must be an array of *nDimsX X nTrials* response matrix describing the past of the emitter variables on each of the *nDims* dimensions for each trial.
%%% - *hY*: must be an array of *nDimsY X nTrials* response matrix describing the past of the receiver variable on each of the *nDims* dimensions for each trial.
%%% - *Y*: must be an array of *nDimsY X nTrials* response matrix describing the response of each of the *nDims* dimensions for each trial.
%%% - *opts*: options used to calculate II (see further notes).

I1 = information(hY,Y,opts,{'I'});
I2 = information([hY; hX],Y,opts,{'I'});

info = I2{1}(1)-I1{1}(1);

end