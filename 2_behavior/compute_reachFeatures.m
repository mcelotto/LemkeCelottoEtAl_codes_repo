function reach_features = compute_reachFeatures_autcorr(rf_params, params, bin_params)
%
% computes and outputs reach_features struct with the following reach features:
%   -> maxPosX = maximum x position of paw during rf_params.prevBins before pellet touch
%   -> maxPosY = maximum y position of paw during rf_params.prevBins before pellet touch
%   -> maxVelX = maximum x velocity of paw during rf_params.prevBins before pellet touch
%   -> maxVelY = maximum y velocity of paw during rf_params.prevBins before pellet touch
%   -> maxAccX = maximum x acceleration of paw during rf_params.prevBins before pellet touch
%   -> maxAccY = maximum y acceleration of paw during rf_params.prevBins before pellet touch
%   -> minAccX = minimum x acceleration of paw during rf_params.prevBins before pellet touch
%   -> minAccY = minimum y acceleration of paw during rf_params.prevBins before pellet touch
%   -> maxVel = maximum velocity (mean x and y) of paw during rf_params.prevBins before pellet touch
%   -> maxAcc = maximum acceleration (mean x and y) of paw during rf_params.prevBins before pellet touch
%   -> minAcc = minimum acceleration (mean x and y) of paw during rf_params.prevBins before pellet touch
%   -> movDur = duration of movement from movement onset to pellet touch
%   -> distTrav = distance traveled by paw in x and y dimension from movement onset to pellet touch
%   -> PCAstab = PCA 'stability' metric
%
% notes:
%   uses "fillmissing" nan-interpolation with linear interpolation
%   and a moving window of 3

% load trajectory data
trajectories = load(rf_params.full_path);
tMin = trajectories.tMin;
tMax = trajectories.tMax;
trajectories = trajectories.trajectories;

% determine number of bins
numBins = length(tMin:tMax)/rf_params.binSize;
if rem(length(tMin:tMax),rf_params.binSize)~=0
    error('trajectory length is not divisible by bin size')
end

for animal = 1:numel(trajectories)
    for day = 1:numel(trajectories{animal})
        
        % load success rate
        successTmp = load([rf_params.data_path '/animal_data/' params.animals{animal} '/success_rate/success_rate_day_' num2str(params.days{animal}(day)) '.mat']);

        % intialize features
        trial_cons_nans = zeros(1,size(trajectories{animal}{day},2));
        trial_rf = zeros(13,size(trajectories{animal}{day},2));
        trial_rf_timing = zeros(11,size(trajectories{animal}{day},2)); % argmax/argmin of static reach features (no mov. onset and traj. length)
        trial_success = zeros(1,size(trajectories{animal}{day},2));
        tPoints = 1+size(trajectories{animal}{day},3)+(5000-tMax)-rf_params.binSize*rf_params.prevBins:size(trajectories{animal}{day},3)+(5000-tMax); % time points used for NaN trial disqualification and max/min features
        % instantaneous features
        trial_rf_inst = zeros(9, size(trajectories{animal}{day},2), size(trajectories{animal}{day},3)/bin_params.oldBin);
        autcorr_rf = zeros(9, size(trajectories{animal}{day},2), 2*size(trial_rf_inst,3)-1);
        
        for trial = 1:size(trajectories{animal}{day},2)
            
            % calculates consecutive NaNs during the specified number of time bins that precede pellet touch
            [~,trial_cons_nans(trial)] = nanCount(squeeze(trajectories{animal}{day}(1,trial,tPoints)));
            
            if trial_cons_nans(trial)<rf_params.nan_thresh
                
                trial_success(trial) = successTmp.trial_success(trial);
                
                for featureIdx = 1:6 % maxPosX, maxPosY, maxVelX, maxVelY, maxAccX, maxAccY
                    [trial_rf(featureIdx,trial), trial_rf_timing(featureIdx,trial)] = max(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(featureIdx,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)));
                    trial_rf_inst(featureIdx,trial,:) = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(featureIdx,trial,:))','linear'),bin_params.oldBin,[]));
                    autcorr_rf(featureIdx,trial,:) = xcov(trial_rf_inst(featureIdx,trial,:),'normalized');
                end
                for featureIdx = 7:8 % minAccX, minAccY
                    [trial_rf(featureIdx,trial), trial_rf_timing(featureIdx,trial)] = min(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(featureIdx-2,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)));
%                     trial_rf_inst(featureIdx,trial,:) = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(featureIdx-2,trial,:))','linear'),bin_params.oldBin,[]));
%                     autcorr_rf(featureIdx,trial,:) = xcov(trial_rf_inst(featureIdx,trial,:),'normalized');
                end
                
                % maxVel (mean of X and Y dimension)
                [trial_rf(9,trial), trial_rf_timing(9,trial)] = max(mean([abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins))); ...
                                              abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(4,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)))]));
                                          
                trial_rf_inst(7,trial,:) = mean([abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,:))','linear'),bin_params.oldBin,[]))); ...
                                              abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(4,trial,:))','linear'),bin_params.oldBin,[])))]);
                
%                 trial_rf_inst(7,trial,:) = mean([nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,:))','linear'),bin_params.oldBin,[])); ...
%                                               nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(4,trial,:))','linear'),bin_params.oldBin,[]))]);
                

                autcorr_rf(7,trial,:) = xcov(trial_rf_inst(7,trial,:),'normalized');
                                          
                % maxAcc (mean of X and Y dimension)
                [trial_rf(10,trial), trial_rf_timing(10,trial)] = max(mean([abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(5,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins))); ...
                    abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(6,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)))]));
                
                trial_rf_inst(8,trial,:) = mean([abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(5,trial,:))','linear'),bin_params.oldBin,[]))); ...
                                              abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(6,trial,:))','linear'),bin_params.oldBin,[])))]);
                                          
                autcorr_rf(8,trial,:) = xcov(trial_rf_inst(8,trial,:),'normalized');
                
                % minAcc (mean of X and Y dimension)
                [trial_rf(11,trial), trial_rf_timing(11,trial)] = min(mean([abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(5,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins))); ...
                    abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(6,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)))]));
                
                % Inst position
                trial_rf_inst(9,trial,:) = mean([abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(1,trial,:))','linear'),bin_params.oldBin,[]))); ...
                                              abs(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(2,trial,:))','linear'),bin_params.oldBin,[])))]);
                
                autcorr_rf(9,trial,:) = xcov(trial_rf_inst(9,trial,:),'normalized');

                % movDur & distTrav
                [~,i] = max(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)));
                trialXvel = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,:))','linear'),rf_params.binSize,numBins));            
                movOnset = max(find(trialXvel(1:numBins-(rf_params.prevBins-i))<=0));
%                 figure; plot(trialXvel); hold on; scatter(movOnset,trialXvel(movOnset),30,[1 0 0],'filled'); title(['trial # ' num2str(trial)]);

                if ~isempty(movOnset)
                    if nanCount(squeeze(trajectories{animal}{day}(3,trial,movOnset*rf_params.binSize:numBins*rf_params.binSize)))<rf_params.nan_thresh_global

                        trial_rf(12,trial) = (numBins-movOnset)*rf_params.binSize;
                        
                        trialXpos = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(1,trial,:))','linear'),rf_params.binSize,numBins));
                        trialYpos = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(2,trial,:))','linear'),rf_params.binSize,numBins));
                        distTrav = [];
                        for nBin = movOnset:numBins-1
                            distTrav = [distTrav pdist([trialXpos(nBin) trialYpos(nBin); trialXpos(nBin+1) trialYpos(nBin+1)],'euclidean')];
                        end
                        trial_rf(13,trial) = sum(distTrav);
                        
                    else
                        trial_rf(12,trial) = NaN;
                        trial_rf(13,trial) = NaN;
                    end
                else
                    trial_rf(12,trial) = NaN;
                    trial_rf(13,trial) = NaN;
                end
                
            else 
                
                trial_success(trial) = nan;
                
                for featureIdx = 1:13
                    trial_rf(featureIdx,trial) = nan;
                end
                
            end
        end
        
        % PCAstab
        PCAfeatures = trial_rf(1:6,:)'; % using max trial values of pos, vel, and acceleration (x and y)
        data_mean = mean(PCAfeatures,'omitnan'); % z-score
        data_std = std(PCAfeatures, 'omitnan');
        PCAfeatures = (PCAfeatures-data_mean)./(data_std);
        [coeff, ~, ~, ~, varExp] = pca(PCAfeatures,'Rows','complete'); % ignore nan values
        PCAtemplate = coeff(:,1:rf_params.num_PCs);
        trial_PCA_stability = sum((PCAfeatures*PCAtemplate).*(varExp(1:rf_params.num_PCs)'/sum(varExp(1:rf_params.num_PCs))),2);
        trial_rf = [trial_rf; trial_PCA_stability'];
        
        featureIDs = {'maxPosX','maxPosY','maxVelX','maxVelY','maxAccX','maxAccY','minAccX','minAccY','maxVel','maxAcc','minAcc','movDur','distTrav','PCAstab'};
        for featureIDx = 1:numel(featureIDs)
            reach_features{animal}(day).(featureIDs{featureIDx}) = trial_rf(featureIDx,:);
            if featureIDx <= 11 % argmax/min over time computed only for a subset of 11 features
                reach_features{animal}(day).featTime.(featureIDs{featureIDx}) = trial_rf_timing(featureIDx,:);
            end
        end
        reach_features{animal}(day).consNaNs = trial_cons_nans;
        reach_features{animal}(day).success_rate = trial_success;
        
        % instantaneous reach features
        featureInstIDs = {'posX','posY','velX','velY','accX','accY','vel','acc','pos'};
        for featureIDx = 1:numel(featureInstIDs)
            reach_features{animal}(day).instant.(featureInstIDs{featureIDx}) = squeeze(trial_rf_inst(featureIDx,:,:));
            reach_features{animal}(day).autcorr.(featureInstIDs{featureIDx}) = squeeze(autcorr_rf(featureIDx,:,:));
        end
    end
end

end
