%% COMPUTE NEURAL ACTIVITY 
% Pre-computes and generates one .mat file per animal with trial-aligned 
% neural activity (both from spikes and LFP) on each day of training for 
% all trials newTimeBin is a fundamental parameter here because the bigger
% it is, the lighter are the output files.

    clear all; clc; close all;

%%% INITIALIZE PARAMETERS

    % data paramaters
    animals = {'T102','T107','T110','T200','T201','SLDD1','SLDD2','T398'};
    days = {[1:8],[2:6],[2:7],[2:10],[4:17],[3:14],[4:13],[3:9]};
    days{3}(3) = []; days{8}(4) = []; % remove missing day  
    path = 'D:\MEGA\data';

    % general paramaters     
    areas = {'M1','DLS'};
    lfp_feature = {'reach_related_lfp','reach_related_lfp_ref'}; % either 'reach_related_lfp' for non-referenced LFP or 'reach_related_lfp_ref' for referenced LFP    
    neural_signals = {'spikes','lfp'};
    make_plot = 0;
    doSpikes = 1;
    doLFP = 1; 

    % temporal parameters    
    tMin = 3501;
    tMax = 6000;
    newTimeBin = 1; %100ms
    bin2sec = 1000/newTimeBin; % original sampling rate is 1000Hz
    tPoints.spikes = tMax-tMin+1;

%%% COMPUTE

    for animal = 1:numel(animals)

        disp(['Animal ', num2str(animal)])
        animal_days = days{animal}; 
        dayCount = 1;
        spikeTrains = cell(1,numel(animal_days));
        LFP = cell(1,numel(animal_days));

        for day = animal_days

            nDays(animal) = numel(animal_days);

            if doSpikes
                % M1 spikeTrains
                load([path '\' animals{animal} '\neural_data_reach\' animals{animal} '_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking')
                m1_reach_spiking = rmfield(m1_reach_spiking,'baseline_spiking');
                nTrials(dayCount) = size(m1_reach_spiking(1).reach_related_spiking,1);
                nUnits.M1(dayCount) = length(m1_reach_spiking);

                spikeTrainsOldBin.M1 = zeros(nUnits.M1(dayCount),tPoints.spikes,nTrials(dayCount));
                for uidx = 1:nUnits.M1(dayCount)
                    spikeTrainsOldBin.M1(uidx,:,:) = m1_reach_spiking(uidx).reach_related_spiking(:,tMin:tMax)';
                end

                spikeTrains{dayCount}.M1 = temporal_rebinning(spikeTrainsOldBin.M1, newTimeBin);
                spikeTrains{dayCount}.M1(spikeTrains{dayCount}.M1>1) = 1;
                tPointsNew.spikes = size(spikeTrains{dayCount}.M1,2);

                clear m1_reach_spiking

                % DLS spikeTrains
                load([path '\' animals{animal} '\neural_data_reach\' animals{animal} '_day_' num2str(day) '_reach_spiking_LFP.mat'],'dls_reach_spiking')
                dls_reach_spiking = rmfield(dls_reach_spiking,'baseline_spiking');

                nUnits.DLS(dayCount) = length(dls_reach_spiking);

                spikeTrainsOldBin.DLS = zeros(nUnits.DLS(dayCount),tPoints.spikes,nTrials(dayCount));
                for uidx = 1:nUnits.DLS(dayCount)
                    spikeTrainsOldBin.DLS(uidx,:,:) = dls_reach_spiking(uidx).reach_related_spiking(:,tMin:tMax)';
                end

                spikeTrains{dayCount}.DLS = temporal_rebinning(spikeTrainsOldBin.DLS, newTimeBin);
                spikeTrains{dayCount}.DLS(spikeTrains{dayCount}.DLS>1) = 1;
                tPointsNew.spikes = size(spikeTrains{dayCount}.DLS,2);

                clear dls_reach_spiking spikeTrainsOldBin
            end

            if doLFP
                % M1 LFP
                load([path '\' animals{animal} '\neural_data_reach\' animals{animal} '_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_lfp')
                m1_reach_lfp = rmfield(m1_reach_lfp,{'baseline_lfp','baseline_lfp_ref'});
                nLFPchan.M1(dayCount) = numel(m1_reach_lfp);
                for LFP_type = 1:2

                    selFeature = lfp_feature{LFP_type};

                    % GET LFP SIGNALS
                    tPoints.lfp = size(m1_reach_lfp(1).(selFeature),2);

                    LFPOldBin.M1 = zeros(nLFPchan.M1(dayCount),tPoints.lfp,nTrials(dayCount));

                    for uidx = 1:nLFPchan.M1(dayCount)
                        LFPOldBin.M1(uidx,:,:) = m1_reach_lfp(uidx).(selFeature)(:,:)';
                    end

                    % Interpolate LFP to 20000 points to, then, resample it
                    % aligning with trajectories
                    orig_time = 1:tPoints.lfp;
                    interp_time = linspace(1,tPoints.lfp,20000);
                    LFP{dayCount}.(selFeature).M1 = interp1(orig_time, permute(LFPOldBin.M1,[2 1 3]), interp_time);

                    LFP{dayCount}.(selFeature).M1 = permute(LFP{dayCount}.(selFeature).M1,[2 1 3]);
                    clear tmpLFP

                    % temporal rebinning (here I average the raw values over the
                    % required temporal bin newTimeBin)
                    LFP{dayCount}.(selFeature).M1 = temporal_rebinning(LFP{dayCount}.(selFeature).M1(:,tMin*2-1:tMax*2,:), newTimeBin*2,'mean');
                    tPointsNew.lfp = size(LFP{dayCount}.(selFeature).M1,2);
                end
                clear m1_reach_lfp

                % DLS LFP
                load([path '\' animals{animal} '\neural_data_reach\' animals{animal} '_day_' num2str(day) '_reach_spiking_LFP.mat'],'dls_reach_lfp')          
                dls_reach_lfp = rmfield(dls_reach_lfp,{'baseline_lfp','baseline_lfp_ref'});

                nLFPchan.DLS(dayCount) = numel(dls_reach_lfp);
                for LFP_type = 1:2
                    selFeature = lfp_feature{LFP_type};
                    % GET LFP SIGNALS
                    tPoints.lfp = size(dls_reach_lfp(1).(selFeature),2);

                    LFPOldBin.DLS = zeros(nLFPchan.DLS(dayCount),tPoints.lfp,nTrials(dayCount));

                    for uidx = 1:nLFPchan.DLS(dayCount)
                        LFPOldBin.DLS(uidx,:,:) = dls_reach_lfp(uidx).(selFeature)(:,:)';
                    end
                    % Interpolate LFP to 20000 points to, then, resample it
                    % aligning with trajectories
                    orig_time = 1:tPoints.lfp;
                    interp_time = linspace(1,tPoints.lfp,20000);
                    LFP{dayCount}.(selFeature).DLS = interp1(orig_time, permute(LFPOldBin.DLS,[2 1 3]), interp_time);

                    LFP{dayCount}.(selFeature).DLS = permute(LFP{dayCount}.(selFeature).DLS,[2 1 3]);
                    clear tmpLFP

                    % temporal rebinning (here I average the raw values over the
                    % required temporal bin newTimeBin)
                    LFP{dayCount}.(selFeature).DLS = temporal_rebinning(LFP{dayCount}.(selFeature).DLS(:,tMin*2-1:tMax*2,:), newTimeBin*2,'mean');
                    tPointsNew.lfp = size(LFP{dayCount}.(selFeature).M1,2);
                end
                clear dls_reach_lfp
            end
            dayCount = dayCount + 1;
        end
        fname = ['neuralActivity_',animals{animal},'_',num2str(newTimeBin),'ms_',date];
        save(['C:\Users\slemke\Documents\InfoFlowProject\data\preprocessed_data\' fname], 'spikeTrains', 'LFP', 'newTimeBin','tMin','tMax','-v7.3');

    end