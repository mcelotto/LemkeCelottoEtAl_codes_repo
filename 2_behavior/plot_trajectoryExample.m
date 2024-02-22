function [] = plot_trajectoryExample(params,rf_params,paths,save_params)

% T102 | 0.0711 cm/pixel | 30 frames/second
% T107 | 0.0711 cm/pixel | 30 frames/second
% T110 | 0.0711 cm/pixel | 30 frames/second
% T200 | 0.0222 cm/pixel | 60 frames/second
% T201 | 0.0254 cm/pixel | 60 frames/second
% SLDD1 | 0.0237 cm/pixel | 75 frames/second
% SLDD2 | 0.0237 cm/pixel | 75 frames/second
% T398 | 0.0127 cm/pixel | 100 frames/second

load(paths.trajectories_forExample); % load trajectory data
animalColors = distinguishable_colors(8);

for animal = 8
    
    early_posX = [];
    early_posY = [];
    early_velX = [];
    early_velY = [];
    for day = 1:params.num_earlylate_days{animal}
        tPoints = 1+size(trajectories{animal}{day},3)+(5000-tMax)-rf_params.binSize*rf_params.prevBins:size(trajectories{animal}{day},3)+(5000-tMax); % time points used for NaN trial disqualification and max/min features
        for trial = 1:size(trajectories{animal}{day},2)
            [~,trial_cons_nans(trial)] = nanCount(squeeze(trajectories{animal}{day}(1,trial,tPoints)));
            if trial_cons_nans(trial)<rf_params.nan_thresh
                
                early_velX = [early_velX; smooth(fillmissing(squeeze(trajectories{animal}{day}(3,trial,:))','movmean',3),10)'];
                early_velY = [early_velY; smooth(fillmissing(squeeze(trajectories{animal}{day}(4,trial,:))','movmean',3),10)'];
                
                [~,i] = max(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)));
                trialXvel = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,1:1000))','linear'),rf_params.binSize,20));            
                movOnset = max(find(trialXvel(1:20-(rf_params.prevBins-i))<=0));
                
                tmpX = smooth(fillmissing(squeeze(trajectories{animal}{day}(1,trial,movOnset*50:1100))','movmean',3),10)';
                tmpY = smooth(fillmissing(squeeze(trajectories{animal}{day}(2,trial,movOnset*50:1100))','movmean',3),10)';
                
                early_posX = [early_posX; nan(1,1100-length(tmpX)) tmpX];
                early_posY = [early_posY; nan(1,1100-length(tmpX)) tmpY];
                
            end
        end
    end

    late_posX = [];
    late_posY = [];
    late_velX = [];
    late_velY = [];
    late_accX = [];
    late_accY = [];
    for day = numel(trajectories{animal})-params.num_earlylate_days{animal}+1:numel(trajectories{animal})
        tPoints = 1+size(trajectories{animal}{day},3)+(5000-tMax)-rf_params.binSize*rf_params.prevBins:size(trajectories{animal}{day},3)+(5000-tMax); % time points used for NaN trial disqualification and max/min features
        for trial = 1:size(trajectories{animal}{day},2)
            [~,trial_cons_nans(trial)] = nanCount(squeeze(trajectories{animal}{day}(1,trial,tPoints)));
            if trial_cons_nans(trial)<rf_params.nan_thresh

                late_velX = [late_velX; smooth(fillmissing(squeeze(trajectories{animal}{day}(3,trial,:))','movmean',3),10)'];
                late_velY = [late_velY; smooth(fillmissing(squeeze(trajectories{animal}{day}(4,trial,:))','movmean',3),10)'];
                
                [~,i] = max(nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,tPoints))','linear'),rf_params.binSize,rf_params.prevBins)));
                trialXvel = nanmean(reshape(fillmissing(squeeze(trajectories{animal}{day}(3,trial,1:1000))','linear'),rf_params.binSize,20));            
                movOnset = max(find(trialXvel(1:20-(rf_params.prevBins-i))<=0));
                
                tmpX = smooth(fillmissing(squeeze(trajectories{animal}{day}(1,trial,movOnset*50:1100))','movmean',3),10)';
                tmpY = smooth(fillmissing(squeeze(trajectories{animal}{day}(2,trial,movOnset*50:1100))','movmean',3),10)';
                
                late_posX = [late_posX; nan(1,1100-length(tmpX)) tmpX];
                late_posY = [late_posY; nan(1,1100-length(tmpX)) tmpY];
                
            end
        end
    end

    min_trial = min([size(early_posX,1) size(late_posX,1)]);
    early_posX = early_posX(1:min_trial,:);
    early_posY = early_posY(1:min_trial,:);
    late_posX = late_posX(end-min_trial+1:end,:);
    late_posY = late_posY(end-min_trial+1:end,:);
    early_velX = early_velX(1:min_trial,:);
    early_velY = early_velY(1:min_trial,:);
    late_velX = late_velX(end-min_trial+1:end,:);
    late_velY = late_velY(end-min_trial+1:end,:);

end

figure;
    subplot(1,2,1); hold on;
        for trial = 1:50%100
            plot(early_posX(trial,:),early_posY(trial,:),'color',[.5 .5 .5])
        end
        plot(nanmean(early_posX(:,701:1100)),nanmean(early_posY(:,701:1100)),'color','k','LineWidth',2)
        scatter(nanmean(early_posX(:,1000)),nanmean(early_posY(:,1000)),30,[1 0 0],'filled')
        xlim([-400 50])
        ylim([-160 75])
        title('early')
	subplot(1,2,2); hold on;
        for trial = 1:50%size(late_posX,1)-99:size(late_posX,1)
            plot(late_posX(trial,:),late_posY(trial,:),'color',[.5 .5 .5])
        end
        plot(nanmean(late_posX(:,701:1100)),nanmean(late_posY(:,701:1100)),'color','k','LineWidth',2)
        scatter(nanmean(late_posX(:,1000)),nanmean(late_posY(:,1000)),30,[1 0 0],'filled')
        xlim([-400 50])
        ylim([-160 75])
        title('late')
        
        if save_params.save == 1
            print([save_params.save_path '\successExample_traj.eps'],'-painters','-depsc');
        end


early_velX_cm = early_velX.*0.0127*1000;
early_velY_cm = early_velY.*0.0127*1000;
late_velX_cm = late_velX.*0.0127*1000;
late_velY_cm = late_velY.*0.0127*1000; 
    
figure;

    subplot(2,1,1); hold on;

        plot(nanmean(early_velX_cm(1:200,601:1100))+(nanstd(early_velX_cm(1:200,601:1100))/sqrt(size(early_velX_cm(1:200,601:1100),1))),'color','k','LineWidth',2)
        plot(nanmean(early_velX_cm(1:200,601:1100))-(nanstd(early_velX_cm(1:200,601:1100))/sqrt(size(early_velX_cm(1:200,601:1100),1))),'color','k','LineWidth',2)
        plot(nanmean(late_velX_cm(1:200,601:1100))+(nanstd(late_velX_cm(1:200,601:1100))/sqrt(size(late_velX_cm(1:200,601:1100),1))),'color','r','LineWidth',2)
        plot(nanmean(late_velX_cm(1:200,601:1100))-(nanstd(late_velX_cm(1:200,601:1100))/sqrt(size(late_velX_cm(1:200,601:1100),1))),'color','r','LineWidth',2)

        xlim([1 500])
        ylim([-10 20])
        title('X dimension')

    subplot(2,1,2); hold on;

        plot(nanmean(early_velY_cm(1:200,601:1100))+(nanstd(early_velY_cm(1:200,601:1100))/sqrt(size(early_velY_cm(1:200,601:1100),1))),'color','k','LineWidth',2)
        plot(nanmean(early_velY_cm(1:200,601:1100))-(nanstd(early_velY_cm(1:200,601:1100))/sqrt(size(early_velY_cm(1:200,601:1100),1))),'color','k','LineWidth',2)
        plot(nanmean(late_velY_cm(1:200,601:1100))+(nanstd(late_velY_cm(1:200,601:1100))/sqrt(size(late_velY_cm(1:200,601:1100),1))),'color','r','LineWidth',2)
        plot(nanmean(late_velY_cm(1:200,601:1100))-(nanstd(late_velY_cm(1:200,601:1100))/sqrt(size(late_velY_cm(1:200,601:1100),1))),'color','r','LineWidth',2)

        xlim([1 500])
        ylim([-10 20])
        title('Y dimension')


%         trialCount = 1;
%         trialPlot = 1;
%         while trialPlot<101
%             if max(early_velX_cm(trialCount,601:1100))<40 && min(early_velX_cm(trialCount,601:1100))>-20
%                 plot(early_velX_cm(trialCount,601:1100),'color',[.5 .5 .5])
%                 trialPlot = trialPlot + 1;
%                 trialCount = trialCount + 1;
%             else
%                 trialCount = trialCount + 1;
%             end
%         end

        
%     subplot(2,2,2); hold on;
%         trialCount = 1;
%         trialPlot = 1;
%         while trialPlot<101
%             if max(late_velX_cm(trialCount,601:1100))<40 && min(late_velX_cm(trialCount,601:1100))>-20
%                 plot(late_velX_cm(trialCount,601:1100),'color',[.5 .5 .5])
%                 trialPlot = trialPlot + 1;
%                 trialCount = trialCount + 1;
%             else
%                 trialCount = trialCount + 1;
%             end
%         end
%         plot(nanmean(late_velX_cm(:,601:1100)),'color','k','LineWidth',2)
%         xlim([1 500])
%         ylim([-20 40])
%         title('late')
        
%     subplot(2,2,3); hold on;
%         trialCount = 1;
%         trialPlot = 1;
%         while trialPlot<101
%             if max(early_velY_cm(trialCount,601:1100))<40 && min(early_velY_cm(trialCount,601:1100))>-20
%                 plot(early_velY_cm(trialCount,601:1100),'color',[.5 .5 .5])
%                 trialPlot = trialPlot + 1;
%                 trialCount = trialCount + 1;
%             else
%                 trialCount = trialCount + 1;
%             end
%         end
%         plot(nanmean(early_velY_cm(:,601:1100)),'color','k','LineWidth',2)
%         xlim([1 500])
%         ylim([-20 40])
%         title('early')
%         
%     subplot(2,2,4); hold on;
%         trialCount = 1;
%         trialPlot = 1;
%         while trialPlot<101
%             if max(late_velY_cm(trialCount,601:1100))<40 && min(late_velY_cm(trialCount,601:1100))>-20
%                 plot(late_velY_cm(trialCount,601:1100),'color',[.5 .5 .5])
%                 trialPlot = trialPlot + 1;
%                 trialCount = trialCount + 1;
%             else
%                 trialCount = trialCount + 1;
%             end
%         end
%         plot(nanmean(late_velY_cm(:,601:1100)),'color','k','LineWidth',2)
%         xlim([1 500])
%         ylim([-20 40])
%         title('late')
        
        if save_params.save == 1
            print([save_params.save_path '\trajVel_vel.eps'],'-painters','-depsc');
        end
        
end
