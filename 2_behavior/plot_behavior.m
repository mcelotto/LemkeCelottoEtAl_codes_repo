function [] = plot_behavior(reach_features,params,save_params)

animalColors = distinguishable_colors(8);
pellet_diameter = {8, 7, 8, 18, 15.5, 13, 15, 20}; % size of pellet in pixels (2.8121 pellet/cm)

for featureIdx = 1:length(params.reachFeatures)

    figure('units','normalized','outerposition',[0 0 1 1]);
    
    % plot reach feature over days

        for animalIdx = 1:length(params.animals)

            if animalIdx<5
                subplot(2,5,animalIdx); hold on;
            else
                subplot(2,5,animalIdx+1); hold on;
            end

            cm_per_px = 1/(2.8121*pellet_diameter{animalIdx});

            if ismember(featureIdx,[3:11])

                dayVals = [];    
                for day = 1:numel(reach_features{animalIdx})
                    scatter(day,nanmean(reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px*1000),30,animalColors(animalIdx,:),'filled'); % convert px/bin to cm/sec
                    dayVals = [dayVals nanmean(reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px*1000)];
                end
                plot(dayVals,'color',animalColors(animalIdx,:))

            elseif ismember(featureIdx,[1 2 13])

                dayVals = [];    
                for day = 1:numel(reach_features{animalIdx})
                    scatter(day,nanmean(reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px),30,animalColors(animalIdx,:),'filled'); % convert px/bin to cm/sec
                    dayVals = [dayVals nanmean(reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px)];
                end
                plot(dayVals,'color',animalColors(animalIdx,:))

            else

                dayVals = [];    
                for day = 1:numel(reach_features{animalIdx})
                    scatter(day,nanmean(reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})),30,animalColors(animalIdx,:),'filled'); % convert px/bin to cm/sec
                    dayVals = [dayVals nanmean(reach_features{animalIdx}(day).(params.reachFeatures{featureIdx}))];
                end
                plot(dayVals,'color',animalColors(animalIdx,:))                

            end

        end

    % plot early late comparison 

    subplot(2,5,[5 10]); hold on;

        if ismember(featureIdx,[3:11])

            feature_earlyLate = [];
    
            for animalIdx = 1:length(params.animals)
            
                cm_per_px = 1/(2.8121*pellet_diameter{animalIdx});
                
                % early
                tmp_early = [];
                for day = 1:params.num_earlylate_days{animalIdx}
                    tmp_early = [tmp_early reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px*1000];
                end
                
                % late
                tmp_late = [];
                for day = numel(reach_features{animalIdx})-params.num_earlylate_days{animalIdx}+1:numel(reach_features{animalIdx})
                    tmp_late = [tmp_late reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px*1000];
                end    
    
                tmp_match = min([length(tmp_early) length(tmp_late)]);
                tmp_early = tmp_early(1:tmp_match);
                tmp_late = tmp_late(end-tmp_match+1:end);
                feature_earlyLate = [feature_earlyLate; nanmean(tmp_early) nanmean(tmp_late)];
    
                plot([1 2],[nanmean(tmp_early) nanmean(tmp_late)],'LineWidth',1','color',animalColors(animalIdx,:));
                scatter(1,nanmean(tmp_early),30,animalColors(animalIdx,:),'filled')
                scatter(2,nanmean(tmp_late),30,animalColors(animalIdx,:),'filled')
    
            end

        elseif ismember(featureIdx,[1 2 13])
            
            feature_earlyLate = [];
    
            for animalIdx = 1:length(params.animals)

                cm_per_px = 1/(2.8121*pellet_diameter{animalIdx});
                
                % early
                tmp_early = [];
                for day = 1:params.num_earlylate_days{animalIdx}
                    tmp_early = [tmp_early reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px];
                end
                
                % late
                tmp_late = [];
                for day = numel(reach_features{animalIdx})-params.num_earlylate_days{animalIdx}+1:numel(reach_features{animalIdx})
                    tmp_late = [tmp_late reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})*cm_per_px];
                end    
    
                tmp_match = min([length(tmp_early) length(tmp_late)]);
                tmp_early = tmp_early(1:tmp_match);
                tmp_late = tmp_late(end-tmp_match+1:end);
                feature_earlyLate = [feature_earlyLate; nanmean(tmp_early) nanmean(tmp_late)];
    
                plot([1 2],[nanmean(tmp_early) nanmean(tmp_late)],'LineWidth',1','color',animalColors(animalIdx,:));
                scatter(1,nanmean(tmp_early),30,animalColors(animalIdx,:),'filled')
                scatter(2,nanmean(tmp_late),30,animalColors(animalIdx,:),'filled')
    
            end

        else
            
            feature_earlyLate = [];
    
            for animalIdx = 1:length(params.animals)

%                 cm_per_px = 1/(2.8121*pellet_diameter{animalIdx});
                
                % early
                tmp_early = [];
                for day = 1:params.num_earlylate_days{animalIdx}
                    tmp_early = [tmp_early reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})];
                end
                
                % late
                tmp_late = [];
                for day = numel(reach_features{animalIdx})-params.num_earlylate_days{animalIdx}+1:numel(reach_features{animalIdx})
                    tmp_late = [tmp_late reach_features{animalIdx}(day).(params.reachFeatures{featureIdx})];
                end    
    
                tmp_match = min([length(tmp_early) length(tmp_late)]);
                tmp_early = tmp_early(1:tmp_match);
                tmp_late = tmp_late(end-tmp_match+1:end);
                feature_earlyLate = [feature_earlyLate; nanmean(tmp_early) nanmean(tmp_late)];
    
                plot([1 2],[nanmean(tmp_early) nanmean(tmp_late)],'LineWidth',1','color',animalColors(animalIdx,:));
                scatter(1,nanmean(tmp_early),30,animalColors(animalIdx,:),'filled')
                scatter(2,nanmean(tmp_late),30,animalColors(animalIdx,:),'filled')
    
            end

        end 

        errorbar(.75,mean(feature_earlyLate(:,1)),std(feature_earlyLate(:,1))/sqrt(numel(feature_earlyLate(:,1))),'color','k','LineWidth',2);
        errorbar(2.25,mean(feature_earlyLate(:,2)),std(feature_earlyLate(:,2))/sqrt(numel(feature_earlyLate(:,2))),'color','k','LineWidth',2);
        xlim([.5 2.5]);
        [h p1 ci stats] = ttest(feature_earlyLate(:,1),feature_earlyLate(:,2));
        [p2 h] = signrank(feature_earlyLate(:,1),feature_earlyLate(:,2));
        title([' paired t-test: ' num2str(p1) ' | paired rank-sum: ' num2str(p2)])

    disp([params.reachFeatures{featureIdx} ...
        ' | Early: ' num2str(mean(feature_earlyLate(:,1))) ' +- ' num2str(std(feature_earlyLate(:,1))/sqrt(numel(feature_earlyLate(:,1)))) ...
        ' Late: ' num2str(mean(feature_earlyLate(:,2))) ' +- ' num2str(std(feature_earlyLate(:,2))/sqrt(numel(feature_earlyLate(:,2)))) ...
        ' | t(' num2str(stats.df) ')=' num2str(stats.tstat) ', P=' num2str(p1) ' P=' num2str(p2)]);

    sgtitle((params.reachFeatures{featureIdx}));

    if featureIdx == 15
        for subPlotIdx = [1:4 6:9]
            subplot(2,5,subPlotIdx)
            ylim([0 1])
        end
        subplot(2,5,[5 10])
        ylim([0 1])
    end

    if save_params.save == 1
        print([save_params.save_path '\' (params.reachFeatures{featureIdx}) '_plot.eps'],'-painters','-depsc');
    end

end

end 