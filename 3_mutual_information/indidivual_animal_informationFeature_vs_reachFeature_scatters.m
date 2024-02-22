load('G:\CurrBiol_revision\data\processed_data\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat');
load('G:\CurrBiol_revision\data\processed_data\realBinTimes.mat')

pellet_diameter = {8, 7, 8, 18, 15.5, 13, 15, 20}; % size of pellet in pixels (2.8121 pellet/cm)
late_minus_early_feature = nan(length(params.animals),6);
feature_count = 1;

for featureIdx = [9 13 10 12 1 2]
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
        end
        late_minus_early_feature(:,feature_count)=feature_earlyLate(:,2)-feature_earlyLate(:,1);
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
        end
        late_minus_early_feature(:,feature_count)=feature_earlyLate(:,2)-feature_earlyLate(:,1);
    else
        feature_earlyLate = [];
        for animalIdx = 1:length(params.animals)
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
        end
        late_minus_early_feature(:,feature_count)=feature_earlyLate(:,2)-feature_earlyLate(:,1);            
    end 
    feature_count = feature_count+1;
end

n_lfp = 1;
m1_dls_late_minus_early_info_mag = nan(2,numel(params.animals),6);
feature_count = 1;
for fidx = [9 13 10 12 1 2]
    for animal = 1:numel(params.animals)
        for aidx = 1:length(params.areas)
            early_late_info = [];
            for early_late = 1:2
                tmp_info = [];
                for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                    infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                    infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                    infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                    [v,i] = max(infQuant);
                    if v>0 && i>=12 && i<=86
                        tmp_info = [tmp_info; infQuant];
                    end
                end
                early_late_info = [early_late_info nanmean(max(tmp_info'))];
            end
            m1_dls_late_minus_early_info_mag(aidx,animal,feature_count) = early_late_info(2)-early_late_info(1);
        end
    end
    feature_count = feature_count + 1;
end

n_lfp = 1;
m1_dls_late_minus_early_info_timing_onset = nan(2,numel(params.animals),6);
m1_dls_late_minus_early_info_timing_peak = nan(2,numel(params.animals),6);
feature_count = 1;
for fidx = [9 13 10 12 1 2]
    for animal = 1:numel(params.animals)
        for aidx = 1:length(params.areas)
            early_late_max = [];
            early_late_onset = [];
            for early_late = 1:2
                tmp_max = [];
                tmp_onset = [];
                for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                    infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                    infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                    infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                    [v,i] = max(infQuant);
                    if v>0 && i>=12 && i<=86
                        baseSD = nanstd(infQuant(1:36));
                        baseMean = nanmean(infQuant(1:36));
                        [v, i] = max(infQuant(12:86)); % -1s to +500ms
                        if v>baseMean+(baseSD*4) & i~=length(12:86)
                            tmp_onset = [tmp_onset min(find(infQuant>v*0.66))];
                            tmp_max = [tmp_max 11+i];
                        end
                    end
                end
                early_late_max = [early_late_max mean(tmp_max)];
                early_late_onset = [early_late_onset mean(tmp_onset)];
            end
            m1_dls_late_minus_early_info_timing_peak(aidx,animal,feature_count) = early_late_max(2)-early_late_max(1);
            m1_dls_late_minus_early_info_timing_onset(aidx,animal,feature_count) = early_late_onset(2)-early_late_onset(1);
        end
    end
    feature_count = feature_count + 1;
end    

% late_minus_early_feature = nan(length(params.animals),6);
% m1_dls_late_minus_early_info_mag = nan(2,numel(params.animals),6);
% m1_dls_late_minus_early_info_timing_onset = nan(2,numel(params.animals),6);
% m1_dls_late_minus_early_info_timing_peak = nan(2,numel(params.animals),6);

for fidx = 1:6

    fig = figure('units','pixels','position',[200 200 800 800]);
        tiledlayout(3,1)
            nexttile; hold on;
                set(fig,'color','w'); 
                set(fig,'InvertHardcopy','off')                
                ax = gca;
                ax.YColor =[0 0 0];
                ax.XColor =[0 0 0];
                ax.Color =[1 1 1];
                ax.YAxis.FontSize = 16;
                ax.XAxis.FontSize = 16;
                ax.TickDir = 'out';
                ax.TickLength = [0.015 0.04];
                ax.LineWidth = 1.5;            
                scatter(late_minus_early_feature(:,fidx),m1_dls_late_minus_early_info_mag(1,:,fidx),50,[0 0 0],'filled')
                [R1,P1] = corrcoef(late_minus_early_feature(:,fidx),m1_dls_late_minus_early_info_mag(1,:,fidx))
                [p1,S1] = polyfit(late_minus_early_feature(:,fidx),m1_dls_late_minus_early_info_mag(1,:,fidx),1)
                [y_fit,delta] = polyval(p1,late_minus_early_feature(:,fidx),S1);
                plot(late_minus_early_feature(:,fidx),y_fit,'k','LineWidth',2)
                scatter(late_minus_early_feature(:,fidx),m1_dls_late_minus_early_info_mag(2,:,fidx),50,[1 0 0],'filled')
                [R2,P2] = corrcoef(late_minus_early_feature(:,fidx),m1_dls_late_minus_early_info_mag(2,:,fidx))
                [p2,S2] = polyfit(late_minus_early_feature(:,fidx),m1_dls_late_minus_early_info_mag(2,:,fidx),1)
                [y_fit,delta] = polyval(p2,late_minus_early_feature(:,fidx),S2);
                plot(late_minus_early_feature(:,fidx),y_fit,'r','LineWidth',2)
                title(['p = ' num2str(round(P1(1,2),2)) ' | p = ' num2str(round(P2(1,2),2))],'FontSize',20)
                xlabel('Feature Value')
                ylabel('Information magnitude')
            nexttile; hold on;
                set(fig,'color','w'); 
                set(fig,'InvertHardcopy','off')                
                ax = gca;
                ax.YColor =[0 0 0];
                ax.XColor =[0 0 0];
                ax.Color =[1 1 1];
                ax.YAxis.FontSize = 16;
                ax.XAxis.FontSize = 16;
                ax.TickDir = 'out';
                ax.TickLength = [0.015 0.04];
                ax.LineWidth = 1.5;  
                tmp_nans = (~isnan(m1_dls_late_minus_early_info_timing_onset(2,:,fidx)) & ~isnan(m1_dls_late_minus_early_info_timing_onset(1,:,fidx)));
                scatter(late_minus_early_feature(tmp_nans,fidx),m1_dls_late_minus_early_info_timing_onset(2,tmp_nans,fidx)-m1_dls_late_minus_early_info_timing_onset(1,tmp_nans,fidx),50,[0 0 0],'filled')
                [R1,P1] = corrcoef(late_minus_early_feature(tmp_nans,fidx),m1_dls_late_minus_early_info_timing_onset(2,tmp_nans,fidx)-m1_dls_late_minus_early_info_timing_onset(1,tmp_nans,fidx))
                [p1,S1] = polyfit(late_minus_early_feature(tmp_nans,fidx),m1_dls_late_minus_early_info_timing_onset(2,tmp_nans,fidx)-m1_dls_late_minus_early_info_timing_onset(1,tmp_nans,fidx),1)
                [y_fit,delta] = polyval(p1,late_minus_early_feature(tmp_nans,fidx),S1);
                plot(late_minus_early_feature(tmp_nans,fidx),y_fit,'k','LineWidth',2)
                tmp_nans = (~isnan(m1_dls_late_minus_early_info_timing_peak(1,:,fidx)) & ~isnan(m1_dls_late_minus_early_info_timing_peak(1,:,fidx)));                
                scatter(late_minus_early_feature(tmp_nans,fidx),m1_dls_late_minus_early_info_timing_peak(2,tmp_nans,fidx)-m1_dls_late_minus_early_info_timing_peak(1,tmp_nans,fidx),50,[1 0 0],'filled')
                [R2,P2] = corrcoef(late_minus_early_feature(tmp_nans,fidx),m1_dls_late_minus_early_info_timing_peak(2,tmp_nans,fidx)-m1_dls_late_minus_early_info_timing_peak(1,tmp_nans,fidx))
                [p2,S2] = polyfit(late_minus_early_feature(tmp_nans,fidx),m1_dls_late_minus_early_info_timing_peak(2,tmp_nans,fidx)-m1_dls_late_minus_early_info_timing_peak(1,tmp_nans,fidx),1)
                [y_fit,delta] = polyval(p2,late_minus_early_feature(tmp_nans,fidx),S2);
                plot(late_minus_early_feature(tmp_nans,fidx),y_fit,'r','LineWidth',2)
                title(['p = ' num2str(round(P1(1,2),2)) ' | p = ' num2str(round(P2(1,2),2))],'FontSize',16)
                xlabel('Feature Value')
                ylabel('Information Timing')

                if fidx==1
                    sgtitle('max velocity','FontSize',30)
                elseif fidx==2
                    sgtitle('traj length','FontSize',30)                   
                elseif fidx==3
                    sgtitle('max acceleration','FontSize',30)                    
                elseif fidx==4
                    sgtitle('movement duration','FontSize',30)                    
                elseif fidx==5
                    sgtitle('max x position','FontSize',30)                    
                elseif fidx==6
                    sgtitle('max y position','FontSize',30)                    
                end
                
                print(['C:\Users\StefanLemke\Desktop\by_animal_reachFeature_vs_infoFeature_Feature' num2str(fidx) '.eps'],'-painters','-depsc');
     
end

    