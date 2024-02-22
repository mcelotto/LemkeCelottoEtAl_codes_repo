function [] = plot_PID_Spikes(PIDpath,params,reach_bins,plot_params,clusterParams)

early_info_all = [];
late_info_all = [];
early_info_all_noZero = [];
late_info_all_noZero = [];

tmpFiles = ls([PIDpath]);
tmpFiles = tmpFiles([5:9 3 4 10],:);

for animal = 1:size(tmpFiles,1)

    load([PIDpath '\' tmpFiles(animal,:)])
    params.reachFeatures = {'maxVel','distTrav'};
    if isempty(find(~cellfun(@isempty,PIDout)))
        continue
    end
    PIDout = PIDout{find(~cellfun(@isempty,PIDout))};

    for fidx = 1:2
    
        % early
        for day = 1:params.num_earlylate_days{animal}
            if day<=length(PIDout.spikes)
                if isfield(PIDout.spikes{day},(params.reachFeatures{fidx}))
                    for unit = 1:size(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared,1)
                        if sum(sum(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:)))>0
                            tmp_info = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:));
                            tmp_infoSh = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).sharedSh(unit,:,:,:));
                            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                            early_info_all_noZero = cat(3,early_info_all_noZero,tmp_info);
                            tmp_info(~tmp_sig) = 0;
                            early_info_all = cat(3,early_info_all,tmp_info);
                        end
                    end
                end
            end
        end
    
        % late
        for day = length(params.days{animal})-params.num_earlylate_days{animal}+1:length(params.days{animal})
            if day<=length(PIDout.spikes)
                if isfield(PIDout.spikes{day},(params.reachFeatures{fidx}))
                    for unit = 1:size(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared,1)
                        if sum(sum(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:)))>0
                            tmp_info = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:));
                            tmp_infoSh = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).sharedSh(unit,:,:,:));
                            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                            late_info_all_noZero = cat(3,late_info_all_noZero,tmp_info);
                            tmp_info(~tmp_sig) = 0;
                            late_info_all = cat(3,late_info_all,tmp_info);
                        end
                    end
                end
            end
        end

    end
end

tmp_delay_early = [];
for pair = 1:size(early_info_all,3)
    [v, i] = max(median(squeeze(early_info_all(reach_bins,:,pair)),1));
    [v2, i2] = max(median(squeeze(early_info_all_noZero(reach_bins,:,pair)),1));
    if v>0 & i2>1 & i2<51
        tmp_delay_early = [tmp_delay_early i];
    else
        tmp_delay_early = [tmp_delay_early nan];
    end
end

tmp_delay_late = [];
for pair = 1:size(late_info_all,3)
    [v, i] = max(median(squeeze(late_info_all(reach_bins,:,pair)),1));
    [v2, i2] = max(median(squeeze(late_info_all_noZero(reach_bins,:,pair)),1));
    if v>0 & i2>1 & i2<51
        tmp_delay_late = [tmp_delay_late i];
    else
        tmp_delay_late = [tmp_delay_late nan];
    end
end

figure;
subplot(2,2,[1 2]); hold on;
histogram(tmp_delay_early(~isnan(tmp_delay_early)),[1:5:51],'normalization','probability','DisplayStyle','Stairs')
histogram(tmp_delay_late(~isnan(tmp_delay_late)),[1:5:51],'normalization','probability','DisplayStyle','Stairs')
xline(26)
ylim([0 0.25])
subplot(2,2,3); hold on;
bar(1,sum(tmp_delay_early<26)/sum(~isnan(tmp_delay_early)))
bar(2,sum(tmp_delay_late<26)/sum(~isnan(tmp_delay_late)))
ylim([0 .8])
subplot(2,2,4); hold on;
bar(1,sum(tmp_delay_early>26)/sum(~isnan(tmp_delay_early)))
bar(2,sum(tmp_delay_late>26)/sum(~isnan(tmp_delay_late)))
ylim([0 .8])
[h p] = kstest2(tmp_delay_early(~isnan(tmp_delay_early)),tmp_delay_late(~isnan(tmp_delay_late)));
sgtitle([num2str(length(tmp_delay_early)) ' | ' num2str(sum(~isnan(tmp_delay_early))) ' | ' num2str(length(tmp_delay_late)) ' | ' num2str(sum(~isnan(tmp_delay_late))) ' | p: ' num2str(p)]);

if plot_params.save == 1
    print(['G:\CurrBiol_revision\figures\figure_4\PID_spikes_histogram_bar_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
end

% 
% 
% 
% 
% 
% 
% for fidx = [1 2]
% 
%     early_info_all = [];
%     late_info_all = [];
%     early_info_all_noZero = [];
%     late_info_all_noZero = [];
% 
%     tmpFiles = ls([PIDpath]);
%     if fidx == 1
%         tmpFiles = tmpFiles([5 6 8 3 4 9],:);
%     elseif fidx == 2
%         tmpFiles = tmpFiles([5 6 7 3 4 9],:);
%     end
%     early_info_all = [];
%     late_info_all = [];
%     early_info_all_noZero = [];
%     late_info_all_noZero = [];
%     for animal = 1:size(tmpFiles,1)
% 
%         load([PIDpath '\' tmpFiles(animal,:)])
%         params.reachFeatures = {'maxVel','distTrav'};
%         if isempty(find(~cellfun(@isempty,PIDout)))
%             continue
%         end
%         PIDout = PIDout{find(~cellfun(@isempty,PIDout))};
% 
%         % early
%         for day = 1:params.num_earlylate_days{animal}
%             if day<=length(PIDout.spikes)
%                 if isfield(PIDout.spikes{day},(params.reachFeatures{fidx}))
%                     for unit = 1:size(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared,1)
%                         if sum(sum(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:)))>0
%                             tmp_info = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:));
%                             tmp_infoSh = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).sharedSh(unit,:,:,:));
%                             tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
%                             tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%                             early_info_all_noZero = cat(3,early_info_all_noZero,tmp_info);
%                             tmp_info(~tmp_sig) = 0;
%                             early_info_all = cat(3,early_info_all,tmp_info);
%                         end
%                     end
%                 end
%             end
%         end
% 
%         % late
%         for day = length(params.days{animal})-params.num_earlylate_days{animal}+1:length(params.days{animal})
%             if day<=length(PIDout.spikes)
%                 if isfield(PIDout.spikes{day},(params.reachFeatures{fidx}))
%                     for unit = 1:size(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared,1)
%                         if sum(sum(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:)))>0
%                             tmp_info = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:));
%                             tmp_infoSh = squeeze(PIDout.spikes{day}.(params.reachFeatures{fidx}).sharedSh(unit,:,:,:));
%                             tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
%                             tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%                             late_info_all_noZero = cat(3,late_info_all_noZero,tmp_info);
%                             tmp_info(~tmp_sig) = 0;
%                             late_info_all = cat(3,late_info_all,tmp_info);
%                         end
%                     end
%                 end
%             end
%         end
%     end
% 
%     tmp_delay_early = [];
%     for pair = 1:size(early_info_all,3)
%         [v, i] = max(median(squeeze(early_info_all(reach_bins,:,pair)),1));
%         [v2, i2] = max(median(squeeze(early_info_all_noZero(reach_bins,:,pair)),1));
%         if v>0 & i>1& i<51
%             tmp_delay_early = [tmp_delay_early i2];
%         else
%             tmp_delay_early = [tmp_delay_early nan];
%         end
%     end
% 
%     tmp_delay_late = [];
%     for pair = 1:size(late_info_all,3)
%         [v, i] = max(median(squeeze(late_info_all(reach_bins,:,pair)),1));
%         [v2, i2] = max(median(squeeze(late_info_all_noZero(reach_bins,:,pair)),1));
%         if v>0 & i>1 & i<51
%             tmp_delay_late = [tmp_delay_late i2];
%         else
%             tmp_delay_late = [tmp_delay_late nan];
%         end
%     end
% 
%     figure;
%     subplot(2,2,[1 2]); hold on;
%     histogram(tmp_delay_early(~isnan(tmp_delay_early)),[1:4:51],'normalization','probability','DisplayStyle','Stairs')
%     histogram(tmp_delay_late(~isnan(tmp_delay_late)),[1:4:51],'normalization','probability','DisplayStyle','Stairs')
%     xline(26)
%     ylim([0 0.25])
%     subplot(2,2,3); hold on;
%     bar(1,sum(tmp_delay_early<26)/sum(~isnan(tmp_delay_early)))
%     bar(2,sum(tmp_delay_late<26)/sum(~isnan(tmp_delay_late)))
%     ylim([0 1])
%     subplot(2,2,4); hold on;
%     bar(1,sum(tmp_delay_early>26)/sum(~isnan(tmp_delay_early)))
%     bar(2,sum(tmp_delay_late>26)/sum(~isnan(tmp_delay_late)))
%     ylim([0 1])
%     [h p] = kstest2(tmp_delay_early(~isnan(tmp_delay_early)),tmp_delay_late(~isnan(tmp_delay_late)));
%     sgtitle([num2str(reach_bins) ' | ' num2str(length(tmp_delay_early)) ' | ' num2str(sum(~isnan(tmp_delay_early))) ' | ' num2str(length(tmp_delay_late)) ' | ' num2str(sum(~isnan(tmp_delay_late))) ' | p: ' num2str(p)]);
% 
%     if plot_params.save == 1
%         print(['G:\CurrBiol_revision\figures\figure_4\PID_spikes_histogram_bar_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
%     end
% 
% end

end
