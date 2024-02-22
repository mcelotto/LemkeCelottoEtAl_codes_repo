function [] = plot_PID_Pooled_Sh_newChan_matchedBeh(PIDpath,clusterParams,params,reach_bins,matched_features,plot_params)
% same as compute_PID_Pooled but with new cluster stat as in Combrisson et
% al. (2022) used for channels selection

animal_colors = distinguishable_colors(numel(params.animals));

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = 1:length(params.reachFeatures)
        
        disp(params.reachFeatures{fidx})
        
        early_info_all = [];
        early_info_Orig_all = [];
        early_info_Sh_all = [];
        early_pairs_per_animal = 0;
        late_info_all = [];
        late_info_Orig_all = [];
        late_info_Sh_all = [];
        late_pairs_per_animal = 0;

        tmpFiles = ls([PIDpath 'PID*']);        
        tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
        
        for animal = 1:size(tmpFiles,1)
            
            load([PIDpath tmpFiles(animal,:)])
            
            %%% EARLY
            early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
            early_info_Orig = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
            early_info_Sh = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
            for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
                tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
                tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
                tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                early_info_Orig(pairs,:,:) = tmp_info;
                early_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
                tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                tmp_info(~tmp_sig) = 0;
                early_info(pairs,:,:) = tmp_info;
            end
            early_info_all = cat(1,early_info_all,early_info);
            early_info_Orig_all = cat(1,early_info_Orig_all,early_info_Orig);
            early_info_Sh_all = cat(1,early_info_Sh_all,early_info_Sh);
            early_pairs_per_animal = [early_pairs_per_animal size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)];

            %%% LATE
            late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
            late_info_Orig = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
            late_info_Sh = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
            for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
                tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
                tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
                tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                late_info_Orig(pairs,:,:) = tmp_info;
                late_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
                tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));                
                tmp_info(~tmp_sig) = 0;
                late_info(pairs,:,:) = tmp_info;
            end
            late_info_all = cat(1,late_info_all,late_info);
            late_info_Orig_all = cat(1,late_info_Orig_all,late_info_Orig);
            late_info_Sh_all = cat(1,late_info_Sh_all,late_info_Sh);
            late_pairs_per_animal = [late_pairs_per_animal size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)];

        end

        if fidx==1
            ylim_info = [0 0.01];
            ylim_info_summary = [0 0.1];
            ylim_info_diff = [-.006 .006];
        elseif fidx==2
            ylim_info = [0 0.02];
            ylim_info_summary = [0 0.15];            
            ylim_info_diff = [-.015 .015];
        end
        
        %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)

            figure('units','normalized','outerposition',[0 0 1 1]);

                subplot(2,2,1); hold on;
                    imagesc(squeeze(mean(early_info_all,1))',ylim_info)
                    plot([0.5 39.5],[25.5 25.5],'color','r')
                    plot([0.5 39.5],[26.5 26.5],'color','r')
                    plot([26 26],[.5 51.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 51.5])
                    yticks(1:2:51)
                    yticklabels(-250:20:250)
                    ylabel('DLS time delay (ms)')
                    title(num2str(size(early_info_all,1)));
                subplot(2,2,2); hold on;
                    imagesc(squeeze(mean(late_info_all,1))',ylim_info)
                    plot([0.5 39.5],[25.5 25.5],'color','r')
                    plot([0.5 39.5],[26.5 26.5],'color','r')
                    plot([26 26],[.5 51.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 51.5])
                    yticks(1:2:51)
                    yticklabels(-250:20:250)
                    ylabel('DLS time delay (ms)')
                    title(num2str(size(late_info_all,1)));                    
                clearvars g
                    x = [1:39];
                    y = [mean(early_info_all(:,:,27:51),3); mean(early_info_all(:,:,1:25),3)];
                    c = [ones(1,size(mean(early_info_all(:,:,27:51),3),1)) 2*ones(1,size(mean(early_info_all(:,:,1:25),3),1))];
                    g(2,1)=gramm('x',x,'y',y,'color',c);
                    g(2,1).stat_summary('type','sem');
                    g(2,1).axe_property('XLim',[6 36]);
                    g(2,1).axe_property('YLim',ylim_info);    
                    g(2,1).set_title(['early']);        
                    y = [mean(late_info_all(:,:,27:51),3); mean(late_info_all(:,:,1:25),3)];
                    c = [ones(1,size(mean(late_info_all(:,:,27:51),3),1)) 2*ones(1,size(mean(late_info_all(:,:,1:25),3),1))];
                    g(2,2)=gramm('x',x,'y',y,'color',c);
                    g(2,2).stat_summary('type','sem');
                    g(2,2).axe_property('XLim',[6 36]);
                    g(2,2).axe_property('YLim',ylim_info);    
                    g(2,2).set_title(['late']);        
                    g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
                end
                    
            figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(2,2,1); hold on;
                    imagesc(squeeze(mean(late_info_all,1))'-squeeze(mean(early_info_all,1))',ylim_info_diff)
                    colormap(redblue)
                    plot([0.5 39.5],[25.5 25.5],'color','r')
                    plot([0.5 39.5],[26.5 26.5],'color','r')
                    plot([26 26],[.5 51.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 51.5])
                    yticks(1:2:51)
                    yticklabels(-250:20:250)
                    ylabel('DLS time delay (ms)')
                subplot(2,2,3); hold on;
                    plot(mean(mean(late_info_all(:,:,27:51),3),1)-mean(mean(early_info_all(:,:,27:51),3),1),'color',[.5 1 .25],'LineWidth',2);
                    plot(mean(mean(late_info_all(:,:,1:25),3),1)-mean(mean(early_info_all(:,:,1:25),3),1),'color',[0 .25 0],'LineWidth',2);
                    xlim([6 36])
                    ylim(ylim_info_diff)
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_difference_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
                end
            
        %% HISTOGRAM

            max_early_DLStoM1 = [];
            for pair = 1:size(early_info_all,1)
                max_early_DLStoM1 = [max_early_DLStoM1 max(squeeze(early_info_all(pair,reach_bins,1:25)),[],'all')];
            end
            max_early_M1toDLS = [];
            for pair = 1:size(early_info_all,1)
                max_early_M1toDLS = [max_early_M1toDLS max(squeeze(early_info_all(pair,reach_bins,27:51)),[],'all')];
            end
            max_late_DLStoM1 = [];
            for pair = 1:size(late_info_all,1)
                max_late_DLStoM1 = [max_late_DLStoM1 max(squeeze(late_info_all(pair,reach_bins,1:25)),[],'all')];
            end
            max_late_M1toDLS = [];
            for pair = 1:size(late_info_all,1)
                max_late_M1toDLS = [max_late_M1toDLS max(squeeze(late_info_all(pair,reach_bins,27:51)),[],'all')];
            end

            figure;
                subplot(2,2,1); hold on;
                    histogram(max_early_M1toDLS,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
                    histogram(max_early_DLStoM1,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                    [p,~,stats] = signrank(max_early_M1toDLS,max_early_DLStoM1);
                    title(p);
                    xlim([0 0.1]);
                subplot(2,2,3);
                    hold on;
                    cdfplot(max_early_M1toDLS);
                    cdfplot(max_early_DLStoM1);
                    xlim([0 0.1]);
                    grid off;
                subplot(2,2,2); hold on;
                    histogram(max_late_M1toDLS,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
                    histogram(max_late_DLStoM1,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                    [p,~,stats] = signrank(max_late_M1toDLS,max_late_DLStoM1);
                    title(p);
                    xlim([0 0.1]);
                subplot(2,2,4);
                    hold on;
                    cdfplot(max_late_M1toDLS);
                    cdfplot(max_late_DLStoM1);
                    xlim([0 0.1]);
                    grid off;

                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '_early_late_Histogram.eps'],'-painters','-depsc');
                end   

        %% PLOT DELAYS

            early_peakDelay = cell(1,length(early_pairs_per_animal)-1);
            late_peakDelay = cell(1,length(early_pairs_per_animal)-1);
            
            for animal = 1:length(early_pairs_per_animal)-1
                                
                tmp_peakDelay = [];
                for pair=1+sum(early_pairs_per_animal(1:animal)):sum(early_pairs_per_animal(1:animal+1))
                    [~,time_i] = max(squeeze(max(early_info_all(pair,:,:),[],3)));                    
                    if ismember(time_i,reach_bins)
                        [~, i] = max(squeeze(max(early_info_all(pair,reach_bins,:),[],2)));
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                early_peakDelay{animal} = tmp_peakDelay;
                
                tmp_peakDelay = [];
                for pair = 1+sum(late_pairs_per_animal(1:animal)):sum(late_pairs_per_animal(1:animal+1))
                    [~,time_i] = max(squeeze(max(late_info_all(pair,:,:),[],3)));                    
                    if ismember(time_i,reach_bins)
                        [~, i] = max(squeeze(max(late_info_all(pair,reach_bins,:),[],2)));                        
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                late_peakDelay{animal} = tmp_peakDelay;
                                
            end

            figure;
                subplot(1,2,1); hold on;
                    histogram([early_peakDelay{:}],[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
                    histogram([late_peakDelay{:}],[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
                subplot(1,2,2); hold on;
                    tmp_e = histcounts([early_peakDelay{:}],[1:5:51],'normalization','probability');
                    tmp_l = histcounts([late_peakDelay{:}],[1:5:51],'normalization','probability');
                    scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
                    plot(tmp_l-tmp_e,'color',[0 0 0])
                    plot([0 11],[0 0],'color','k')
                    plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
                    xticks([1:1:10])     
                    xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_delays_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
                end          

            change_pos = [];
            change_neg = [];
            for animal = 1:length(early_pairs_per_animal)-1
                change_pos = [change_pos sum(late_peakDelay{animal}>=27 & late_peakDelay{animal}<=41)/length(late_peakDelay{animal})-sum(early_peakDelay{animal}>=27 & early_peakDelay{animal}<=41)/length(early_peakDelay{animal})];
                change_neg = [change_neg sum(late_peakDelay{animal}>=6 & late_peakDelay{animal}<=21)/length(late_peakDelay{animal})-sum(early_peakDelay{animal}>=6 & early_peakDelay{animal}<=21)/length(early_peakDelay{animal})];
            end         
            
            figure;
                hold on;
                for animal = 1:length(early_pairs_per_animal)-1
                    plot([1 2],[change_neg(animal) change_pos(animal)],'color',animal_colors(animal,:),'LineWidth',1)
                end
                errorbar(1,nanmean(change_neg),nanstd(change_neg)/sqrt(sum(~isnan(change_neg))),'LineWidth',2,'color','k');            
                errorbar(2,nanmean(change_pos),nanstd(change_pos)/sqrt(sum(~isnan(change_pos))),'LineWidth',2,'color','k');
                plot([.5 2.5],[0 0],'color','k','LineStyle','--')
                xlim([.5 2.5])
                [h p1, ~, stats] = ttest(change_pos,change_neg);
                disp(['all animal channel pair change: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                                

                title([params.reachFeatures{fidx} ' | ' num2str(p1)])
            if plot_params.save == 1
                print([plot_params.save_path '\PID_sharedInfo_delays_by_animal_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
            end    


%% OTHER PLOTS
% 
%         %%% SUMMARY STATS
% 
%             max_early_neg = [];
%             for pair = 1:size(early_info_Orig_all,1)
%                 max_early_neg = [max_early_neg max(squeeze(early_info_Orig_all(pair,reach_bins,1:25)),[],'all')];
%             end
%             max_early_neg_Sh = [];
%             for pair = 1:size(early_info_Sh_all,1)
%                 max_early_neg_Sh = [max_early_neg_Sh max(squeeze(early_info_Sh_all(pair,reach_bins,1:25)),[],'all')];
%             end       
%             max_early_pos = [];
%             for pair = 1:size(early_info_Orig_all,1)
%                 max_early_pos = [max_early_pos max(squeeze(early_info_Orig_all(pair,reach_bins,27:51)),[],'all')];
%             end
%             max_early_pos_Sh = [];
%             for pair = 1:size(early_info_Sh_all,1)
%                 max_early_pos_Sh = [max_early_pos_Sh max(squeeze(early_info_Sh_all(pair,reach_bins,27:51)),[],'all')];
%             end
% 
%             max_late_neg = [];
%             for pair = 1:size(late_info_Orig_all,1)
%                 max_late_neg = [max_late_neg max(squeeze(late_info_Orig_all(pair,reach_bins,1:25)),[],'all')];
%             end
%             max_late_neg_Sh = [];
%             for pair = 1:size(late_info_Orig_all,1)
%                 max_late_neg_Sh = [max_late_neg_Sh max(squeeze(late_info_Sh_all(pair,reach_bins,1:25)),[],'all')];
%             end       
%             max_late_pos = [];
%             for pair = 1:size(late_info_Orig_all,1)
%                 max_late_pos = [max_late_pos max(squeeze(late_info_Orig_all(pair,reach_bins,27:51)),[],'all')];
%             end
%             max_late_pos_Sh = [];
%             for pair = 1:size(late_info_Orig_all,1)
%                 max_late_pos_Sh = [max_late_pos_Sh max(squeeze(late_info_Sh_all(pair,reach_bins,27:51)),[],'all')];
%             end
% 
%             figure('units','normalized','outerposition',[0 0 1 1]);
%                 clearvars g    
%                 
%                 Y = [max_early_neg max_early_neg_Sh max_early_pos max_early_pos_Sh];   
%                 X = cell(1,length(Y));
%                 X(1:(length(max_early_neg)+length(max_early_neg_Sh))) = {'DLS -> M1'};
%                 X(1+(length(max_early_neg)+length(max_early_neg_Sh)):end) = {'M1 -> DLS'};
%                 ShReal = [ones(1,length(max_early_neg)) 2*ones(1,length(max_early_neg_Sh)) ones(1,length(max_early_pos)) 2*ones(1,length(max_early_pos_Sh))];        
%                 g(1,1)=gramm('x',X,'y',Y,'color',ShReal); 
%                 g(1,1).stat_boxplot();
%                 [h p1,~,stats] = ttest(max_early_neg,max_early_neg_Sh);
%                 disp(['DLS->M1 naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                
%                 [h p2,~,stats] = ttest(max_early_pos,max_early_pos_Sh);
%                 disp(['M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p2)]);                                
%                 neg_tmp = max_early_neg-max_early_neg_Sh;
%                 pos_tmp = max_early_pos-max_early_pos_Sh;
%                 [h p3,~,stats] = ttest2(neg_tmp,pos_tmp);
%                 disp(['DLS->M1 vs. M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p3)]);                                                
%                 g(1,1).set_title(['M1->DLS p=' num2str(p2) ' | DLS->M1 p=' num2str(p1) ' | Real M1->DLS vs DLS->M1 p=' num2str(p3)]);
%                 g(1,1).axe_property('YLim',ylim_info_summary);
% 
%                 Y = [max_late_neg max_late_neg_Sh max_late_pos max_late_pos_Sh];   
%                 X = cell(1,length(Y));
%                 X(1:(length(max_late_neg)+length(max_late_neg_Sh))) = {'DLS -> M1'};
%                 X(1+(length(max_late_neg)+length(max_late_neg_Sh)):end) = {'M1 -> DLS'};
%                 ShReal = [ones(1,length(max_late_neg)) 2*ones(1,length(max_late_neg_Sh)) ones(1,length(max_late_pos)) 2*ones(1,length(max_late_pos_Sh))];        
%                 g(1,2)=gramm('x',X,'y',Y,'color',ShReal); 
%                 g(1,2).stat_boxplot();
%                 [h p1,~,stats] = ttest(max_late_neg,max_late_neg_Sh);
%                 disp(['DLS->M1 naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                
%                 [h p2,~,stats] = ttest(max_late_pos,max_late_pos_Sh);
%                 disp(['M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p2)]);                                                
%                 neg_tmp = max_late_neg-max_late_neg_Sh;
%                 pos_tmp = max_late_pos-max_late_pos_Sh;
%                 [h p3,~,stats] = ttest2(neg_tmp,pos_tmp);
%                 disp(['DLS->M1 vs. M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p3)]);                                                                
%                 g(1,2).set_title(['M1->DLS p=' num2str(p2) ' | DLS->M1 p=' num2str(p1) ' | Real M1->DLS vs DLS->M1 p=' num2str(p3)]);
%                 g(1,2).axe_property('YLim',ylim_info_summary);    
% 
%                 g.set_title(['EARLY (left) and LATE (right) | ' params.reachFeatures{fidx}]);                    
%                 g.draw();
%                 
%                 if plot_params.save == 1
%                     print([plot_params.save_path '\PID_sharedInfo_summary_stats_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
%                 end    
%                 
%             figure('units','normalized','outerposition',[0 0 1 1]);
%                 subplot(3,2,1); hold on;
%                     imagesc(squeeze(mean(early_info_all,1))',[0 0.01])
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(3,2,2); hold on;
%                     plot(squeeze(mean(mean(early_info_all,1),2))-std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2);
%                     plot(squeeze(mean(mean(early_info_all,1),2))+std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2);
%                     for ny = 5:5:45
%                         if ny == 25
%                             plot([xlim],[ny ny],'color','k','LineWidth',2)
%                         else
%                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
%                         end
%                     end
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(3,2,3); hold on;
%                     imagesc(squeeze(mean(late_info_all,1))',[0 0.01])
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(3,2,4); hold on;
%                     plot(squeeze(mean(mean(late_info_all,1),2))-std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
%                     plot(squeeze(mean(mean(late_info_all,1),2))+std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
%                     for ny = 5:5:45
%                         if ny == 25
%                             plot([xlim],[ny ny],'color','k','LineWidth',2)
%                         else
%                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
%                         end
%                     end
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(3,2,5); hold on;
%                     plot(squeeze(mean(mean(early_info_all,1),2))-std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2);
%                     plot(squeeze(mean(mean(early_info_all,1),2))+std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2); 
%                     plot(squeeze(mean(mean(late_info_all,1),2))-std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
%                     plot(squeeze(mean(mean(late_info_all,1),2))+std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
%                     for ny = 5:5:45
%                         if ny == 25
%                             plot([xlim],[ny ny],'color','k','LineWidth',2)
%                         else
%                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
%                         end
%                     end
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')            
%                 subplot(3,2,6); hold on;
%                     plot(squeeze(mean(mean(late_info_all,1),2))-squeeze(mean(mean(early_info_all,1),2)),[1:51],'color','k','LineWidth',2);
%                     for ny = 5:5:45
%                         if ny == 25
%                             plot([xlim],[ny ny],'color','k','LineWidth',2)
%                         else
%                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
%                         end
%                     end
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')                   

        %%% ORIG VS. SHUFFLED. VS. SIG
        
%             figure; 
%                 subplot(2,3,1); hold on;
%                     imagesc(squeeze(mean(early_info_Orig_all,1))'/size(early_info_Orig_all,1));
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(2,3,2); hold on;
%                     imagesc(squeeze(mean(early_info_Sh_all,1))'/size(early_info_Sh_all,1));
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(2,3,3); hold on;
%                     imagesc(squeeze(mean(early_info_all,1))'/length(find(squeeze(sum(sum(early_info_all,2),3)))))
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(2,3,4); hold on;
%                     imagesc(squeeze(mean(late_info_Orig_all,1))'/size(late_info_Orig_all,1));
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(2,3,5); hold on;
%                     imagesc(squeeze(mean(late_info_Sh_all,1))'/size(late_info_Sh_all,1));
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')
%                 subplot(2,3,6); hold on;
%                     imagesc(squeeze(mean(late_info_all,1))'/length(find(squeeze(sum(sum(late_info_all,2),3)))))
%                     plot([0.5 39.5],[25.5 25.5],'color','r')
%                     plot([0.5 39.5],[26.5 26.5],'color','r')
%                     plot([26 26],[.5 51.5],'color','k')
%                     colorbar;
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                     ylim([.5 51.5])
%                     yticks(1:2:51)
%                     yticklabels(-250:20:250)
%                     ylabel('DLS time delay (ms)')

        %%% COMPARE EARLY/LATE INFO SIG
            
%             early_sig_chans = find(squeeze(sum(sum(early_info_all,2),3)));
%             late_sig_chans = find(squeeze(sum(sum(late_info_all,2),3)));
%             
%             figure
%                 subplot(2,3,[1 2]); hold on;
%                     plot(mean(mean(early_info_all(early_sig_chans,:,27:51),3),1)+(std(mean(early_info_all(early_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(early_info_all(early_sig_chans,:,27:51),3),1)-(std(mean(early_info_all(early_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(late_info_all(late_sig_chans,:,27:51),3),1)+(std(mean(late_info_all(late_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
%                     plot(mean(mean(late_info_all(late_sig_chans,:,27:51),3),1)-(std(mean(late_info_all(late_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
%                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
%                     title(['early/late sig M1 -> DLS shared information | ' params.reachFeatures{fidx}])
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                 subplot(2,3,[4 5]); hold on;
%                     plot(mean(mean(early_info_all(early_sig_chans,:,1:25),3),1)+(std(mean(early_info_all(early_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(early_info_all(early_sig_chans,:,1:25),3),1)-(std(mean(early_info_all(early_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(late_info_all(late_sig_chans,:,1:25),3),1)+(std(mean(late_info_all(late_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
%                     plot(mean(mean(late_info_all(late_sig_chans,:,1:25),3),1)-(std(mean(late_info_all(late_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
%                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
%                     title(['early/late sig DLS -> M1 shared information | ' params.reachFeatures{fidx}])
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                 subplot(2,3,3); hold on;
%                     errorbar(1,mean(max(mean(early_info_all(early_sig_chans,6:36,27:51),3)')),std(max(mean(early_info_all(early_sig_chans,6:36,27:51),3)'))/sqrt(length(early_sig_chans)),'LineWidth',2,'color',[.5 .5 .5])
%                     errorbar(2,mean(max(mean(late_info_all(late_sig_chans,6:36,27:51),3)')),std(max(mean(late_info_all(late_sig_chans,6:36,27:51),3)'))/sqrt(length(late_sig_chans)),'LineWidth',2,'color','b')
%                     [h p] = ttest2(max(mean(early_info_all(early_sig_chans,6:36,27:51),3)'),max(mean(late_info_all(late_sig_chans,6:36,27:51),3)'));
%                     title(num2str(p))
%                     xlim([.5 2.5])
%                 subplot(2,3,6); hold on;
%                     errorbar(1,mean(max(mean(early_info_all(early_sig_chans,6:36,1:25),3)')),std(max(mean(early_info_all(early_sig_chans,6:36,1:25),3)'))/sqrt(length(early_sig_chans)),'LineWidth',2,'color',[.5 .5 .5])
%                     errorbar(2,mean(max(mean(late_info_all(late_sig_chans,6:36,1:25),3)')),std(max(mean(late_info_all(late_sig_chans,6:36,1:25),3)'))/sqrt(length(late_sig_chans)),'LineWidth',2,'color','b')
%                     [h p] = ttest2(max(mean(early_info_all(early_sig_chans,6:36,1:25),3)'),max(mean(late_info_all(late_sig_chans,6:36,1:25),3)'));
%                     title(num2str(p))
%                     xlim([.5 2.5])
               
        %%% COMPARE EARLY/LATE INFO ALL

%             figure
%                 subplot(2,3,[1 2]); hold on;
%                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)+(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)-(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)+(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)-(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
%                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
%                     title(['early/late sig M1 -> DLS shared information | ' params.reachFeatures{fidx}])
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                 subplot(2,3,[4 5]); hold on;
%                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)+(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)-(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)+(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)-(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
%                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
%                     title(['early/late sig DLS -> M1 shared information | ' params.reachFeatures{fidx}])
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                 subplot(2,3,3); hold on;
%                     errorbar(1,mean(max(mean(early_info_Orig_all(:,:,27:51),3)')),std(max(mean(early_info_Orig_all(:,:,27:51),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[.5 .5 .5])
%                     errorbar(2,mean(max(mean(late_info_Orig_all(:,:,27:51),3)')),std(max(mean(late_info_Orig_all(:,:,27:51),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color','b')
%                     [h p] = ttest2(max(mean(early_info_Orig_all(:,:,27:51),3)'),max(mean(late_info_Orig_all(:,:,27:51),3)'));
%                     title(num2str(p))
%                     xlim([.5 2.5])
%                 subplot(2,3,6); hold on;
%                     errorbar(1,mean(max(mean(early_info_Orig_all(:,:,1:25),3)')),std(max(mean(early_info_Orig_all(:,:,1:25),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[.5 .5 .5])
%                     errorbar(2,mean(max(mean(late_info_Orig_all(:,:,1:25),3)')),std(max(mean(late_info_Orig_all(:,:,1:25),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color','b')
%                     [h p] = ttest2(max(mean(early_info_Orig_all(:,:,1:25),3)'),max(mean(late_info_Orig_all(:,:,1:25),3)'));
%                     title(num2str(p))
%                     xlim([.5 2.5])

        %%% COMPARE M1->DLS/DLS->M1 INFO ALL
    
%             figure
%                 subplot(2,3,[1 2]); hold on;
%                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)+(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
%                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)-(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
%                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)+(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
%                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)-(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
%                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
%                     title(['early sig shared information | ' params.reachFeatures{fidx}])
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                 subplot(2,3,[4 5]); hold on;
%                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)+(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)-(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)+(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
%                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)-(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
%                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
%                     title(['late sig shared information | ' params.reachFeatures{fidx}])
%                     xlim([6 36])
%                     xticks([6 16 26 36])
%                     xticklabels([-1 -.5 0 .5])
%                     xlabel('time from pellet touch (ms)')
%                 subplot(2,3,3); hold on;
%                     errorbar(1,mean(max(mean(early_info_Orig_all(:,:,27:51),3)')),std(max(mean(early_info_Orig_all(:,:,27:51),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[.5 1 .25])
%                     errorbar(2,mean(max(mean(early_info_Orig_all(:,:,1:25),3)')),std(max(mean(early_info_Orig_all(:,:,1:25),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[0 .25 0])
%                     [h p] = ttest2(max(mean(early_info_Orig_all(:,:,27:51),3)'),max(mean(early_info_Orig_all(:,:,1:25),3)'));
%                     title(num2str(p))
%                     xlim([.5 2.5])
%                 subplot(2,3,6); hold on;
%                     errorbar(1,mean(max(mean(late_info_Orig_all(:,:,27:51),3)')),std(max(mean(late_info_Orig_all(:,:,27:51),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color',[.5 1 .25])
%                     errorbar(2,mean(max(mean(late_info_Orig_all(:,:,1:25),3)')),std(max(mean(late_info_Orig_all(:,:,1:25),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color',[0 .25 0])
%                     [h p] = ttest2(max(mean(late_info_Orig_all(:,:,27:51),3)'),max(mean(late_info_Orig_all(:,:,1:25),3)'));
%                     title(num2str(p))
%                     xlim([.5 2.5])
       
    end
end

end

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end
