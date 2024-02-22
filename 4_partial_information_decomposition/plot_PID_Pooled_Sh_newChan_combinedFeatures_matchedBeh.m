function [] = plot_PID_Pooled_Sh_newChan_combinedFeatures_matchedBeh(PIDpath,clusterParams,params,reach_bins,matched_features,plot_params)
% same as compute_PID_Pooled but with new cluster stat as in Combrisson et
% al. (2022) used for channels selection

animal_colors = distinguishable_colors(numel(params.animals));

for n_lfp = 1:length(params.lfpFeatures)

    for n_matched = 1:length(matched_features)

        tmpFiles = ls([PIDpath '*' matched_features{n_matched} '*']);
        tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
        params.reachFeatures = {'maxVel','distTrav'};

        early_info_all = cell(length(params.reachFeatures),size(tmpFiles,1));
        early_pairs_per_animal = cell(length(params.reachFeatures),size(tmpFiles,1));
        late_info_all = cell(length(params.reachFeatures),size(tmpFiles,1));
        late_pairs_per_animal = cell(length(params.reachFeatures),size(tmpFiles,1));
        for animal = 1:size(tmpFiles,1)
            load([PIDpath tmpFiles(animal,:)])  
            params.reachFeatures = {'maxVel','distTrav'};
            for fidx = 1:length(params.reachFeatures)
                disp(params.reachFeatures{fidx})
                %%% EARLY
                early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
                for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
                    tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
                    tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                    % tmp_info(~tmp_sig) = 0;
                    early_info(pairs,:,:) = tmp_info;
                end
                early_info_all{fidx,animal} = early_info;
                early_pairs_per_animal{fidx,animal} = size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1);
                %%% LATE
                late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
                for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
                    tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
                    tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                    % tmp_info(~tmp_sig) = 0;
                    late_info(pairs,:,:) = tmp_info;
                end
                late_info_all{fidx,animal} = late_info;
                late_pairs_per_animal{fidx,animal} = size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1);
            end
        end

        %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)
    
%         if n_matched==1
%             ylim_feature = 0.018;
%             ylim_feature2 = 0.018;            
%         elseif n_matched==2
%             ylim_feature = 0.018;
%             ylim_feature2 = 0.018;    
%         elseif n_matched==3
%             ylim_feature2 = 0.022;
%             ylim_feature2 = 0.022;
%         end

        figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(2,2,1); hold on;
                chan_per_feature = []; 
                tmp_infoe = [];
                for n1 = 1:size(early_info_all,1)
                    tmp_length = 0;
                    for n2 = 1:size(early_info_all,2)
                        tmp_infoe = cat(1,tmp_infoe,early_info_all{n1,n2});
                        tmp_length = tmp_length + size(early_info_all{n1,n2},1);
                    end
                    chan_per_feature = [chan_per_feature tmp_length];
                end        
                imagesc(squeeze(mean(tmp_infoe,1))',[0 0.016])
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
                title(num2str(chan_per_feature));
            subplot(2,2,2); hold on;
                chan_per_feature = []; 
                tmp_infol = [];
                for n1 = 1:size(late_info_all,1)
                    tmp_length = 0;
                    for n2 = 1:size(late_info_all,2)
                        tmp_infol = cat(1,tmp_infol,late_info_all{n1,n2});
                        tmp_length = tmp_length + size(late_info_all{n1,n2},1);
                    end
                    chan_per_feature = [chan_per_feature tmp_length];
                end      
                imagesc(squeeze(mean(tmp_infol,1))',[0 0.016])
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
                title(num2str(chan_per_feature));
            clearvars g
                x = [1:39];
                y = [mean(tmp_infoe(:,:,27:41),3); mean(tmp_infoe(:,:,11:25),3)];
                c = [ones(1,size(mean(tmp_infoe(:,:,27:41),3),1)) 2*ones(1,size(mean(tmp_infoe(:,:,11:25),3),1))];
                g(2,1)=gramm('x',x,'y',y,'color',c);
                g(2,1).stat_summary('type','sem');
                g(2,1).axe_property('XLim',[6 36]);
                g(2,1).axe_property('YLim',[0 0.016]);    
                g(2,1).set_title(['early']);        
                y = [mean(tmp_infol(:,:,27:41),3); mean(tmp_infol(:,:,11:25),3)];
                c = [ones(1,size(mean(tmp_infol(:,:,27:41),3),1)) 2*ones(1,size(mean(tmp_infol(:,:,11:25),3),1))];
                g(2,2)=gramm('x',x,'y',y,'color',c);
                g(2,2).stat_summary('type','sem');
                g(2,2).axe_property('XLim',[6 36]);
                g(2,2).axe_property('YLim',[0 0.016]);    
                g(2,2).set_title(['late']);        
                g.draw();
            if plot_params.save == 1
                print([plot_params.save_path '\PID_info_' matched_features{n_matched} '.eps'],'-painters','-depsc');
            end

%         figure('units','normalized','outerposition',[0 0 1 1]);
%             subplot(2,2,1); hold on;
%                 imagesc(squeeze(mean(tmp_infol,1))'-squeeze(mean(tmp_infoe,1))',[-.007 .021])
%                 tmp1 = [linspace(0,1,64) ones(1,256-64)];
%                 tmp2 = [linspace(0,1,64) linspace(1,0,256-64)];
%                 tmp3 = [ones(1,64) linspace(1,0,256-64)];
%                 tmp_colormap = [tmp1' tmp2' tmp3'];
%                 colormap(tmp_colormap)
%                 plot([0.5 39.5],[25.5 25.5],'color','r')
%                 plot([0.5 39.5],[26.5 26.5],'color','r')
%                 plot([26 26],[.5 51.5],'color','k')
%                 colorbar;
%                 xlim([6 36])
%                 xticks([6 16 26 36])
%                 xticklabels([-1 -.5 0 .5])
%                 xlabel('time from pellet touch (ms)')
%                 ylim([.5 51.5])
%                 yticks(1:2:51)
%                 yticklabels(-250:20:250)
%                 ylabel('DLS time delay (ms)')
%             subplot(2,2,3); hold on;
%                 plot(mean(mean(tmp_infol(:,:,27:51),3),1)-mean(mean(tmp_infoe(:,:,27:51),3),1),'color',[.5 1 .25],'LineWidth',2);
%                 plot(mean(mean(tmp_infol(:,:,1:25),3),1)-mean(mean(tmp_infoe(:,:,1:25),3),1),'color',[0 .25 0],'LineWidth',2);
%                 xlim([6 36])
%                 ylim([-.015 .015])
%                 xticks([6 16 26 36])
%                 xticklabels([-1 -.5 0 .5])
%                 xlabel('time from pellet touch (ms)')
%             if plot_params.save == 1
%                 print([plot_params.save_path '\PID_info_diff_' matched_features{n_matched} '.eps'],'-painters','-depsc');
%             end
    
        %% CDFs
        
        tmp_infoe = [];
        for n1 = 1:size(early_info_all,1)
            for n2 = 1:size(early_info_all,2)
                tmp_infoe = cat(1,tmp_infoe,early_info_all{n1,n2});
            end
        end
        tmp_infol = [];
        for n1 = 1:size(late_info_all,1)
            for n2 = 1:size(late_info_all,2)
                tmp_infol = cat(1,tmp_infol,late_info_all{n1,n2});
            end
        end
        
        max_early_DLStoM1 = [];
        for pair = 1:size(tmp_infoe,1)
            max_early_DLStoM1 = [max_early_DLStoM1 max(squeeze(tmp_infoe(pair,reach_bins,1:25)),[],'all')];
        end
        max_early_M1toDLS = [];
        for pair = 1:size(tmp_infoe,1)
            max_early_M1toDLS = [max_early_M1toDLS max(squeeze(tmp_infoe(pair,reach_bins,27:51)),[],'all')];
        end
        max_late_DLStoM1 = [];
        for pair = 1:size(tmp_infol,1)
            max_late_DLStoM1 = [max_late_DLStoM1 max(squeeze(tmp_infol(pair,reach_bins,1:25)),[],'all')];
        end
        max_late_M1toDLS = [];
        for pair = 1:size(tmp_infol,1)
            max_late_M1toDLS = [max_late_M1toDLS max(squeeze(tmp_infol(pair,reach_bins,27:51)),[],'all')];
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
            print([plot_params.save_path '\PID_cdfs_' matched_features{n_matched} '.eps'],'-painters','-depsc');
        end
    
        %% PLOT DELAYS
    
        early_peakDelay = cell(length(params.reachFeatures),size(early_info_all,2));
        late_peakDelay = cell(length(params.reachFeatures),size(late_info_all,2));
    
        figure;    
        for fidx = 1:length(params.reachFeatures)
    
            for animal = 1:size(early_info_all,2)
                tmp_peakDelay = [];
                for pair=1:size(early_info_all{fidx,animal},1)
                    [~,time_i] = max(squeeze(max(early_info_all{fidx,animal}(pair,:,:),[],3)));
                    if ismember(time_i,reach_bins)
                        [~, i] = max(squeeze(max(early_info_all{fidx,animal}(pair,reach_bins,:),[],2)));
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                early_peakDelay{fidx,animal} = tmp_peakDelay;
                tmp_peakDelay = [];
                for pair=1:size(late_info_all{fidx,animal},1)
                    [~,time_i] = max(squeeze(max(late_info_all{fidx,animal}(pair,:,:),[],3)));
                    if ismember(time_i,reach_bins)
                        [~, i] = max(squeeze(max(late_info_all{fidx,animal}(pair,reach_bins,:),[],2)));
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                late_peakDelay{fidx,animal} = tmp_peakDelay;
            end
    
            subplot(1,2,1); hold on;
                histogram([early_peakDelay{fidx,:}],[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
                histogram([late_peakDelay{fidx,:}],[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
            subplot(1,2,2); hold on;
                tmp_e = histcounts([early_peakDelay{fidx,:}],[1:5:51],'normalization','probability');
                tmp_l = histcounts([late_peakDelay{fidx,:}],[1:5:51],'normalization','probability');
                scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
                plot(tmp_l-tmp_e,'color',[0 0 0])
                plot([0 11],[0 0],'color','k')
                plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
                xticks([1:1:10])
                xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
        end
        if plot_params.save == 1
            print([plot_params.save_path '\PID_histogram_' matched_features{n_matched} '.eps'],'-painters','-depsc');
        end    
    
        change_pos = [];
        change_neg = [];
        for animal = 1:length(early_pairs_per_animal)-1
            tmp_pos = [];
            tmp_neg = [];
            for fidx = 1:length(params.reachFeatures)
                tmp_pos = [tmp_pos sum(late_peakDelay{fidx,animal}>=27)/sum(~isnan((late_peakDelay{fidx,animal})))-sum(early_peakDelay{fidx,animal}>=27)/sum(~isnan((early_peakDelay{fidx,animal})))];
                tmp_neg = [tmp_neg sum(late_peakDelay{fidx,animal}<=25)/sum(~isnan((late_peakDelay{fidx,animal})))-sum(early_peakDelay{fidx,animal}<=25)/sum(~isnan((early_peakDelay{fidx,animal})))];
            end
            change_pos = [change_pos nanmean(tmp_pos)];
            change_neg = [change_neg nanmean(tmp_neg)];
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
            title([num2str(p1)])
        if plot_params.save == 1
            print([plot_params.save_path '\PID_by_animal_' matched_features{n_matched} '.eps'],'-painters','-depsc');
        end          

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
