function [] = plot_instPID_Pooled_Sh(PIDpath,clusterParams,params,feat2plot,reach_bins,plot_params,recIdx)
% same as compute_PID_Pooled but with new cluster stat as in Combrisson et
% al. (2022) used for channels selection

animal_colors = distinguishable_colors(numel(params.animals));
PID_time = linspace(-1.25,0.65,39); % time points (relative to pellet touch) corresponding to time frames in the PIDout structure
[~,kinFrame] = min(abs(PID_time-recIdx/100));

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = feat2plot
        
        disp(params.instReachFeatures{fidx})
        
        early_info_all = [];
        early_info_Orig_all = [];
        early_info_Sh_all = [];
        early_pairs_per_animal = 0;
        late_info_all = [];
        late_info_Orig_all = [];
        late_info_Sh_all = [];
        late_pairs_per_animal = 0;

        tmpFiles = dir([PIDpath 'PID*']);        
        tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
        
        for animal = 1:size(tmpFiles,1)
            
            load([PIDpath tmpFiles(animal).name],'PIDout')
            
            %%% EARLY
            early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51); % for submitted PID analysis last dim = 51
            early_info_Orig = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            early_info_Sh = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)
                tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared(pairs,:,:));
                tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).sharedSh(pairs,:,:,:));
                tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                early_info_Orig(pairs,:,:) = tmp_info;
                early_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
                tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%                 tmp_info(~tmp_sig) = 0;
                early_info(pairs,:,:) = tmp_info;
            end
            early_info_all = cat(1,early_info_all,early_info);
            early_info_Orig_all = cat(1,early_info_Orig_all,early_info_Orig);
            early_info_Sh_all = cat(1,early_info_Sh_all,early_info_Sh);
            early_pairs_per_animal = [early_pairs_per_animal size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)];

            %%% LATE
            late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            late_info_Orig = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            late_info_Sh = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)
                tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared(pairs,:,:));
                tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).sharedSh(pairs,:,:,:));
                tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                late_info_Orig(pairs,:,:) = tmp_info;
                late_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
                tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));                
%                 tmp_info(~tmp_sig) = 0;
                late_info(pairs,:,:) = tmp_info;
            end
            late_info_all = cat(1,late_info_all,late_info);
            late_info_Orig_all = cat(1,late_info_Orig_all,late_info_Orig);
            late_info_Sh_all = cat(1,late_info_Sh_all,late_info_Sh);
            late_pairs_per_animal = [late_pairs_per_animal size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)];

        end

        if fidx==7
            ylim_info = [0 0.012];
            ylim_info_summary = [0 0.13];
            ylim_info_diff = [-.006 .006];
        elseif fidx==9
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
                    
                    plot([kinFrame kinFrame],[.5 25.5],'--','color',[1,1,1],'linewidth',2)
                    plot([kinFrame kinFrame+4],[25.5 51.5],'--','color',[1,1,1],'linewidth',2)
                    plot([kinFrame kinFrame],[25.5 51.5],'--','color',[227, 113, 227]/255,'linewidth',2)
                    plot([kinFrame+4 kinFrame],[.5 25.5],'--','color',[227, 113, 227]/255,'linewidth',2)

                    colorbar;
                    xlim([11 28])
                    xticks([11 16 21 26]) 
                    xticklabels([-.75 -.5 -.25 0])
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
                    
                    plot([kinFrame kinFrame],[.5 25.5],'--','color',[1,1,1],'linewidth',2)
                    plot([kinFrame kinFrame+4],[25.5 51.5],'--','color',[1,1,1],'linewidth',2)
                    plot([kinFrame kinFrame],[25.5 51.5],'--','color',[227, 113, 227]/255,'linewidth',2)
                    plot([kinFrame+4 kinFrame],[.5 25.5],'--','color',[227, 113, 227]/255,'linewidth',2)

                    colorbar;
                    xlim([11 28])
                    xticks([11 16 21 26]) 
                    xticklabels([-.75 -.5 -.25 0])
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
                    g.set_title(['shared info about ' params.instReachFeatures{fidx} ' | kin lag ' num2str(10*recIdx) 'ms']);

                    g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_' params.instReachFeatures{fidx} '.eps'],'-painters','-depsc');
                end
                    
            figure('units','normalized','outerposition',[0 0 1 1]);
                subplot(2,2,1); hold on;
                    imagesc(squeeze(mean(late_info_all,1))'-squeeze(mean(early_info_all,1))',ylim_info_diff)
                    colormap(redblue)
                    plot([0.5 39.5],[25.5 25.5],'color','r')
                    plot([0.5 39.5],[26.5 26.5],'color','r')
                    plot([26 26],[.5 51.5],'color','k')
                    colorbar;
                    xlim([11 28])
                    xticks([11 16 21 26]) 
                    xticklabels([-.75 -.5 -.25 0])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 51.5])
                    yticks(1:2:51)
                    yticklabels(-250:20:250)
                    ylabel('DLS time delay (ms)')
                subplot(2,2,3); hold on;
                    plot(mean(mean(late_info_all(:,:,27:51),3),1)-mean(mean(early_info_all(:,:,27:51),3),1),'color',[.5 1 .25],'LineWidth',2);
                    plot(mean(mean(late_info_all(:,:,1:25),3),1)-mean(mean(early_info_all(:,:,1:25),3),1),'color',[0 .25 0],'LineWidth',2);
                    xlim([11 28])
                    ylim(ylim_info_diff)
                    xticks([11 16 21 26]) 
                    xticklabels([-.75 -.5 -.25 0])
                    xlabel('time from pellet touch (ms)')
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_difference_' params.instReachFeatures{fidx} '.eps'],'-painters','-depsc');
                end
                
               %%% PLOT Time-Delay max stat (MC July '23)
                if strcmp(plot_params.time_delay,'maxmax')
                    
                    [info_M1_DLS_early,delay_idx_M1_DLS_early] = max(early_info_all(:,:,27:51),[],3);
                    [info_DLS_M1_early,delay_idx_DLS_M1_early] = max(early_info_all(:,:,1:25),[],3);
                    [info_M1_DLS_late,delay_idx_M1_DLS_late] = max(late_info_all(:,:,27:51),[],3);
                    [info_DLS_M1_late,delay_idx_DLS_M1_late] = max(late_info_all(:,:,1:25),[],3);

                    % take maximum over time per pair
                    [peaks_M1_DLS_early,time_idx_M1_DLS_early] = max(info_M1_DLS_early(:,reach_bins),[],2);
                    [peaks_DLS_M1_early,time_idx_DLS_M1_early] = max(info_DLS_M1_early(:,reach_bins),[],2);
                    [peaks_M1_DLS_late,time_idx_M1_DLS_late] = max(info_M1_DLS_late(:,reach_bins),[],2);
                    [peaks_DLS_M1_late,time_idx_DLS_M1_late] = max(info_DLS_M1_late(:,reach_bins),[],2);
                    %sgtitle([params.reachFeatures{fidx},' '])
                elseif strcmp(plot_params.time_delay,'meanmean')
                    peaks_M1_DLS_early = squeeze(mean(mean(early_info_all(:,20:23,27:41),2),3));
                    peaks_DLS_M1_early = squeeze(mean(mean(early_info_all(:,20:23,11:25),2),3));
                    peaks_M1_DLS_late = squeeze(mean(mean(late_info_all(:,23:26,27:41),2),3));
                    peaks_DLS_M1_late = squeeze(mean(mean(late_info_all(:,23:26,11:25),2),3));
                    %sgtitle([params.reachFeatures{fidx},' '])
                elseif strcmp(plot_params.time_delay,'maxDmeanT')
                    peaks_M1_DLS_early = squeeze(mean(max(early_info_all(:,20:23,27:41),[],3),2));
                    peaks_DLS_M1_early = squeeze(mean(max(early_info_all(:,20:23,11:25),[],3),2));
                    peaks_M1_DLS_late = squeeze(mean(max(late_info_all(:,23:26,27:41),[],3),2));
                    peaks_DLS_M1_late = squeeze(mean(max(late_info_all(:,23:26,11:25),[],3),2));
                end
 
                figure()
                pellet_touch_idx = find(reach_bins==26);
                kinFrame_idx = find(reach_bins==kinFrame);
                x_ticks_labels=((1:2:numel(reach_bins))-pellet_touch_idx)*50;
                
                subplot(2,2,1)
                histogram(time_idx_M1_DLS_early,'normalization','probability')
                hold on
                histogram(time_idx_DLS_M1_late,'normalization','probability')
                [p] = ranksum(time_idx_M1_DLS_early,time_idx_DLS_M1_late);
                ylabel('% channel pairs')
                xlim([1,numel(reach_bins)])
                xticks([1:2:numel(reach_bins)]) 
                xticklabels([x_ticks_labels])
                xlabel('time of the receiver [ms]')
                xline(pellet_touch_idx,'k');
                xline(kinFrame_idx,'r');
                title(['Main directions, ' num2str(p)])
                legend('M1->DLS naive','DLS->M1 skilled','time pellet touch','velocity time')

                subplot(2,2,2)
                histogram(time_idx_DLS_M1_early,'normalization','probability')
                hold on
                histogram(time_idx_M1_DLS_late,'normalization','probability')
                [p] = ranksum(time_idx_DLS_M1_early,time_idx_M1_DLS_late);
                ylabel('% channel pairs')
                xlim([1,numel(reach_bins)])
                xticks([1:2:numel(reach_bins)]) 
                xticklabels([x_ticks_labels])
                xlabel('time of the receiver [ms]')
                xline(pellet_touch_idx,'k');
                xline(kinFrame_idx,'r');
                title(['Secondary directions, ' num2str(p)])
                legend('DLS->M1 naive','M1->DLS skilled','time pellet touch','velocity time','Location','northwest')
                
                subplot(2,2,3)
                histogram(time_idx_M1_DLS_early,'normalization','probability')
                hold on
                histogram(time_idx_DLS_M1_early,'normalization','probability')
                [p] = signrank(time_idx_M1_DLS_early,time_idx_DLS_M1_early);
                ylabel('% channel pairs')
                xlim([1,numel(reach_bins)])
                xticks([1:2:numel(reach_bins)]) 
                xticklabels([x_ticks_labels])
                xlabel('time of the receiver [ms]')
                xline(pellet_touch_idx,'k');
                xline(kinFrame_idx,'r');
                title(['Naive, ' num2str(p)])
                legend('M1->DLS','DLS->M1','time pellet touch','velocity time')

                subplot(2,2,4)
                histogram(time_idx_M1_DLS_late,'normalization','probability')
                hold on
                histogram(time_idx_DLS_M1_late,'normalization','probability')
                [p] = signrank(time_idx_M1_DLS_late,time_idx_DLS_M1_late);
                ylabel('% channel pairs')
                xlim([1,numel(reach_bins)])
                xticks([1:2:numel(reach_bins)]) 
                xticklabels([x_ticks_labels])
                xlabel('time of the receiver [ms]')
                xline(pellet_touch_idx,'k');
                xline(kinFrame_idx,'r');
                title(['Skilled, ' num2str(p)])
                legend('M1->DLS','DLS->M1','time pellet touch','velocity time')

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
                    print([plot_params.save_path '\PID_sharedInfo_' params.instReachFeatures{fidx} '_early_late_Histogram.eps'],'-painters','-depsc');
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
            early_delays = [early_peakDelay{:}];
            late_delays = [late_peakDelay{:}];
                subplot(1,2,1); hold on;
                    histogram(early_delays(~isnan(early_delays)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
                    histogram(late_delays(~isnan(late_delays)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
                    [~,p_diff]=ttest2(early_delays,late_delays);
                    title(['p=' num2str(p_diff,2)])
                    
                subplot(1,2,2); hold on;
                    tmp_e = histcounts(early_delays(~isnan(early_delays)),[1:5:51],'normalization','probability');
                    tmp_l = histcounts(late_delays(~isnan(late_delays)),[1:5:51],'normalization','probability');
                    scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
                    plot(tmp_l-tmp_e,'color',[0 0 0])
                    plot([0 11],[0 0],'color','k')
                    plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
                    xticks([1:1:10])     
                    xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_sharedInfo_delays_' params.instReachFeatures{fidx} '.eps'],'-painters','-depsc');
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

                title([params.instReachFeatures{fidx} ' kin lag ' num2str(10*recIdx) 'ms | ' num2str(p1)])
            if plot_params.save == 1
                print([plot_params.save_path '\PID_sharedInfo_delays_by_animal_' params.instReachFeatures{fidx} '.eps'],'-painters','-depsc');
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
