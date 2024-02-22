function [] = plot_FIT_indPoints_Sh(FITpath,clusterParams,params,reach_bins,plot_params)

animal_colors = distinguishable_colors(numel(params.animals));

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = 1:length(params.reachFeatures)
        
        early_info_all_M1receiver = [];
        early_info_Orig_all_M1receiver = [];
        early_info_Sh_all_M1receiver = [];
        early_info_all_DLSreceiver = [];
        early_info_Orig_all_DLSreceiver = [];
        early_info_Sh_all_DLSreceiver = [];
        early_pairs_per_animal = 0;
        
        late_info_all_M1receiver = [];
        late_info_Orig_all_M1receiver = [];
        late_info_Sh_all_M1receiver = [];
        late_info_all_DLSreceiver = [];
        late_info_Orig_all_DLSreceiver = [];
        late_info_Sh_all_DLSreceiver = [];
        late_pairs_per_animal = 0;

        tmpFiles = ls([FITpath 'FIT*']);        
        tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
        
        for animal = 1:size(tmpFiles,1)
            
            load([FITpath tmpFiles(animal,:)])
            
            %%% EARLY
            
                early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                early_info_Orig_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                early_info_Sh_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);
                early_info_Orig_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);
                early_info_Sh_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);

                for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1)

                    tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver(pairs,:,:));
                    tmp_info = temporal_rebinning(tmp_info',10,'movmean',5);
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiverSh(pairs,:,:,:));
                    tmp_infoSh = temporal_rebinning(permute(tmp_infoSh,[2 1 3]),10,'movmean',5);
                    tmp_sig = clusterStat(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    early_info_Orig_M1receiver(pairs,:,:) = tmp_info';
                    early_info_Sh_M1receiver(pairs,:,:) = squeeze(mean(tmp_infoSh,3))';
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                    tmp_info(~tmp_sig) = 0;
                    early_info_M1receiver(pairs,:,:) = tmp_info';

                    tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver(pairs,:,:));
                    tmp_info = temporal_rebinning(tmp_info',10,'movmean',5);                    
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh(pairs,:,:,:));
                    tmp_infoSh = temporal_rebinning(permute(tmp_infoSh,[2 1 3]),10,'movmean',5);                    
                    tmp_sig = clusterStat(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    early_info_Orig_DLSreceiver(pairs,:,:) = tmp_info';
                    early_info_Sh_DLSreceiver(pairs,:,:) = squeeze(mean(tmp_infoSh,3))';
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                    tmp_info(~tmp_sig) = 0;
                    early_info_DLSreceiver(pairs,:,:) = tmp_info';

                end

                early_info_all_M1receiver = cat(1,early_info_all_M1receiver,early_info_M1receiver);
                early_info_Orig_all_M1receiver = cat(1,early_info_Orig_all_M1receiver,early_info_Orig_M1receiver);
                early_info_Sh_all_M1receiver = cat(1,early_info_Sh_all_M1receiver,early_info_Sh_M1receiver);
                early_info_all_DLSreceiver = cat(1,early_info_all_DLSreceiver,early_info_DLSreceiver);
                early_info_Orig_all_DLSreceiver = cat(1,early_info_Orig_all_DLSreceiver,early_info_Orig_DLSreceiver);
                early_info_Sh_all_DLSreceiver = cat(1,early_info_Sh_all_DLSreceiver,early_info_Sh_DLSreceiver);
                early_pairs_per_animal = [early_pairs_per_animal size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1)];

            %%% LATE

                late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                late_info_Orig_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                late_info_Sh_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);
                late_info_Orig_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);
                late_info_Sh_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);

                for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1)

                    tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver(pairs,:,:));
                    tmp_info = temporal_rebinning(tmp_info',10,'movmean',5);
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiverSh(pairs,:,:,:));
                    tmp_infoSh = temporal_rebinning(permute(tmp_infoSh,[2 1 3]),10,'movmean',5);
                    tmp_sig = clusterStat(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    late_info_Orig_M1receiver(pairs,:,:) = tmp_info';
                    late_info_Sh_M1receiver(pairs,:,:) = squeeze(mean(tmp_infoSh,3))';
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                    tmp_info(~tmp_sig) = 0;
                    late_info_M1receiver(pairs,:,:) = tmp_info';

                    tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver(pairs,:,:));
                    tmp_info = temporal_rebinning(tmp_info',10,'movmean',5);
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh(pairs,:,:,:));
                    tmp_infoSh = temporal_rebinning(permute(tmp_infoSh,[2 1 3]),10,'movmean',5);
                    tmp_sig = clusterStat(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    late_info_Orig_DLSreceiver(pairs,:,:) = tmp_info';
                    late_info_Sh_DLSreceiver(pairs,:,:) = squeeze(mean(tmp_infoSh,3))';
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                    tmp_info(~tmp_sig) = 0;
                    late_info_DLSreceiver(pairs,:,:) = tmp_info';

                end

                late_info_all_M1receiver = cat(1,late_info_all_M1receiver,late_info_M1receiver);
                late_info_Orig_all_M1receiver = cat(1,late_info_Orig_all_M1receiver,late_info_Orig_M1receiver);
                late_info_Sh_all_M1receiver = cat(1,late_info_Sh_all_M1receiver,late_info_Sh_M1receiver);
                late_info_all_DLSreceiver = cat(1,late_info_all_DLSreceiver,late_info_DLSreceiver);
                late_info_Orig_all_DLSreceiver = cat(1,late_info_Orig_all_DLSreceiver,late_info_Orig_DLSreceiver);
                late_info_Sh_all_DLSreceiver = cat(1,late_info_Sh_all_DLSreceiver,late_info_Sh_DLSreceiver);
                late_pairs_per_animal = [late_pairs_per_animal size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1)];

        end
        
        if fidx==1
            ylim_info_M1receiver = [0 0.0025];
            ylim_info_DLSreceiver = [0 0.0015];
            ylim_info1 = [0 0.001];
            ylim_info2 = [0 0.001];
        elseif fidx==2
            ylim_info_M1receiver = [0 0.005];
            ylim_info_DLSreceiver = [0 0.0025];
            ylim_info1 = [0 0.0015];
            ylim_info2 = [0 0.003];
        end
        
        %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)

            figure('units','normalized','outerposition',[0 0 1 1]);

                subplot(3,2,1); hold on;
                    imagesc(squeeze(mean(early_info_all_M1receiver,1))',ylim_info_M1receiver)
                    plot([26 26],[.5 25.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 25.5])
                    yticks(1:2:25)
                    yticklabels(-250:20:-10)
                    ylabel('DLS time delay (ms)')
                    title(num2str(size(early_info_all_M1receiver,1)));
                subplot(3,2,2); hold on;
                    imagesc(squeeze(mean(late_info_all_M1receiver,1))',ylim_info_M1receiver)
                    plot([26 26],[.5 25.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 25.5])
                    yticks(1:2:25)
                    yticklabels(-250:20:-10)
                    ylabel('DLS time delay (ms)')
                    title(num2str(size(late_info_all_M1receiver,1))); 

                subplot(3,2,3); hold on;
                    imagesc(squeeze(mean(early_info_all_DLSreceiver,1))',ylim_info_DLSreceiver)
                    plot([26 26],[.5 25.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 25.5])
                    yticks(1:2:25)
                    yticklabels(-250:20:-10)
                    ylabel('M1 time delay (ms)')
                    title(num2str(size(early_info_all_DLSreceiver,1)));
                subplot(3,2,4); hold on;
                    imagesc(squeeze(mean(late_info_all_DLSreceiver,1))',ylim_info_DLSreceiver)
                    plot([26 26],[.5 25.5],'color','k')
                    colorbar;
                    xlim([6 36])
                    xticks([6 16 26 36])
                    xticklabels([-1 -.5 0 .5])
                    xlabel('time from pellet touch (ms)')
                    ylim([.5 25.5])
                    yticks(1:2:25)
                    yticklabels(-250:20:-10)
                    ylabel('M1 time delay (ms)')
                    title(num2str(size(late_info_all_DLSreceiver,1))); 

                    clearvars g
                    x = [1:39];
                    y = [mean(early_info_all_M1receiver,3); mean(early_info_all_DLSreceiver,3)];
                    c = [ones(1,size(mean(early_info_all_M1receiver,3),1)) 2*ones(1,size(mean(early_info_all_DLSreceiver,3),1))];
                    g(3,1)=gramm('x',x,'y',y,'color',c);
                    g(3,1).stat_summary('type','sem');
                    g(3,1).axe_property('XLim',[6 36]);
                    g(3,1).axe_property('YLim',ylim_info1);
                    g(3,1).set_title(['early']);
                    y = [mean(late_info_all_M1receiver,3); mean(late_info_all_DLSreceiver,3)];
                    c = [ones(1,size(mean(late_info_all_M1receiver,3),1)) 2*ones(1,size(mean(late_info_all_DLSreceiver,3),1))];
                    g(3,2)=gramm('x',x,'y',y,'color',c);
                    g(3,2).stat_summary('type','sem');
                    g(3,2).axe_property('XLim',[6 36]);
                    g(3,2).axe_property('YLim',ylim_info2);
                    g(3,2).set_title(['late']);
                    g.draw();
                    if plot_params.save == 1
                        print([plot_params.save_path '\FIT_sharedInfo_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
                    end
 
        %% HISTOGRAM

            max_early_DLStoM1 = [];
            for pair = 1:size(early_info_all_M1receiver,1)
                max_early_DLStoM1 = [max_early_DLStoM1 max(squeeze(early_info_all_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_early_M1toDLS = [];
            for pair = 1:size(early_info_all_DLSreceiver,1)
                max_early_M1toDLS = [max_early_M1toDLS max(squeeze(early_info_all_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end
            max_late_DLStoM1 = [];
            for pair = 1:size(late_info_all_M1receiver,1)
                max_late_DLStoM1 = [max_late_DLStoM1 max(squeeze(late_info_all_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_late_M1toDLS = [];
            for pair = 1:size(late_info_all_DLSreceiver,1)
                max_late_M1toDLS = [max_late_M1toDLS max(squeeze(late_info_all_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end

            figure;
                subplot(2,2,1); hold on;
                    histogram(max_early_M1toDLS,[0:0.005:.05],'Normalization','Probability','DisplayStyle','Stairs');    
                    histogram(max_late_M1toDLS,[0:0.005:.05],'Normalization','Probability','DisplayStyle','Stairs');
                    [p,~,stats] = ranksum(max_early_M1toDLS,max_late_M1toDLS);
                    title(p);
                    xlim([0 0.05]);
                subplot(2,2,3);
                    hold on;
                    cdfplot(max_early_M1toDLS);
                    cdfplot(max_late_M1toDLS);
                    xlim([0 0.05]);
                    grid off;
                subplot(2,2,2); hold on;
                    histogram(max_early_DLStoM1,[0:0.005:.05],'Normalization','Probability','DisplayStyle','Stairs');    
                    histogram(max_late_DLStoM1,[0:0.005:.05],'Normalization','Probability','DisplayStyle','Stairs');
                    [p,~,stats] = ranksum(max_early_DLStoM1,max_late_DLStoM1);
                    title(p);
                    xlim([0 0.05]);
                subplot(2,2,4);
                    hold on;
                    cdfplot(max_early_DLStoM1);
                    cdfplot(max_late_DLStoM1);
                    xlim([0 0.05]);
                    grid off;

                if plot_params.save == 1
                    print([plot_params.save_path '\FIT_sharedInfo_' params.reachFeatures{fidx} '_early_late_Histogram.eps'],'-painters','-depsc');
                end          
                    
        %% PLOT DELAYS

            early_peakDelay = cell(1,length(early_pairs_per_animal)-1);
            late_peakDelay = cell(1,length(early_pairs_per_animal)-1);
            
            for animal = 1:length(early_pairs_per_animal)-1
                                
                tmp_peakDelay = [];
                for pair=1+sum(early_pairs_per_animal(1:animal)):sum(early_pairs_per_animal(1:animal+1))                    
                    tmp_info = squeeze(early_info_all_M1receiver(pair,:,:))';
                    tmp_info2 = squeeze(early_info_all_DLSreceiver(pair,:,:))';
                    tmp_info = [tmp_info; tmp_info2(flip(1:25),:)];
                    [~,time_i] = max(max(tmp_info,[],1));                    
                    if ismember(time_i,reach_bins)
                        [~, i] = max(max(tmp_info,[],2));
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                early_peakDelay{animal} = tmp_peakDelay;
                
                tmp_peakDelay = [];
                for pair=1+sum(late_pairs_per_animal(1:animal)):sum(late_pairs_per_animal(1:animal+1))                    
                    tmp_info = squeeze(late_info_all_M1receiver(pair,:,:))';
                    tmp_info2 = squeeze(late_info_all_DLSreceiver(pair,:,:))';
                    tmp_info = [tmp_info; tmp_info2(flip(1:25),:)];
                    [~,time_i] = max(max(tmp_info,[],1));                    
                    if ismember(time_i,reach_bins)
                        [~, i] = max(max(tmp_info,[],2));
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
%                 if plot_params.save == 1
%                     print([plot_params.save_path '\FIT_sharedInfo_delays_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
%                 end          

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
                [h p1] = ttest(change_pos,change_neg)
                title([params.reachFeatures{fidx} ' | ' num2str(p1)])
%             if plot_params.save == 1
%                 print([plot_params.save_path '\FIT_sharedInfo_delays_by_animal_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
%             end    
               
        %% SUMMARY STATS

            max_early_M1receiver = [];
            for pair = 1:size(early_info_Orig_all_M1receiver,1)
                max_early_M1receiver = [max_early_M1receiver max(squeeze(early_info_Orig_all_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_early_DLSreceiver = [];
            for pair = 1:size(early_info_Orig_all_DLSreceiver,1)
                max_early_DLSreceiver = [max_early_DLSreceiver max(squeeze(early_info_Orig_all_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end            
            max_late_M1receiver = [];
            for pair = 1:size(late_info_Orig_all_M1receiver,1)
                max_late_M1receiver = [max_late_M1receiver max(squeeze(late_info_Orig_all_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_late_DLSreceiver = [];
            for pair = 1:size(late_info_Orig_all_DLSreceiver,1)
                max_late_DLSreceiver = [max_late_DLSreceiver max(squeeze(late_info_Orig_all_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end    
            
            max_early_M1receiver_Sh = [];
            for pair = 1:size(early_info_Sh_all_M1receiver,1)
                max_early_M1receiver_Sh = [max_early_M1receiver_Sh max(squeeze(early_info_Sh_all_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_early_DLSreceiver_Sh = [];
            for pair = 1:size(early_info_Sh_all_DLSreceiver,1)
                max_early_DLSreceiver_Sh = [max_early_DLSreceiver_Sh max(squeeze(early_info_Sh_all_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end            
            max_late_M1receiver_Sh = [];
            for pair = 1:size(late_info_Sh_all_M1receiver,1)
                max_late_M1receiver_Sh = [max_late_M1receiver_Sh max(squeeze(late_info_Sh_all_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_late_DLSreceiver_Sh = [];
            for pair = 1:size(late_info_Sh_all_DLSreceiver,1)
                max_late_DLSreceiver_Sh = [max_late_DLSreceiver_Sh max(squeeze(late_info_Sh_all_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end    

            figure('units','normalized','outerposition',[0 0 1 1]);
                clearvars g    
                
                Y = [max_early_M1receiver max_early_M1receiver_Sh max_late_M1receiver max_late_M1receiver_Sh];   
                X = cell(1,length(Y));
                X(1:(length(max_early_M1receiver)+length(max_early_M1receiver_Sh))) = {'early M1 receiver'};
                X(1+(length(max_early_M1receiver)+length(max_early_M1receiver_Sh)):end) = {'late M1 receiver'};
                ShReal = [ones(1,length(max_early_M1receiver)) 2*ones(1,length(max_early_M1receiver_Sh)) ones(1,length(max_late_M1receiver)) 2*ones(1,length(max_late_M1receiver_Sh))];        
                g(1,1)=gramm('x',X,'y',Y,'color',ShReal); 
                g(1,1).stat_boxplot();
                [h p1] = ttest(max_early_M1receiver,max_early_M1receiver_Sh)
                [h p2] = ttest(max_late_M1receiver,max_late_M1receiver_Sh)
                early_tmp = max_early_M1receiver-max_early_M1receiver_Sh;
                late_tmp = max_late_M1receiver-max_late_M1receiver_Sh;
                [h p3,~,stats] = ttest2(early_tmp,late_tmp)
                g(1,1).set_title(['late M1 receiver p=' num2str(p2) ' | early M1 receiver p=' num2str(p1) ' | Real early vs late M1 receiver p=' num2str(p3)]);
                g(1,1).axe_property('YLim',[0 0.05]);

                Y = [max_early_DLSreceiver max_early_DLSreceiver_Sh max_late_DLSreceiver max_late_DLSreceiver_Sh];   
                X = cell(1,length(Y));
                X(1:(length(max_early_DLSreceiver)+length(max_early_DLSreceiver_Sh))) = {'early DLS receiver'};
                X(1+(length(max_early_DLSreceiver)+length(max_early_DLSreceiver_Sh)):end) = {'late DLS receiver'};
                ShReal = [ones(1,length(max_early_DLSreceiver)) 2*ones(1,length(max_early_DLSreceiver_Sh)) ones(1,length(max_late_DLSreceiver)) 2*ones(1,length(max_late_DLSreceiver_Sh))];        
                g(1,2)=gramm('x',X,'y',Y,'color',ShReal); 
                g(1,2).stat_boxplot();
                [h p1] = ttest(max_early_DLSreceiver,max_early_DLSreceiver_Sh)
                [h p2] = ttest(max_late_DLSreceiver,max_late_DLSreceiver_Sh)
                early_tmp = max_early_DLSreceiver-max_early_DLSreceiver_Sh;
                late_tmp = max_late_DLSreceiver-max_late_DLSreceiver_Sh;
                [h p3,~,stats] = ttest2(early_tmp,late_tmp)
                g(1,2).set_title(['late DLS receiver p=' num2str(p2) ' | early DLS receiver p=' num2str(p1) ' | Real early vs late DLS receiver p=' num2str(p3)]);
                g(1,2).axe_property('YLim',[0 0.05]);

                g.set_title(['M1 (left) and DLS(right) | ' params.reachFeatures{fidx}])                    
                g.draw();
                
                if plot_params.save == 1
                    print([plot_params.save_path '\FIT_sharedInfo_summary_stats_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
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
