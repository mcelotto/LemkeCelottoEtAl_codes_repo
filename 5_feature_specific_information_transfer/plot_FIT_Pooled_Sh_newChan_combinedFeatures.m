function [] = plot_FIT_Pooled_Sh_newChan_combinedFeatures(FITpath,clusterParams,params,reach_bins,plot_params)

animal_colors = distinguishable_colors(numel(params.animals));

for n_lfp = 1:length(params.lfpFeatures)

    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);

    early_info_all_M1receiver = cell(length(params.reachFeatures),size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(length(params.reachFeatures),size(tmpFiles,1));
    late_info_all_M1receiver = cell(length(params.reachFeatures),size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(length(params.reachFeatures),size(tmpFiles,1));
        
    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        for fidx = 1:length(params.reachFeatures)
            disp(params.reachFeatures{fidx})
            %%% EARLY
                early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);
                for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1)
                    tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver(pairs,:,:));
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiverSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
    %                 tmp_info(~tmp_sig) = 0;
                    early_info_M1receiver(pairs,:,:) = tmp_info;
                    tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver(pairs,:,:));
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
    %                 tmp_info(~tmp_sig) = 0;
                    early_info_DLSreceiver(pairs,:,:) = tmp_info;
                end
                early_info_all_M1receiver{fidx,animal} = early_info_M1receiver;
                early_info_all_DLSreceiver{fidx,animal} = early_info_DLSreceiver;
            %%% LATE
                late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1),39,25);
                late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver,1),39,25);
                for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver,1)
                    tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiver(pairs,:,:));
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_M1receiverSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
    %                 tmp_info(~tmp_sig) = 0;
                    late_info_M1receiver(pairs,:,:) = tmp_info;
                    tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiver(pairs,:,:));
                    tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).FIT_negd_DLSreceiverSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
    %                 tmp_info(~tmp_sig) = 0;
                    late_info_DLSreceiver(pairs,:,:) = tmp_info;
                end
                late_info_all_M1receiver{fidx,animal} = late_info_M1receiver;
                late_info_all_DLSreceiver{fidx,animal} = late_info_DLSreceiver;
        end
    end

    %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)
    
        figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(3,2,1); hold on;
                chan_per_feature = [];
                tmp_infoe_M1receiver = [];
                for n1 = 1:size(early_info_all_M1receiver,1)
                    tmp_length = 0;
                    for n2 = 1:size(early_info_all_M1receiver,2)
                        tmp_infoe_M1receiver = cat(1,tmp_infoe_M1receiver,early_info_all_M1receiver{n1,n2});
                        tmp_length = tmp_length + size(early_info_all_M1receiver{n1,n2},1);
                    end
                    chan_per_feature = [chan_per_feature tmp_length];
                end
                imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 6e-3])
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
                title(['early m1 receiver ' num2str(chan_per_feature)]);
            subplot(3,2,2); hold on;
                chan_per_feature = [];
                tmp_infoe_DLSreceiver = [];
                for n1 = 1:size(early_info_all_DLSreceiver,1)
                    tmp_length = 0;
                    for n2 = 1:size(early_info_all_DLSreceiver,2)
                        tmp_infoe_DLSreceiver = cat(1,tmp_infoe_DLSreceiver,early_info_all_DLSreceiver{n1,n2});
                        tmp_length = tmp_length + size(early_info_all_DLSreceiver{n1,n2},1);
                    end
                    chan_per_feature = [chan_per_feature tmp_length];
                end
                imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 2.25e-3])
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
                title(['early dls receiver ' num2str(chan_per_feature)]);
            subplot(3,2,3); hold on;
                chan_per_feature = [];
                tmp_infol_M1receiver = [];
                for n1 = 1:size(late_info_all_M1receiver,1)
                    tmp_length = 0;
                    for n2 = 1:size(late_info_all_M1receiver,2)
                        tmp_infol_M1receiver = cat(1,tmp_infol_M1receiver,late_info_all_M1receiver{n1,n2});
                        tmp_length = tmp_length + size(late_info_all_M1receiver{n1,n2},1);
                    end
                    chan_per_feature = [chan_per_feature tmp_length];
                end
                imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 6e-3])
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
                title(['late m1 receiver ' num2str(chan_per_feature)]);
            subplot(3,2,4); hold on;
                chan_per_feature = [];
                tmp_infol_DLSreceiver = [];
                for n1 = 1:size(late_info_all_DLSreceiver,1)
                    tmp_length = 0;
                    for n2 = 1:size(late_info_all_DLSreceiver,2)
                        tmp_infol_DLSreceiver = cat(1,tmp_infol_DLSreceiver,late_info_all_DLSreceiver{n1,n2});
                        tmp_length = tmp_length + size(late_info_all_DLSreceiver{n1,n2},1);
                    end
                    chan_per_feature = [chan_per_feature tmp_length];
                end
                imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 2.25e-3])
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
                title(['late dls receiver ' num2str(chan_per_feature)]);
            clearvars g
                x = [1:39];
                y = [mean(tmp_infoe_M1receiver,3); mean(tmp_infol_M1receiver,3)];
                c = [ones(1,size(mean(tmp_infoe_M1receiver,3),1)) 2*ones(1,size(mean(tmp_infol_M1receiver,3),1))];
                g(3,1)=gramm('x',x,'y',y,'color',c);
                g(3,1).stat_summary('type','sem');
                g(3,1).axe_property('XLim',[6 36]);
                g(3,1).axe_property('YLim',[0 0.003]);
                g(3,1).set_title(['M1']);
                y = [mean(tmp_infoe_DLSreceiver,3); mean(tmp_infol_DLSreceiver,3)];
                c = [ones(1,size(mean(tmp_infoe_DLSreceiver,3),1)) 2*ones(1,size(mean(tmp_infol_DLSreceiver,3),1))];
                g(3,2)=gramm('x',x,'y',y,'color',c);
                g(3,2).stat_summary('type','sem');
                g(3,2).axe_property('XLim',[6 36]);
                g(3,2).axe_property('YLim',[0 0.002]);
                g(3,2).set_title(['DLS']);
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\FIT_imagesc.eps'],'-painters','-depsc');
                end

        %% HISTOGRAM

            max_early_DLStoM1 = [];
            for pair = 1:size(tmp_infoe_M1receiver,1)
                max_early_DLStoM1 = [max_early_DLStoM1 max(squeeze(tmp_infoe_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_early_M1toDLS = [];
            for pair = 1:size(tmp_infoe_DLSreceiver,1)
                max_early_M1toDLS = [max_early_M1toDLS max(squeeze(tmp_infoe_DLSreceiver(pair,reach_bins,:)),[],'all')];
            end
            max_late_DLStoM1 = [];
            for pair = 1:size(tmp_infol_M1receiver,1)
                max_late_DLStoM1 = [max_late_DLStoM1 max(squeeze(tmp_infol_M1receiver(pair,reach_bins,:)),[],'all')];
            end
            max_late_M1toDLS = [];
            for pair = 1:size(tmp_infol_DLSreceiver,1)
                max_late_M1toDLS = [max_late_M1toDLS max(squeeze(tmp_infol_DLSreceiver(pair,reach_bins,:)),[],'all')];
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
                    print([plot_params.save_path '\FIT_histogram.eps'],'-painters','-depsc');
                end          
            
    end
end

