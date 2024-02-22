function [] = plot_mutualInformation_LFP_indPoints(MIout,params,MIsigtiming,clusterParams,f2plot,plot_params)

% sig_chans = cell(numel(params.animals),length(params.lfpFeatures),length(params.reachFeatures));
% for animal = 1:numel(params.animals)
%     for n_lfp = 1:length(params.lfpFeatures)
%         for fidx = 1:length(params.reachFeatures)
%             tmp_sig_chans = cell(length(params.areas),2);
%             for aidx = 1:length(params.areas)
%                 for early_late = 1:2
%                     sig_chan = [];
%                     for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                         
%                         infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
%                         infQuant = temporal_rebinning(infQuant,5,'movmean',2);
%                         
%                         infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
%                         infQuantSh = temporal_rebinning(infQuantSh,5,'movmean',2);
%                         
%                         tmp_sig = clusterStat(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
%                         infQuant = infQuant-mean(squeeze(infQuantSh),2)';                        
%                         infQuant(~tmp_sig) = 0;
%                         [~,i] = max(infQuant);
% 
%                         sigTime_start = find(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==MIsigtiming(1));
%                         sigTime_end = find(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==MIsigtiming(2));
%                         sig_chan = [sig_chan (i>=sigTime_start && i<=sigTime_end)];
%                     end
%                     
%                     tmp_sig_chans{aidx,early_late} = sig_chan;
%                 end
%             end
%             sig_chans{animal,n_lfp,fidx} = tmp_sig_chans;
%         end
%     end
% end

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = f2plot
        
        info_mean_sig.M1 = cell(1,2);
        info_mean_sig.DLS = cell(1,2);
        for animal = 1:numel(params.animals)            
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                         if sig_chans{animal,n_lfp,fidx}{aidx,early_late}(chan)==1      
                            
                            infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                            infQuant = temporal_rebinning(infQuant,5,'movmean',2);
                            infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                            infQuantSh = temporal_rebinning(infQuantSh,5,'movmean',2);
                            
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_mean_sig.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
                            
%                         end
                    end
                end
            end
        end

        load('F:\final_IIT\data\processed_data\realBinTimes.mat');
        timeBins4LFP = round(newBinTimes/10)+150;

        m1_early_info = [];
        m1_early_info_all = [];
        for unit = 1:size(info_mean_sig.M1{1},1)
            tmp_info = info_mean_sig.M1{1}(unit,:);
            m1_early_info_all = [m1_early_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info];
            end
        end

        m1_late_info = [];
        m1_late_info_all = [];
        for unit = 1:size(info_mean_sig.M1{2},1)
            tmp_info = info_mean_sig.M1{2}(unit,:);
            m1_late_info_all = [m1_late_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info];
            end
        end

        dls_early_info = [];
        dls_early_info_all = [];
        for unit = 1:size(info_mean_sig.DLS{1},1)
            tmp_info = info_mean_sig.DLS{1}(unit,:);
            dls_early_info_all = [dls_early_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info];
            end
        end

        dls_late_info = [];
        dls_late_info_all = [];
        for unit = 1:size(info_mean_sig.DLS{2},1)
            tmp_info = info_mean_sig.DLS{2}(unit,:);
            dls_late_info_all = [dls_late_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info];
            end
        end

        disp(['LFP | ' params.reachFeatures{fidx} ' | M1 | early: ' num2str(size(m1_early_info,1)) ' out of ' num2str(size(m1_early_info_all,1)) ' - ' num2str(100*size(m1_early_info,1)/size(m1_early_info_all,1))  ...
                                                      '% | late: ' num2str(size(m1_late_info,1)) ' out of ' num2str(size(m1_late_info_all,1)) ' - ' num2str(100*size(m1_late_info,1)/size(m1_late_info_all,1))  '%'])
        disp(['LFP | ' params.reachFeatures{fidx} ' | DLS | early: ' num2str(size(dls_early_info,1)) ' out of ' num2str(size(dls_early_info_all,1)) ' - ' num2str(100*size(dls_early_info,1)/size(dls_early_info_all,1))  ...
                                                      '% | late: ' num2str(size(dls_late_info,1)) ' out of ' num2str(size(dls_late_info_all,1)) ' - ' num2str(100*size(dls_late_info,1)/size(dls_late_info_all,1))  '%'])

    %%% PLOT INFO CASCADES AND MEAN / SEM

        figure('units','normalized','outerposition',[0 0 1 1])
            subplot(3,2,1); hold on;
                [~, i_early_info] = max(m1_early_info');
                [~, i_early_info] = sort(i_early_info);
                imagesc([m1_early_info(i_early_info,:)],[0 .15]);
                colorbar;
                plot([1 98],[length(i_early_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_early_info)])                
            subplot(3,2,3); hold on;
                [~, i_late_info] = max(m1_late_info');
                [~, i_late_info] = sort(i_late_info);
                imagesc([m1_late_info(i_late_info,:)],[0 .15]);
                plot([1 98],[length(i_late_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_late_info)])           
            subplot(3,2,2); hold on;
                [~, i_early_info] = max(dls_early_info');
                [~, i_early_info] = sort(i_early_info);
                imagesc([dls_early_info(i_early_info,:)],[0 .15]);
                plot([1 98],[length(i_early_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_early_info)])           
            subplot(3,2,4); hold on;
                [~, i_late_info] = max(dls_late_info');
                [~, i_late_info] = sort(i_late_info);
                imagesc([dls_late_info(i_late_info,:)],[0 .15]);
                plot([1 98],[length(i_late_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_late_info)])           
            subplot(3,2,5); hold on;
                plot(nanmean(m1_early_info)+nanstd(m1_early_info)/sqrt(size(m1_early_info,1)),'k')
                plot(nanmean(m1_early_info)-nanstd(m1_early_info)/sqrt(size(m1_early_info,1)),'k')   
                plot(nanmean(m1_late_info)+nanstd(m1_late_info)/sqrt(size(m1_late_info,1)),'r')
                plot(nanmean(m1_late_info)-nanstd(m1_late_info)/sqrt(size(m1_late_info,1)),'r')
                ylim([0 0.04])
                xlim([12 86])
            subplot(3,2,6); hold on;
                plot(nanmean(dls_early_info)+nanstd(dls_early_info)/sqrt(size(dls_early_info,1)),'k')
                plot(nanmean(dls_early_info)-nanstd(dls_early_info)/sqrt(size(dls_early_info,1)),'k')   
                plot(nanmean(dls_late_info)+nanstd(dls_late_info)/sqrt(size(dls_late_info,1)),'r')
                plot(nanmean(dls_late_info)-nanstd(dls_late_info)/sqrt(size(dls_late_info,1)),'r')
                ylim([0 0.04])
                xlim([12 86])

        figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(m1_early_info); abs(m1_late_info)];
                c = [ones(1,size(m1_early_info,1)) 2*ones(1,size(m1_late_info,1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 .05]);           
                g(1,1).set_title('m1 lfp information');        
                y = [abs(dls_early_info); abs(dls_late_info)];
                c = [ones(1,size(dls_early_info,1)) 2*ones(1,size(dls_late_info,1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 .05]);     
                g(1,2).set_title('dls lfp information');        
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_mean_sem_indPoints.eps'],'-painters','-depsc');
                end     

    %%% PEAK SHIFT INFO

        edge_m1_early = [];
        for unit = 1:size(m1_early_info,1)
            baseSD = 2*std(m1_early_info(unit,:));
            baseMean = mean(m1_early_info(unit,:));
            edge_m1_early = [edge_m1_early min(find(m1_early_info(unit,:)>baseMean+baseSD))];
        end
        edge_m1_early = edge_m1_early(edge_m1_early>33 & edge_m1_early<85);      

        edge_m1_late = [];
        for unit = 1:size(m1_late_info,1)
            baseSD = 2*std(m1_late_info(unit,:));
            baseMean = mean(m1_late_info(unit,:));
            edge_m1_late = [edge_m1_late min(find(m1_late_info(unit,:)>baseMean+baseSD))];
        end
        edge_m1_late = edge_m1_late(edge_m1_late>33 & edge_m1_late<85);     

        edge_dls_early = [];
        for unit = 1:size(dls_early_info,1)
            baseSD = 2*std(dls_early_info(unit,:));
            baseMean = mean(dls_early_info(unit,:));
            edge_dls_early = [edge_dls_early min(find(dls_early_info(unit,:)>baseMean+baseSD))];
        end
        edge_dls_early = edge_dls_early(edge_dls_early>33 & edge_dls_early<85);     

        edge_dls_late = [];
        for unit = 1:size(dls_late_info,1)
            baseSD = 2*std(dls_late_info(unit,:));
            baseMean = mean(dls_late_info(unit,:));
            edge_dls_late = [edge_dls_late min(find(dls_late_info(unit,:)>baseMean+baseSD))];
        end
        edge_dls_late = edge_dls_late(edge_dls_late>33 & edge_dls_late<85);                 

        figure('units','normalized','outerposition',[0 0 1 1]);
            clearvars g;
            peak1_i = edge_m1_early;
            peak2_i = edge_m1_late;
            grouping = [ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))];
            ID = cell(1,length([ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))]));
            ID(1:length(peak1_i)) = {'Naive'};
            ID(1+length(peak1_i):length(grouping)) = {'Learned'};
            g(1,1)=gramm('x',ID,'y',[peak1_i peak2_i],'color',grouping);
            g(1,1).stat_boxplot();
            g(1,1).coord_flip();
            g(1,1).axe_property('YLim',[12 86]);
            [h, p, ci, stats] = ttest2(peak1_i,peak2_i);
            [p2, ~, ~] = ranksum(peak1_i,peak2_i);
            g.set_title(['m1 shift | ' params.reachFeatures{fidx} ' | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
            g.draw();
            if plot_params.save == 1
                print([plot_params.save_path '\LFP_information_M1_' params.reachFeatures{fidx} '_timingShift_indPoints.eps'],'-painters','-depsc');
            end
    
        figure('units','normalized','outerposition',[0 0 1 1]);
            clearvars g;
            peak1_i = edge_dls_early;
            peak2_i = edge_dls_late;
            grouping = [ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))];
            ID = cell(1,length([ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))]));
            ID(1:length(peak1_i)) = {'Naive'};
            ID(1+length(peak1_i):length(grouping)) = {'Learned'};
            g(1,1)=gramm('x',ID,'y',[peak1_i peak2_i],'color',grouping);
            g(1,1).stat_boxplot();
            g(1,1).coord_flip();
            g(1,1).axe_property('YLim',[12 86]);
            [h, p, ci, stats] = ttest2(peak1_i,peak2_i);
            [p2, ~, ~] = ranksum(peak1_i,peak2_i);
            g.set_title(['dls shift | ' params.reachFeatures{fidx} ' | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
            g.draw();
            if plot_params.save == 1
                print([plot_params.save_path '\LFP_information_DLS_' params.reachFeatures{fidx} '_timingShift_indPoints.eps'],'-painters','-depsc');
            end

    %%% NEW NAIVE SKILLED COMPARISON HISTOGRAM/CDF

        figure;
            subplot(2,2,1); hold on;
                histogram(max(m1_early_info_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(m1_late_info_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(m1_early_info_all'),max(m1_late_info_all'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,3);
                hold on;
                cdfplot(max(m1_early_info_all'));
                cdfplot(max(m1_late_info_all'));
                xlim([0 0.25]);
                grid off;
            subplot(2,2,2); hold on;
                histogram(max(dls_early_info_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(dls_late_info_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(dls_early_info_all'),max(dls_late_info_all'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,4);
                hold on;
                cdfplot(max(dls_early_info_all'));
                cdfplot(max(dls_late_info_all'));
                xlim([0 0.25]);
                grid off;
            if plot_params.save == 1
                print([plot_params.save_path '\LFP_information_new_skilled_' params.reachFeatures{fidx} '_indPoints.eps'],'-painters','-depsc');
            end

    end
end

end