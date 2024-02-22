function [] = plot_mutualInformation_successFail(MIsuccessFail,MIorig,clusterParams,f2plot,save_params,params,newBinTimes)

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = f2plot
        pause(0.1);
    
    %%% LOAD DATA

        info_orig.M1 = cell(1,2);
        info_orig.DLS = cell(1,2);
        info_success.M1 = cell(1,2);
        info_success.DLS = cell(1,2);
        info_fail.M1 = cell(1,2);
        info_fail.DLS = cell(1,2);        
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIorig.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)

                        infQuant = MIorig.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MIorig.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_orig.(params.areas{aidx}){early_late} = [info_orig.(params.areas{aidx}){early_late}; infQuant];

                        infQuant = MIsuccessFail.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess(chan,:);
                        infQuantSh = MIsuccessFail.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoShSuccess(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_success.(params.areas{aidx}){early_late} = [info_success.(params.areas{aidx}){early_late}; infQuant];

                        infQuant = MIsuccessFail.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail(chan,:);
                        infQuantSh = MIsuccessFail.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoShFail(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_fail.(params.areas{aidx}){early_late} = [info_fail.(params.areas{aidx}){early_late}; infQuant];

                    end
                end
            end
        end
       
        timeBins4LFP = round(newBinTimes/10)+150;

        m1_early_info = [];
        m1_early_info_all = [];
        for unit = 1:size(info_orig.M1{1},1)
            tmp_info = info_orig.M1{1}(unit,:);
            m1_early_info_all = [m1_early_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info];
            end
        end

        m1_early_info_success = [];
        m1_early_info_success_all = [];
        for unit = 1:size(info_success.M1{1},1)
            tmp_info = info_success.M1{1}(unit,:);
            m1_early_info_success_all = [m1_early_info_success_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info_success = [m1_early_info_success; tmp_info];
            end
        end

        m1_early_info_fail = [];
        m1_early_info_fail_all = [];
        for unit = 1:size(info_fail.M1{1},1)
            tmp_info = info_fail.M1{1}(unit,:);
            m1_early_info_fail_all = [m1_early_info_fail_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info_fail = [m1_early_info_fail; tmp_info];
            end
        end        

        m1_late_info = [];
        m1_late_info_all = [];
        for unit = 1:size(info_orig.M1{2},1)
            tmp_info = info_orig.M1{2}(unit,:);
            m1_late_info_all = [m1_late_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info];
            end
        end

        m1_late_info_success = [];
        m1_late_info_success_all = [];
        for unit = 1:size(info_success.M1{2},1)
            tmp_info = info_success.M1{2}(unit,:);
            m1_late_info_success_all = [m1_late_info_success_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info_success = [m1_late_info_success; tmp_info];
            end
        end

        m1_late_info_fail = [];
        m1_late_info_fail_all = [];
        for unit = 1:size(info_fail.M1{2},1)
            tmp_info = info_fail.M1{2}(unit,:);
            m1_late_info_fail_all = [m1_late_info_fail_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info_fail = [m1_late_info_fail; tmp_info];
            end
        end        

        dls_early_info = [];
        dls_early_info_all = [];
        for unit = 1:size(info_orig.DLS{1},1)
            tmp_info = info_orig.DLS{1}(unit,:);
            dls_early_info_all = [dls_early_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info];
            end
        end

        dls_early_info_success = [];
        dls_early_info_success_all = [];
        for unit = 1:size(info_success.DLS{1},1)
            tmp_info = info_success.DLS{1}(unit,:);
            dls_early_info_success_all = [dls_early_info_success_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info_success = [dls_early_info_success; tmp_info];
            end
        end

        dls_early_info_fail = [];
        dls_early_info_fail_all = [];
        for unit = 1:size(info_fail.DLS{1},1)
            tmp_info = info_fail.DLS{1}(unit,:);
            dls_early_info_fail_all = [dls_early_info_fail_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info_fail = [dls_early_info_fail; tmp_info];
            end
        end

        dls_late_info = [];
        dls_late_info_all = [];
        for unit = 1:size(info_orig.DLS{2},1)
            tmp_info = info_orig.DLS{2}(unit,:);
            dls_late_info_all = [dls_late_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info];
            end
        end

        dls_late_info_success = [];
        dls_late_info_success_all = [];
        for unit = 1:size(info_success.DLS{2},1)
            tmp_info = info_success.DLS{2}(unit,:);
            dls_late_info_success_all = [dls_late_info_success_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info_success = [dls_late_info_success; tmp_info];
            end
        end

        dls_late_info_fail = [];
        dls_late_info_fail_all = [];
        for unit = 1:size(info_fail.DLS{2},1)
            tmp_info = info_fail.DLS{2}(unit,:);
            dls_late_info_fail_all = [dls_late_info_fail_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info_fail = [dls_late_info_fail; tmp_info];
            end
        end

%%% ORGINAL
        
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
                ylim([0 0.05])
                xlim([12 86])
            subplot(3,2,6); hold on;
                plot(nanmean(dls_early_info)+nanstd(dls_early_info)/sqrt(size(dls_early_info,1)),'k')
                plot(nanmean(dls_early_info)-nanstd(dls_early_info)/sqrt(size(dls_early_info,1)),'k')   
                plot(nanmean(dls_late_info)+nanstd(dls_late_info)/sqrt(size(dls_late_info,1)),'r')
                plot(nanmean(dls_late_info)-nanstd(dls_late_info)/sqrt(size(dls_late_info,1)),'r')
                ylim([0 0.05])
                xlim([12 86])
                if save_params.save == 1
                    print([save_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_cascade_plot_mean_sem.eps'],'-painters','-depsc');
                end   

        figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(m1_early_info); abs(m1_late_info)];
                c = [ones(1,size(m1_early_info,1)) 2*ones(1,size(m1_late_info,1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 0.05]);           
                g(1,1).set_title('m1 lfp information');        
                y = [abs(dls_early_info); abs(dls_late_info)];
                c = [ones(1,size(dls_early_info,1)) 2*ones(1,size(dls_late_info,1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 0.05]);     
                g(1,2).set_title('dls lfp information');        
                g.draw();
                if save_params.save == 1
                    print([save_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_mean_sem.eps'],'-painters','-depsc');
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_M1_' params.reachFeatures{fidx} '_timingShift.eps'],'-painters','-depsc');
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_DLS_' params.reachFeatures{fidx} '_timingShift.eps'],'-painters','-depsc');
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_new_skilled_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
            end

%%% SUCCESS

    %%% PLOT INFO CASCADES AND MEAN / SEM

        figure('units','normalized','outerposition',[0 0 1 1])
            subplot(3,2,1); hold on;
                [~, i_early_info] = max(m1_early_info_success');
                [~, i_early_info] = sort(i_early_info);
                imagesc([m1_early_info_success(i_early_info,:)],[0 .15]);
                colorbar;
                plot([1 98],[length(i_early_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_early_info)])                
            subplot(3,2,3); hold on;
                [~, i_late_info] = max(m1_late_info_success');
                [~, i_late_info] = sort(i_late_info);
                imagesc([m1_late_info_success(i_late_info,:)],[0 .15]);
                plot([1 98],[length(i_late_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_late_info)])           
            subplot(3,2,2); hold on;
                [~, i_early_info] = max(dls_early_info_success');
                [~, i_early_info] = sort(i_early_info);
                imagesc([dls_early_info_success(i_early_info,:)],[0 .15]);
                plot([1 98],[length(i_early_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_early_info)])           
            subplot(3,2,4); hold on;
                [~, i_late_info] = max(dls_late_info_success');
                [~, i_late_info] = sort(i_late_info);
                imagesc([dls_late_info_success(i_late_info,:)],[0 .15]);
                plot([1 98],[length(i_late_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_late_info)])           
            subplot(3,2,5); hold on;
                plot(nanmean(m1_early_info_success)+nanstd(m1_early_info_success)/sqrt(size(m1_early_info_success,1)),'k')
                plot(nanmean(m1_early_info_success)-nanstd(m1_early_info_success)/sqrt(size(m1_early_info_success,1)),'k')   
                plot(nanmean(m1_late_info_success)+nanstd(m1_late_info_success)/sqrt(size(m1_late_info_success,1)),'r')
                plot(nanmean(m1_late_info_success)-nanstd(m1_late_info_success)/sqrt(size(m1_late_info_success,1)),'r')
                ylim([0 0.05])
                xlim([12 86])
            subplot(3,2,6); hold on;
                plot(nanmean(dls_early_info_success)+nanstd(dls_early_info_success)/sqrt(size(dls_early_info_success,1)),'k')
                plot(nanmean(dls_early_info_success)-nanstd(dls_early_info_success)/sqrt(size(dls_early_info_success,1)),'k')   
                plot(nanmean(dls_late_info_success)+nanstd(dls_late_info_success)/sqrt(size(dls_late_info_success,1)),'r')
                plot(nanmean(dls_late_info_success)-nanstd(dls_late_info_success)/sqrt(size(dls_late_info_success,1)),'r')
                ylim([0 0.05])
                xlim([12 86])
                if save_params.save == 1
                    print([save_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_cascade_plot_mean_sem.eps'],'-painters','-depsc');
                end   

        figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(m1_early_info_success); abs(m1_late_info_success)];
                c = [ones(1,size(m1_early_info_success,1)) 2*ones(1,size(m1_late_info_success,1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 0.05]);           
                g(1,1).set_title('m1 lfp information');        
                y = [abs(dls_early_info_success); abs(dls_late_info_success)];
                c = [ones(1,size(dls_early_info_success,1)) 2*ones(1,size(dls_late_info_success,1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 0.05]);     
                g(1,2).set_title('dls lfp information');        
                g.draw();
                if save_params.save == 1
                    print([save_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_mean_sem.eps'],'-painters','-depsc');
                end   

    %%% PEAK SHIFT INFO

        edge_m1_early = [];
        for unit = 1:size(m1_early_info_success,1)
            baseSD = 2*std(m1_early_info_success(unit,:));
            baseMean = mean(m1_early_info_success(unit,:));
            edge_m1_early = [edge_m1_early min(find(m1_early_info_success(unit,:)>baseMean+baseSD))];
        end
        edge_m1_early = edge_m1_early(edge_m1_early>33 & edge_m1_early<85);      

        edge_m1_late = [];
        for unit = 1:size(m1_late_info_success,1)
            baseSD = 2*std(m1_late_info_success(unit,:));
            baseMean = mean(m1_late_info_success(unit,:));
            edge_m1_late = [edge_m1_late min(find(m1_late_info_success(unit,:)>baseMean+baseSD))];
        end
        edge_m1_late = edge_m1_late(edge_m1_late>33 & edge_m1_late<85);     

        edge_dls_early = [];
        for unit = 1:size(dls_early_info_success,1)
            baseSD = 2*std(dls_early_info_success(unit,:));
            baseMean = mean(dls_early_info_success(unit,:));
            edge_dls_early = [edge_dls_early min(find(dls_early_info_success(unit,:)>baseMean+baseSD))];
        end
        edge_dls_early = edge_dls_early(edge_dls_early>33 & edge_dls_early<85);     

        edge_dls_late = [];
        for unit = 1:size(dls_late_info_success,1)
            baseSD = 2*std(dls_late_info_success(unit,:));
            baseMean = mean(dls_late_info_success(unit,:));
            edge_dls_late = [edge_dls_late min(find(dls_late_info_success(unit,:)>baseMean+baseSD))];
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_M1_' params.reachFeatures{fidx} '_timingShift.eps'],'-painters','-depsc');
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_DLS_' params.reachFeatures{fidx} '_timingShift.eps'],'-painters','-depsc');
            end

    %%% NEW NAIVE SKILLED COMPARISON HISTOGRAM/CDF

        figure;
            subplot(2,2,1); hold on;
                histogram(max(m1_early_info_success_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(m1_late_info_success_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(m1_early_info_success_all'),max(m1_late_info_success_all'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,3);
                hold on;
                cdfplot(max(m1_early_info_success_all'));
                cdfplot(max(m1_late_info_success_all'));
                xlim([0 0.25]);
                grid off;
            subplot(2,2,2); hold on;
                histogram(max(dls_early_info_success_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(dls_late_info_success_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(dls_early_info_success_all'),max(dls_late_info_success_all'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,4);
                hold on;
                cdfplot(max(dls_early_info_success_all'));
                cdfplot(max(dls_late_info_success_all'));
                xlim([0 0.25]);
                grid off;
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_new_skilled_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
            end


%%% FAIL

    %%% PLOT INFO CASCADES AND MEAN / SEM

        figure('units','normalized','outerposition',[0 0 1 1])
            subplot(3,2,1); hold on;
                [~, i_early_info] = max(m1_early_info_fail');
                [~, i_early_info] = sort(i_early_info);
                imagesc([m1_early_info_fail(i_early_info,:)],[0 .15]);
                colorbar;
                plot([1 98],[length(i_early_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_early_info)])                
            subplot(3,2,3); hold on;
                [~, i_late_info] = max(m1_late_info_fail');
                [~, i_late_info] = sort(i_late_info);
                imagesc([m1_late_info_fail(i_late_info,:)],[0 .15]);
                plot([1 98],[length(i_late_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_late_info)])           
            subplot(3,2,2); hold on;
                [~, i_early_info] = max(dls_early_info_fail');
                [~, i_early_info] = sort(i_early_info);
                imagesc([dls_early_info_fail(i_early_info,:)],[0 .15]);
                plot([1 98],[length(i_early_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_early_info)])           
            subplot(3,2,4); hold on;
                [~, i_late_info] = max(dls_late_info_fail');
                [~, i_late_info] = sort(i_late_info);
                imagesc([dls_late_info_fail(i_late_info,:)],[0 .15]);
                plot([1 98],[length(i_late_info)],'r','LineWidth',1)
                xlim([12 86])
                ylim([1 length(i_late_info)])           
            subplot(3,2,5); hold on;
                plot(nanmean(m1_early_info_fail)+nanstd(m1_early_info_fail)/sqrt(size(m1_early_info_fail,1)),'k')
                plot(nanmean(m1_early_info_fail)-nanstd(m1_early_info_fail)/sqrt(size(m1_early_info_fail,1)),'k')   
                plot(nanmean(m1_late_info_fail)+nanstd(m1_late_info_fail)/sqrt(size(m1_late_info_fail,1)),'r')
                plot(nanmean(m1_late_info_fail)-nanstd(m1_late_info_fail)/sqrt(size(m1_late_info_fail,1)),'r')
                ylim([0 0.05])
                xlim([12 86])
            subplot(3,2,6); hold on;
                plot(nanmean(dls_early_info_fail)+nanstd(dls_early_info_fail)/sqrt(size(dls_early_info_fail,1)),'k')
                plot(nanmean(dls_early_info_fail)-nanstd(dls_early_info_fail)/sqrt(size(dls_early_info_fail,1)),'k')   
                plot(nanmean(dls_late_info_fail)+nanstd(dls_late_info_fail)/sqrt(size(dls_late_info_fail,1)),'r')
                plot(nanmean(dls_late_info_fail)-nanstd(dls_late_info_fail)/sqrt(size(dls_late_info_fail,1)),'r')
                ylim([0 0.05])
                xlim([12 86])
                if save_params.save == 1
                    print([save_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_cascade_plot_mean_sem.eps'],'-painters','-depsc');
                end   

        figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(m1_early_info_fail); abs(m1_late_info_fail)];
                c = [ones(1,size(m1_early_info_fail,1)) 2*ones(1,size(m1_late_info_fail,1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 0.05]);           
                g(1,1).set_title('m1 lfp information');        
                y = [abs(dls_early_info_fail); abs(dls_late_info_fail)];
                c = [ones(1,size(dls_early_info_fail,1)) 2*ones(1,size(dls_late_info_fail,1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 0.05]);     
                g(1,2).set_title('dls lfp information');        
                g.draw();
                if save_params.save == 1
                    print([save_params.save_path '\LFP_information_' params.reachFeatures{fidx} '_mean_sem.eps'],'-painters','-depsc');
                end   

    %%% PEAK SHIFT INFO

        edge_m1_early = [];
        for unit = 1:size(m1_early_info_fail,1)
            baseSD = 2*std(m1_early_info_fail(unit,:));
            baseMean = mean(m1_early_info_fail(unit,:));
            edge_m1_early = [edge_m1_early min(find(m1_early_info_fail(unit,:)>baseMean+baseSD))];
        end
        edge_m1_early = edge_m1_early(edge_m1_early>33 & edge_m1_early<85);      

        edge_m1_late = [];
        for unit = 1:size(m1_late_info_fail,1)
            baseSD = 2*std(m1_late_info_fail(unit,:));
            baseMean = mean(m1_late_info_fail(unit,:));
            edge_m1_late = [edge_m1_late min(find(m1_late_info_fail(unit,:)>baseMean+baseSD))];
        end
        edge_m1_late = edge_m1_late(edge_m1_late>33 & edge_m1_late<85);     

        edge_dls_early = [];
        for unit = 1:size(dls_early_info_fail,1)
            baseSD = 2*std(dls_early_info_fail(unit,:));
            baseMean = mean(dls_early_info_fail(unit,:));
            edge_dls_early = [edge_dls_early min(find(dls_early_info_fail(unit,:)>baseMean+baseSD))];
        end
        edge_dls_early = edge_dls_early(edge_dls_early>33 & edge_dls_early<85);     

        edge_dls_late = [];
        for unit = 1:size(dls_late_info_fail,1)
            baseSD = 2*std(dls_late_info_fail(unit,:));
            baseMean = mean(dls_late_info_fail(unit,:));
            edge_dls_late = [edge_dls_late min(find(dls_late_info_fail(unit,:)>baseMean+baseSD))];
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_M1_' params.reachFeatures{fidx} '_timingShift.eps'],'-painters','-depsc');
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
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_DLS_' params.reachFeatures{fidx} '_timingShift.eps'],'-painters','-depsc');
            end

    %%% NEW NAIVE SKILLED COMPARISON HISTOGRAM/CDF

        figure;
            subplot(2,2,1); hold on;
                histogram(max(m1_early_info_fail_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(m1_late_info_fail_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(m1_early_info_fail_all'),max(m1_late_info_fail_all'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,3);
                hold on;
                cdfplot(max(m1_early_info_fail_all'));
                cdfplot(max(m1_late_info_fail_all'));
                xlim([0 0.25]);
                grid off;
            subplot(2,2,2); hold on;
                histogram(max(dls_early_info_fail_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(dls_late_info_fail_all'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(dls_early_info_fail_all'),max(dls_late_info_fail_all'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,4);
                hold on;
                cdfplot(max(dls_early_info_fail_all'));
                cdfplot(max(dls_late_info_fail_all'));
                xlim([0 0.25]);
                grid off;
            if save_params.save == 1
                print([save_params.save_path '\LFP_information_new_skilled_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
            end


% 
%     disp(['LFP | ' params.reachFeatures{fidx} ' | M1 | early: ' num2str(size(m1_early_info,1)) ' out of ' num2str(size(m1_early_info_all,1)) ' - ' num2str(100*size(m1_early_info,1)/size(m1_early_info_all,1))  ...
%                                                   '% | late: ' num2str(size(m1_late_info,1)) ' out of ' num2str(size(m1_late_info_all,1)) ' - ' num2str(100*size(m1_late_info,1)/size(m1_late_info_all,1))  '%'])
%     disp(['LFP | ' params.reachFeatures{fidx} ' | DLS | early: ' num2str(size(dls_early_info,1)) ' out of ' num2str(size(dls_early_info_all,1)) ' - ' num2str(100*size(dls_early_info,1)/size(dls_early_info_all,1))  ...
%                                                   '% | late: ' num2str(size(dls_late_info,1)) ' out of ' num2str(size(dls_late_info_all,1)) ' - ' num2str(100*size(dls_late_info,1)/size(dls_late_info_all,1))  '%'])
% 


    end
end
end

