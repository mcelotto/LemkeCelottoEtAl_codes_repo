function [] = plot_mutualInformation_LFP_comparison_newChan_allFeatures(MIout,LFP_early_late,params,clusterParams,save_params,newBinTimes)

allFeatureInfo_M1_naive = cell(2,length(params.reachFeatures)-1);
allFeatureInfo_M1_skilled = cell(2,length(params.reachFeatures)-1);
allFeatureInfo_DLS_naive = cell(2,length(params.reachFeatures)-1);
allFeatureInfo_DLS_skilled = cell(2,length(params.reachFeatures)-1);

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = 1:length(params.reachFeatures)-1 % not computed for success rate

        info_mean_all.M1 = cell(1,2);
        info_mean_all.DLS = cell(1,2);        
        info_mean_sig.M1 = cell(1,2);
        info_mean_sig.DLS = cell(1,2);
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        info_mean_sig.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
                    end
                end
            end
        end

    m1_early_info_all = [];
    m1_early_info = [];
    for unit = 1:size(info_mean_sig.M1{1},1)
        tmp_info = info_mean_sig.M1{1}(unit,:);
        tmp_info2 = info_mean_all.M1{1}(unit,:);
        [v,i] = max(tmp_info);        
        m1_early_info_all = [m1_early_info_all; tmp_info];
        if v>0.001 && i>=12 && i<=86
            m1_early_info = [m1_early_info; tmp_info];
        end
    end

    m1_late_info_all = [];
    m1_late_info = [];
    for unit = 1:size(info_mean_sig.M1{2},1)
        tmp_info = info_mean_sig.M1{2}(unit,:);
        tmp_info2 = info_mean_all.M1{2}(unit,:);
        [v,i] = max(tmp_info);
        m1_late_info_all = [m1_late_info_all; tmp_info];        
        if v>0.001 && i>=12 && i<=86
            m1_late_info = [m1_late_info; tmp_info];
        end
    end

    dls_early_info_all = [];
    dls_early_info = [];
    for unit = 1:size(info_mean_sig.DLS{1},1)
        tmp_info = info_mean_sig.DLS{1}(unit,:);
        tmp_info2 = info_mean_all.DLS{1}(unit,:);
        [v,i] = max(tmp_info);        
        dls_early_info_all = [dls_early_info_all; tmp_info];
        if v>0.001 && i>=12 && i<=86
            dls_early_info = [dls_early_info; tmp_info];
        end
    end

    dls_late_info_all = [];
    dls_late_info = [];
    for unit = 1:size(info_mean_sig.DLS{2},1)
        tmp_info = info_mean_sig.DLS{2}(unit,:);
        tmp_info2 = info_mean_all.DLS{2}(unit,:);
        [v,i] = max(tmp_info);        
        dls_late_info_all = [dls_late_info_all; tmp_info];
        if v>0.001 && i>=12 && i<=86
            dls_late_info = [dls_late_info; tmp_info];
        end
    end     
    
    allFeatureInfo_M1_naive{1,fidx} = m1_early_info_all;
    allFeatureInfo_M1_naive{2,fidx} = m1_early_info;
    allFeatureInfo_M1_skilled{1,fidx} = m1_late_info_all;
    allFeatureInfo_M1_skilled{2,fidx} = m1_late_info;
       
    allFeatureInfo_DLS_naive{1,fidx} = dls_early_info_all;
    allFeatureInfo_DLS_naive{2,fidx} = dls_early_info;
    allFeatureInfo_DLS_skilled{1,fidx} = dls_late_info_all;
    allFeatureInfo_DLS_skilled{2,fidx} = dls_late_info;
      
    disp(['LFP | ' params.reachFeatures{fidx} ' | M1 | early: ' num2str(size(m1_early_info,1)) ' out of ' num2str(size(m1_early_info_all,1)) ' - ' num2str(100*size(m1_early_info,1)/size(m1_early_info_all,1))  ...
                                                  '% | late: ' num2str(size(m1_late_info,1)) ' out of ' num2str(size(m1_late_info_all,1)) ' - ' num2str(100*size(m1_late_info,1)/size(m1_late_info_all,1))  '%'])
    disp(['LFP | ' params.reachFeatures{fidx} ' | DLS | early: ' num2str(size(dls_early_info,1)) ' out of ' num2str(size(dls_early_info_all,1)) ' - ' num2str(100*size(dls_early_info,1)/size(dls_early_info_all,1))  ...
                                                  '% | late: ' num2str(size(dls_late_info,1)) ' out of ' num2str(size(dls_late_info_all,1)) ' - ' num2str(100*size(dls_late_info,1)/size(dls_late_info_all,1))  '%'])
    end
end

figure;
    tiledlayout(2,1)
    nexttile;
        hold on;
        count = 1;
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            swarmchart(count*ones(1,size(allFeatureInfo_M1_naive{2,fidx},1)),max(allFeatureInfo_M1_naive{2,fidx}(:,12:86)'),10,[0 0 0],'filled')
            count = count + 3;
        end
        count = 2;
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            swarmchart(count*ones(1,size(allFeatureInfo_M1_skilled{2,fidx},1)),max(allFeatureInfo_M1_skilled{2,fidx}(:,12:86)'),10,[1 0 0],'filled')
            count = count + 3;
        end
        ylim([0 0.25])
    nexttile;
        hold on;
        count = 1;
        tmp_mean = [];
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            errorbar(count,mean(max(allFeatureInfo_M1_naive{1,fidx}(:,12:86)')),std(max(allFeatureInfo_M1_naive{1,fidx}(:,12:86)'))/sqrt(length(max(allFeatureInfo_M1_naive{1,fidx}(:,12:86)'))),'color',[0 0 0])
            tmp_mean = [tmp_mean mean(max(allFeatureInfo_M1_naive{1,fidx}(:,12:86)'))]
            count = count + 1;
        end
        plot(tmp_mean,'color',[0 0 0])
        count = 1;
        tmp_mean = [];        
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            errorbar(count,mean(max(allFeatureInfo_M1_skilled{1,fidx}(:,12:86)')),std(max(allFeatureInfo_M1_skilled{1,fidx}(:,12:86)'))/sqrt(length(max(allFeatureInfo_M1_skilled{1,fidx}(:,12:86)'))),'color',[1 0 0])
            tmp_mean = [tmp_mean mean(max(allFeatureInfo_M1_skilled{1,fidx}(:,12:86)'))]            
            count = count + 1;
        end
        plot(tmp_mean,'color',[1 0 0])
        ylim([0 0.07])
    sgtitle('M1 LFP')
    if save_params.save == 1
        print([save_params.save_path '\allFeature_Info_M1_LFP.eps'],'-painters','-depsc');
    end

figure;
    tiledlayout(2,1)
    nexttile;
        hold on;
        count = 1;
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            swarmchart(count*ones(1,size(allFeatureInfo_DLS_naive{2,fidx},1)),max(allFeatureInfo_DLS_naive{2,fidx}(:,12:86)'),10,[0 0 0],'filled')
            count = count + 3;
        end
        count = 2;
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            swarmchart(count*ones(1,size(allFeatureInfo_DLS_skilled{2,fidx},1)),max(allFeatureInfo_DLS_skilled{2,fidx}(:,12:86)'),10,[1 0 0],'filled')
            count = count + 3;
        end
        ylim([0 0.2])
    nexttile;
        hold on;
        count = 1;
        tmp_mean = [];
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            errorbar(count,mean(max(allFeatureInfo_DLS_naive{1,fidx}(:,12:86)')),std(max(allFeatureInfo_DLS_naive{1,fidx}(:,12:86)'))/sqrt(length(max(allFeatureInfo_DLS_naive{1,fidx}(:,12:86)'))),'color',[0 0 0])
            tmp_mean = [tmp_mean mean(max(allFeatureInfo_DLS_naive{1,fidx}(:,12:86)'))]
            count = count + 1;
        end
        plot(tmp_mean,'color',[0 0 0])
        count = 1;
        tmp_mean = [];        
        for fidx = [13 9 3 4 10 6 5 12 2 1] 
            errorbar(count,mean(max(allFeatureInfo_DLS_skilled{1,fidx}(:,12:86)')),std(max(allFeatureInfo_DLS_skilled{1,fidx}(:,12:86)'))/sqrt(length(max(allFeatureInfo_DLS_skilled{1,fidx}(:,12:86)'))),'color',[1 0 0])
            tmp_mean = [tmp_mean mean(max(allFeatureInfo_DLS_skilled{1,fidx}(:,12:86)'))]            
            count = count + 1;
        end
        plot(tmp_mean,'color',[1 0 0])
        ylim([0 0.05])
    sgtitle('DLS LFP')
    if save_params.save == 1
        print([save_params.save_path '\allFeature_Info_DLS_LFP.eps'],'-painters','-depsc');
    end

%% NAIVE VS SKILLED

early_m1 = [];
late_m1 = [];
early_dls = [];
late_dls = [];
for fidx = 1:13
    early_m1 = [early_m1; max(allFeatureInfo_M1_naive{1,fidx}')];
    late_m1 = [late_m1; max(allFeatureInfo_M1_skilled{1,fidx}')];
    early_dls = [early_dls; max(allFeatureInfo_DLS_naive{1,fidx}')];
    late_dls = [late_dls; max(allFeatureInfo_DLS_skilled{1,fidx}')];
end

figure;
    subplot(1,2,1); hold on;
        cdfplot(mean(early_m1));
        cdfplot(mean(late_m1));
        xlim([0 0.2]);
        grid off;
        [p h] = ranksum(mean(early_m1),mean(late_m1));
        title(['M1 | num early chans ' num2str(length(mean(early_m1))) ' and late chans ' num2str(length(mean(late_m1))) ' |  p = ' num2str(p)])
    subplot(1,2,2); hold on;
        cdfplot(mean(early_dls));
        cdfplot(mean(late_dls));
        xlim([0 0.2]);
        grid off;
        [p h] = ranksum(mean(early_dls),mean(late_dls));
        title(['DLS | num early chans ' num2str(length(mean(early_dls))) ' and late chans ' num2str(length(mean(late_dls))) ' |  p = ' num2str(p)])
    if save_params.save == 1
        print([save_params.save_path '\lfp_info_niave_skilled_spikes.eps'],'-painters','-depsc');
    end   


%     figure;
% 
%         %%% M1 MEAN INFO
%         subplot(2,2,1); hold on;
%         
%             all_m1_early_mean = [];
%             all_m1_early_sem = [];
%             for n = 1:14
%                 all_m1_early_mean = [all_m1_early_mean mean(max(allFeatureInfo_M1_naive{1,n}'))];
%                 all_m1_early_sem = [all_m1_early_sem std(max(allFeatureInfo_M1_naive{1,n}'))/sqrt(length(max(allFeatureInfo_M1_naive{1,n}')))];
%             end
%             all_m1_late_mean = [];
%             all_m1_late_sem = [];
%             for n = 1:14
%                 all_m1_late_mean = [all_m1_late_mean mean(max(allFeatureInfo_M1_skilled{1,n}'))];
%                 all_m1_late_sem = [all_m1_late_sem std(max(allFeatureInfo_M1_skilled{1,n}'))/sqrt(length(max(allFeatureInfo_M1_skilled{1,n}')))];
%             end
%     
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 errorbar(count,all_m1_early_mean(n),all_m1_early_sem(n),'color','k')
%                 count = count + 1;
%             end
%             plot(all_m1_early_mean([13 9 3 4 10 6 5 12 2 1]),'k')
%             
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 errorbar(count,all_m1_late_mean(n),all_m1_late_sem(n),'color','r')
%                 count = count + 1;
%             end
%             plot(all_m1_late_mean([13 9 3 4 10 6 5 12 2 1]),'r')
%             ylim([0 0.07])
%             title('M1 MEAN INFO')
% 
%         %%% DLS MEAN INFO
%         subplot(2,2,2); hold on;
%         
%             all_dls_early_mean = [];
%             all_dls_early_sem = [];
%             for n = 1:14
%                 all_dls_early_mean = [all_dls_early_mean mean(max(allFeatureInfo_DLS_naive{1,n}'))];
%                 all_dls_early_sem = [all_dls_early_sem std(max(allFeatureInfo_DLS_naive{1,n}'))/sqrt(length(max(allFeatureInfo_DLS_naive{1,n}')))];
%             end
%             all_dls_late_mean = [];
%             all_dls_late_sem = [];
%             for n = 1:14
%                 all_dls_late_mean = [all_dls_late_mean mean(max(allFeatureInfo_DLS_skilled{1,n}'))];
%                 all_dls_late_sem = [all_dls_late_sem std(max(allFeatureInfo_DLS_skilled{1,n}'))/sqrt(length(max(allFeatureInfo_DLS_skilled{1,n}')))];
%             end
%     
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 errorbar(count,all_dls_early_mean(n),all_dls_early_sem(n),'color','k')
%                 count = count + 1;
%             end
%             plot(all_dls_early_mean([13 9 3 4 10 6 5 12 2 1]),'k')
%             
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 errorbar(count,all_dls_late_mean(n),all_dls_late_sem(n),'color','r')
%                 count = count + 1;
%             end
%             plot(all_dls_late_mean([13 9 3 4 10 6 5 12 2 1]),'r')
%             ylim([0 0.07])            
%             title('DLS MEAN INFO')
% 
%         %%% M1 SIG INFO CHANNELS
%         subplot(2,2,3); hold on;
% 
%             all_m1_early_sigPerc = [];
%             for n = 1:14
%                 all_m1_early_sigPerc = [all_m1_early_sigPerc size(allFeatureInfo_M1_naive{2,n},1)/size(allFeatureInfo_M1_naive{1,n},1)];
%             end
%             all_m1_late_sigPerc = [];
%             for n = 1:14
%                 all_m1_late_sigPerc = [all_m1_late_sigPerc size(allFeatureInfo_M1_skilled{2,n},1)/size(allFeatureInfo_M1_skilled{1,n},1)];
%             end
%     
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 scatter(count,all_m1_early_sigPerc(n),30,[0 0 0])
%                 count = count + 1;
%             end
%             plot(all_m1_early_sigPerc([13 9 3 4 10 6 5 12 2 1]),'k')
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 scatter(count,all_m1_late_sigPerc(n),30,[1 0 0])
%                 count = count + 1;
%             end
%             plot(all_m1_late_sigPerc([13 9 3 4 10 6 5 12 2 1]),'r')
%             ylim([0 .7])
%             title('M1 SIG CHAN %')
% 
%         %%% DLS SIG INFO CHANNELS
%         subplot(2,2,4); hold on;
% 
%             all_dls_early_sigPerc = [];
%             for n = 1:14
%                 all_dls_early_sigPerc = [all_dls_early_sigPerc size(allFeatureInfo_DLS_naive{2,n},1)/size(allFeatureInfo_DLS_naive{1,n},1)];
%             end
%             all_dls_late_sigPerc = [];
%             for n = 1:14
%                 all_dls_late_sigPerc = [all_dls_late_sigPerc size(allFeatureInfo_DLS_skilled{2,n},1)/size(allFeatureInfo_DLS_skilled{1,n},1)];
%             end
% 
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 scatter(count,all_dls_early_sigPerc(n),30,[0 0 0])
%                 count = count + 1;
%             end
%             plot(all_dls_early_sigPerc([13 9 3 4 10 6 5 12 2 1]),'k')
%             count = 1;
%             for n = [13 9 3 4 10 6 5 12 2 1] 
%                 scatter(count,all_dls_late_sigPerc(n),30,[1 0 0])
%                 count = count + 1;
%             end
%             plot(all_dls_late_sigPerc([13 9 3 4 10 6 5 12 2 1]),'r')
%             ylim([0 .7])            
%             title('DLS SIG CHAN %')


end

