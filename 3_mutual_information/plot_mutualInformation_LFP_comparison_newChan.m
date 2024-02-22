function [] = plot_mutualInformation_LFP_comparison_newChan(MIout,LFP_early_late,params,clusterParams,f2plot,plot_params,newBinTimes)
%% ACTIVITY CASCADES
    
    timeBins4LFP = round(newBinTimes/10)+150;
    n_lfp = 1;
    
        raw_LFP.M1 = cell(1,2);
        raw_LFP.M1{1}=cell(1,1);
        raw_LFP.M1{2}=cell(1,1);
        raw_LFP.DLS = cell(1,2);
        raw_LFP.DLS{1}=cell(1,1);
        raw_LFP.DLS{2}=cell(1,1);
    
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{1}).info,1)
                        raw_LFP.(params.areas{aidx}){early_late}{end+1} =  zscore(squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:)))';
    %                     figure; imagesc(zscore(squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:)))')
                    end
                end
            end
        end
    
        m1_early_LFP = [];
        for chan = 1:length(raw_LFP.M1{1})-1
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.M1{1}{chan+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            m1_early_LFP = [m1_early_LFP; zscore(tmp_LFP)];
        end
    
        m1_late_LFP = [];
        for chan = 1:length(raw_LFP.M1{2})-1
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.M1{2}{chan+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            m1_late_LFP = [m1_late_LFP; zscore(tmp_LFP)];
        end
    
        dls_early_LFP = [];
        for chan = 1:length(raw_LFP.DLS{1})-1
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.DLS{1}{chan+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            dls_early_LFP = [dls_early_LFP; zscore(tmp_LFP)];
        end
    
        dls_late_LFP = [];
        for chan = 1:length(raw_LFP.DLS{2})-1
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.DLS{2}{chan+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            dls_late_LFP = [dls_late_LFP; zscore(tmp_LFP)];
        end
    
    %    figure('units','normalized','outerposition',[0 0 1 1])
            figure;
                subplot(4,2,1); hold on;
                    [vmax, imax] = max(m1_early_LFP');
                    [vmin, imin] = min(m1_early_LFP');
                    tmp_neg = m1_early_LFP(abs(vmin)>=abs(vmax),:);
                    tmp_pos = m1_early_LFP(abs(vmin)<abs(vmax),:);
                    imin = imin(abs(vmin)>=abs(vmax));
                    imax = imax(abs(vmin)<abs(vmax));
                    [~, imax] = sort(imax);
                    [~, imin] = sort(imin);
                    imagesc([tmp_neg(imin,:); tmp_pos(imax,:)],[-5 5]);
                    plot([61.5 61.5],[1 length([imin imax])],'k','LineWidth',1,'LineStyle','--')
                    xlim([12 86])
                    ylim([1 length([imin imax])])
                    title(num2str(length([imin imax])))
                subplot(4,2,2); hold on;
                    [vmax, imax] = max(m1_late_LFP');
                    [vmin, imin] = min(m1_late_LFP');
                    tmp_neg = m1_late_LFP(abs(vmin)>=abs(vmax),:);
                    tmp_pos = m1_late_LFP(abs(vmin)<abs(vmax),:);
                    imin = imin(abs(vmin)>=abs(vmax));
                    imax = imax(abs(vmin)<abs(vmax));
                    [~, imax] = sort(imax);
                    [~, imin] = sort(imin);
                    imagesc([tmp_neg(imin,:); tmp_pos(imax,:)],[-5 5]);
                    plot([61.5 61.5],[1 length([imin imax])],'k','LineWidth',1,'LineStyle','--')
                    xlim([12 86])
                    ylim([1 length([imin imax])])
                    title(num2str(length([imin imax])))                    
                subplot(4,2,3); hold on;
                    [vmax, imax] = max(dls_early_LFP');
                    [vmin, imin] = min(dls_early_LFP');
                    tmp_neg = dls_early_LFP(abs(vmin)>=abs(vmax),:);
                    tmp_pos = dls_early_LFP(abs(vmin)<abs(vmax),:);
                    imin = imin(abs(vmin)>=abs(vmax));
                    imax = imax(abs(vmin)<abs(vmax));
                    [~, imax] = sort(imax);
                    [~, imin] = sort(imin);
                    imagesc([tmp_neg(imin,:); tmp_pos(imax,:)],[-5 5]);
                    plot([61.5 61.5],[1 length([imin imax])],'k','LineWidth',1,'LineStyle','--')
                    xlim([12 86])
                    ylim([1 length([imin imax])])
                    title(num2str(length([imin imax])))                    
                subplot(4,2,4); hold on;
                    [vmax, imax] = max(dls_late_LFP');
                    [vmin, imin] = min(dls_late_LFP');
                    tmp_neg = dls_late_LFP(abs(vmin)>=abs(vmax),:);
                    tmp_pos = dls_late_LFP(abs(vmin)<abs(vmax),:);
                    imin = imin(abs(vmin)>=abs(vmax));
                    imax = imax(abs(vmin)<abs(vmax));
                    [~, imax] = sort(imax);
                    [~, imin] = sort(imin);
                    imagesc([tmp_neg(imin,:); tmp_pos(imax,:)],[-5 5]);
                    plot([61.5 61.5],[1 length([imin imax])],'k','LineWidth',1,'LineStyle','--')
                    xlim([12 86])
                    ylim([1 length([imin imax])])
                    title(num2str(length([imin imax])))                    
                subplot(4,2,5); hold on;
                    plot(nanmean(abs(m1_early_LFP))+nanstd(abs(m1_early_LFP))/sqrt(size(abs(m1_early_LFP),1)),'k')
                    plot(nanmean(abs(m1_early_LFP))-nanstd(abs(m1_early_LFP))/sqrt(size(abs(m1_early_LFP),1)),'k')
                    ylim([0 2])
                    xlim([12 86])
                subplot(4,2,6); hold on;
                    plot(nanmean(abs(m1_late_LFP))+nanstd(abs(m1_late_LFP))/sqrt(size(abs(m1_late_LFP),1)),'k')
                    plot(nanmean(abs(m1_late_LFP))-nanstd(abs(m1_late_LFP))/sqrt(size(abs(m1_late_LFP),1)),'k')
                     ylim([0 2])
                    xlim([12 86])
                subplot(4,2,7); hold on;
                    plot(nanmean(abs(dls_early_LFP))+nanstd(abs(dls_early_LFP))/sqrt(size(abs(dls_early_LFP),1)),'r')
                    plot(nanmean(abs(dls_early_LFP))-nanstd(abs(dls_early_LFP))/sqrt(size(abs(dls_early_LFP),1)),'r')
                    ylim([0 2])
                    xlim([12 86])
                subplot(4,2,8); hold on;
                    plot(nanmean(abs(dls_late_LFP))+nanstd(abs(dls_late_LFP))/sqrt(size(abs(dls_late_LFP),1)),'r')
                    plot(nanmean(abs(dls_late_LFP))-nanstd(abs(dls_late_LFP))/sqrt(size(abs(dls_late_LFP),1)),'r')
                     ylim([0 2])
                    xlim([12 86])

                if plot_params.save == 1
                    print([plot_params.save_path '\lfp_activity_cascade_plots.eps'],'-painters','-depsc');
                end   
 
%% ACTIVITY TIMING

%     baseline_period = 1:18;
%     sd_thresh = 4;
%     onset_prc = 0.66;

    baseline_period = 1:36;
    sd_thresh = 4;
    onset_prc = 0.66;
    
    m1_early_onset = [];
    m1_early_peak = [];
    for unit = 1:size(m1_early_LFP,1)
        baseSD = std(m1_early_LFP(unit,baseline_period));
        baseMean = mean(m1_early_LFP(unit,baseline_period));
        [v, i] = max(m1_early_LFP(unit,12:86)); % -1s to +500ms
        if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
%             disp([num2str(v) ' | ' num2str(11+i)])
            m1_early_onset = [m1_early_onset min(find(m1_early_LFP(unit,:)>v*onset_prc))];
            m1_early_peak = [m1_early_peak 11+i];
        end
    end
    
    m1_late_onset = [];
    m1_late_peak = [];
    for unit = 1:size(m1_late_LFP,1)
        baseSD = std(m1_late_LFP(unit,baseline_period));
        baseMean = mean(m1_late_LFP(unit,baseline_period));
        [v, i] = max(m1_late_LFP(unit,12:86)); % -1s to +500ms
        if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
            m1_late_onset = [m1_late_onset min(find(m1_late_LFP(unit,:)>v*onset_prc))];
            m1_late_peak = [m1_late_peak 11+i];
        end
    end
    
    dls_early_onset = [];
    dls_early_peak = [];
    for unit = 1:size(dls_early_LFP,1)
        baseSD = std(dls_early_LFP(unit,baseline_period));
        baseMean = mean(dls_early_LFP(unit,baseline_period));
        [v, i] = max(dls_early_LFP(unit,12:86)); % -1s to +500ms
        if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
            dls_early_onset = [dls_early_onset min(find(dls_early_LFP(unit,:)>v*onset_prc))];
            dls_early_peak = [dls_early_peak 11+i];
        end
    end
    
    dls_late_onset = [];
    dls_late_peak = [];
    for unit = 1:size(dls_late_LFP,1)
        baseSD = std(dls_late_LFP(unit,baseline_period));
        baseMean = mean(dls_late_LFP(unit,baseline_period));
        [v, i] = max(dls_late_LFP(unit,12:86)); % -1s to +500ms
        if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
            dls_late_onset = [dls_late_onset min(find(dls_late_LFP(unit,:)>v*onset_prc))];
            dls_late_peak = [dls_late_peak 11+i];
        end
    end

    figure;

        subplot(2,2,1); hold on;
            im1 = m1_early_onset;
            idls = dls_early_onset;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            xlim([1 38])
            ylim([0 0.15])
            [h, p, ci, stats] = ttest2(im1,idls);
            [p2, ~, ~] = ranksum(im1,idls);
            title(['onset | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        subplot(2,2,2); hold on;
            im1 = m1_late_onset;
            idls = dls_late_onset;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            xlim([1 38])
            ylim([0 0.15])
            [h, p, ci, stats] = ttest2(im1,idls);
            [p2, ~, ~] = ranksum(im1,idls);
            title(['onset | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        subplot(2,2,3); hold on;
            im1 = m1_early_peak;
            idls = dls_early_peak;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            xlim([1 38])
            ylim([0 0.15])
            [h, p, ci, stats] = ttest2(im1,idls);
            [p2, ~, ~] = ranksum(im1,idls);      
            title(['peak | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        subplot(2,2,4); hold on;
            im1 = m1_late_peak;
            idls = dls_late_peak;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            xlim([1 38])
            ylim([0 0.15])
            [h, p, ci, stats] = ttest2(im1,idls);
            [p2, ~, ~] = ranksum(im1,idls);
            title(['peak | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        if plot_params.save == 1
            print([plot_params.save_path '\lfp_activity_timing.eps'],'-painters','-depsc');
        end

%% INFO CASCADES MAX VEL/TRAJ DIST
                
    info_mean_sig.M1 = cell(2,2);
    info_mean_sig.DLS = cell(2,2);
    info_mean_all.M1 = cell(2,2);
    info_mean_all.DLS = cell(2,2);    
    for animal = 1:numel(params.animals)
        for aidx = 1:length(params.areas)
            for early_late = 1:2
                for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{1}).info,1)
                    feature = 1;
                    for fidx = [9 13]
                        infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){early_late,feature} = [info_mean_all.(params.areas{aidx}){early_late,feature}; infQuant];                        
                        infQuant(~tmp_sig) = 0;                        
                        info_mean_sig.(params.areas{aidx}){early_late,feature} = [info_mean_sig.(params.areas{aidx}){early_late,feature}; infQuant];
                        feature = feature+1;
                    end
                end
            end
        end
    end

    m1_early_info = [];
    for unit = 1:size(info_mean_sig.M1{1,1},1)
        tmp_info = [info_mean_sig.M1{1,1}(unit,:); info_mean_sig.M1{1,2}(unit,:)];
        tmp_info2 = [info_mean_all.M1{1,1}(unit,:); info_mean_all.M1{1,2}(unit,:)];
        [v,i] = max(mean(tmp_info));
        if v>0.001 && i>=12 && i<=86
            m1_early_info = [m1_early_info; mean(tmp_info2)];
        end
    end

    m1_late_info = [];
    for unit = 1:size(info_mean_sig.M1{2,1},1)
        tmp_info = [info_mean_sig.M1{2,1}(unit,:); info_mean_sig.M1{2,2}(unit,:)];
        tmp_info2 = [info_mean_all.M1{2,1}(unit,:); info_mean_all.M1{2,2}(unit,:)];
        [v,i] = max(mean(tmp_info));
        if v>0.001 && i>=12 && i<=86
            m1_late_info = [m1_late_info; mean(tmp_info2)];
        end
    end

    dls_early_info = [];
    for unit = 1:size(info_mean_sig.DLS{1,1},1)
        tmp_info = [info_mean_sig.DLS{1,1}(unit,:); info_mean_sig.DLS{1,2}(unit,:)];
        tmp_info2 = [info_mean_all.DLS{1,1}(unit,:); info_mean_all.DLS{1,2}(unit,:)];
        [v,i] = max(mean(tmp_info));
        if v>0.001 && i>=12 && i<=86
            dls_early_info = [dls_early_info; mean(tmp_info2)];
        end
    end

    dls_late_info = [];
    for unit = 1:size(info_mean_sig.DLS{2,1},1)
        tmp_info = [info_mean_sig.DLS{2,1}(unit,:); info_mean_sig.DLS{2,2}(unit,:)];
        tmp_info2 = [info_mean_all.DLS{2,1}(unit,:); info_mean_all.DLS{2,2}(unit,:)];
        [v,i] = max(mean(tmp_info));
        if v>0.001 && i>=12 && i<=86
            dls_late_info = [dls_late_info; mean(tmp_info2)];                        
        end
    end  

    figure;
        subplot(4,2,1); hold on;
            [~, i_early_info] = max(m1_early_info');
            [~, i_early_info] = sort(i_early_info);
            imagesc([m1_early_info(i_early_info,:)],[0 .15]);
            plot([1 98],[length(i_early_info)],'r','LineWidth',1)
            xlim([12 86])
            ylim([1 length(i_early_info)])       
            title(num2str(length(i_early_info)))            
        subplot(4,2,2); hold on;
            [~, i_late_info] = max(m1_late_info');
            [~, i_late_info] = sort(i_late_info);
            imagesc([m1_late_info(i_late_info,:)],[0 .15]);
            plot([1 98],[length(i_late_info)],'r','LineWidth',1)
            xlim([12 86])
            ylim([1 length(i_late_info)])       
            title(num2str(length(i_late_info)))                        
        subplot(4,2,3); hold on;
            [~, i_early_info] = max(dls_early_info');
            [~, i_early_info] = sort(i_early_info);
            imagesc([dls_early_info(i_early_info,:)],[0 .15]);
            plot([1 98],[length(i_early_info)],'r','LineWidth',1)
            xlim([12 86])
            ylim([1 length(i_early_info)]) 
            title(num2str(length(i_early_info)))                                    
        subplot(4,2,4); hold on;
            [~, i_late_info] = max(dls_late_info');
            [~, i_late_info] = sort(i_late_info);
            imagesc([dls_late_info(i_late_info,:)],[0 .15]);
            plot([1 98],[length(i_late_info)],'r','LineWidth',1)
            xlim([12 86])
            ylim([1 length(i_late_info)])    
            title(num2str(length(i_late_info)))                                    
        subplot(4,2,5); hold on;
            plot(nanmean(m1_early_info)+nanstd(m1_early_info)/sqrt(size(m1_early_info,1)),'k')
            plot(nanmean(m1_early_info)-nanstd(m1_early_info)/sqrt(size(m1_early_info,1)),'k')  
            ylim([0 0.04])
            xlim([12 86])
        subplot(4,2,6); hold on;
            plot(nanmean(m1_late_info)+nanstd(m1_late_info)/sqrt(size(m1_late_info,1)),'k')
            plot(nanmean(m1_late_info)-nanstd(m1_late_info)/sqrt(size(m1_late_info,1)),'k')
            ylim([0 0.04])
            xlim([12 86])
        subplot(4,2,7); hold on;
            plot(nanmean(dls_early_info)+nanstd(dls_early_info)/sqrt(size(dls_early_info,1)),'r')
            plot(nanmean(dls_early_info)-nanstd(dls_early_info)/sqrt(size(dls_early_info,1)),'r')   
            ylim([0 0.04])
            xlim([12 86])
        subplot(4,2,8); hold on;
            plot(nanmean(dls_late_info)+nanstd(dls_late_info)/sqrt(size(dls_late_info,1)),'r')
            plot(nanmean(dls_late_info)-nanstd(dls_late_info)/sqrt(size(dls_late_info,1)),'r')
            ylim([0 0.04])
            xlim([12 86])

        if plot_params.save == 1
            print([plot_params.save_path '\lfp_information_cascade_plots.eps'],'-painters','-depsc');
        end   

        
%         save('G:\CurrBiol_revision\figures\figure_2\movOnset_comparison\PT_LFP_info.mat','m1_early_info','m1_late_info','dls_early_info','dls_late_info')

%% INFO TIMING MAX VEL/TRAJ DIST

%     baseline_period = 1:18;
%     sd_thresh = 4;
%     onset_prc = 0.66;

    baseline_period = 1:36;
    sd_thresh = 4;
    onset_prc = 0.66;

    %%% INFO

        m1_early_info = [];
        for unit = 1:size(info_mean_sig.M1{1,1},1)
            tmp_info = [info_mean_sig.M1{1,1}(unit,:)];
            tmp_info2 = [info_mean_all.M1{1,1}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info2];
            end
        end
        m1_early_info2 = [];
        for unit = 1:size(info_mean_sig.M1{1,2},1)
            tmp_info = [info_mean_sig.M1{1,2}(unit,:)];
            tmp_info2 = [info_mean_all.M1{1,2}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                m1_early_info2 = [m1_early_info2; tmp_info2];
            end
        end
    
        dls_early_info = [];
        for unit = 1:size(info_mean_sig.DLS{1,1},1)
            tmp_info = [info_mean_sig.DLS{1,1}(unit,:)];
            tmp_info2 = [info_mean_all.DLS{1,1}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info2];
            end
        end
        dls_early_info2 = [];
        for unit = 1:size(info_mean_sig.DLS{1,2},1)
            tmp_info = [info_mean_sig.DLS{1,2}(unit,:)];
            tmp_info2 = [info_mean_all.DLS{1,2}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                dls_early_info2 = [dls_early_info2; tmp_info2];
            end
        end
        m1_late_info = [];
        for unit = 1:size(info_mean_sig.M1{2,1},1)
            tmp_info = [info_mean_sig.M1{2,1}(unit,:)];
            tmp_info2 = [info_mean_all.M1{2,1}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info2];
            end
        end
        m1_late_info2 = [];
        for unit = 1:size(info_mean_sig.M1{2,2},1)
            tmp_info = [info_mean_sig.M1{2,2}(unit,:)];
            tmp_info2 = [info_mean_all.M1{2,2}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                m1_late_info2 = [m1_late_info2; tmp_info2];
            end
        end
    
        dls_late_info = [];
        for unit = 1:size(info_mean_sig.DLS{2,1},1)
            tmp_info = [info_mean_sig.DLS{2,1}(unit,:)];
            tmp_info2 = [info_mean_all.DLS{2,1}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info2];
            end
        end
        dls_late_info2 = [];
        for unit = 1:size(info_mean_sig.DLS{2,2},1)
            tmp_info = [info_mean_sig.DLS{2,2}(unit,:)];
            tmp_info2 = [info_mean_all.DLS{2,2}(unit,:)];
            [v,i] = max(tmp_info);
            if v>0.001 && i>=12 && i<=86
                dls_late_info2 = [dls_late_info2; tmp_info2];
            end
        end      

    %%% FEATURE 1
    
        m1_early_onset = [];
        m1_early_peak = [];
        for unit = 1:size(m1_early_info,1)
            baseSD = nanstd(m1_early_info(unit,baseline_period));
            baseMean = nanmean(m1_early_info(unit,baseline_period));
            [v, i] = max(m1_early_info(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                m1_early_onset = [m1_early_onset min(find(m1_early_info(unit,:)>v*onset_prc))];
                m1_early_peak = [m1_early_peak 11+i];
            end
        end
        
        m1_late_onset = [];
        m1_late_peak = [];
        for unit = 1:size(m1_late_info,1)
            baseSD = nanstd(m1_late_info(unit,baseline_period));
            baseMean = nanmean(m1_late_info(unit,baseline_period));
            [v, i] = max(m1_late_info(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                m1_late_onset = [m1_late_onset min(find(m1_late_info(unit,:)>v*onset_prc))];
                m1_late_peak = [m1_late_peak 11+i];
            end
        end
        
        dls_early_onset = [];
        dls_early_peak = [];
        for unit = 1:size(dls_early_info,1)
            baseSD = nanstd(dls_early_info(unit,baseline_period));
            baseMean = nanmean(dls_early_info(unit,baseline_period));
            [v, i] = max(dls_early_info(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                dls_early_onset = [dls_early_onset min(find(dls_early_info(unit,:)>v*onset_prc))];
                dls_early_peak = [dls_early_peak 11+i];
            end
        end
        
        dls_late_onset = [];
        dls_late_peak = [];
        for unit = 1:size(dls_late_info,1)
            baseSD = nanstd(dls_late_info(unit,baseline_period));
            baseMean = nanmean(dls_late_info(unit,baseline_period));
            [v, i] = max(dls_late_info(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                dls_late_onset = [dls_late_onset min(find(dls_late_info(unit,:)>v*onset_prc))];
                dls_late_peak = [dls_late_peak 11+i];
            end
        end

    %%% FEATURE 2
    
        m1_early_onset2 = [];
        m1_early_peak2 = [];
        for unit = 1:size(m1_early_info2,1)
            baseSD = nanstd(m1_early_info2(unit,baseline_period));
            baseMean = nanmean(m1_early_info2(unit,baseline_period));
            [v, i] = max(m1_early_info2(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                m1_early_onset2 = [m1_early_onset2 min(find(m1_early_info2(unit,:)>v*onset_prc))];
                m1_early_peak2 = [m1_early_peak2 11+i];
            end
        end
        
        m1_late_onset2 = [];
        m1_late_peak2 = [];
        for unit = 1:size(m1_late_info2,1)
            baseSD = nanstd(m1_late_info2(unit,baseline_period));
            baseMean = nanmean(m1_late_info2(unit,baseline_period));
            [v, i] = max(m1_late_info2(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                m1_late_onset2 = [m1_late_onset2 min(find(m1_late_info2(unit,:)>v*onset_prc))];
                m1_late_peak2 = [m1_late_peak2 11+i];
            end
        end
        
        dls_early_onset2 = [];
        dls_early_peak2 = [];
        for unit = 1:size(dls_early_info2,1)
            baseSD = nanstd(dls_early_info2(unit,baseline_period));
            baseMean = nanmean(dls_early_info2(unit,baseline_period));
            [v, i] = max(dls_early_info2(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                dls_early_onset2 = [dls_early_onset2 min(find(dls_early_info2(unit,:)>v*onset_prc))];
                dls_early_peak2 = [dls_early_peak2 11+i];
            end
        end
        
        dls_late_onset2 = [];
        dls_late_peak2 = [];
        for unit = 1:size(dls_late_info2,1)
            baseSD = nanstd(dls_late_info2(unit,baseline_period));
            baseMean = nanmean(dls_late_info2(unit,baseline_period));
            [v, i] = max(dls_late_info2(unit,12:86)); % -1s to +500ms
            if v>baseMean+(baseSD*sd_thresh) & i~=length(12:86)
                dls_late_onset2 = [dls_late_onset2 min(find(dls_late_info2(unit,:)>v*onset_prc))];
                dls_late_peak2 = [dls_late_peak2 11+i];
            end
        end

    figure;

        subplot(2,2,1); hold on;
            im1 = m1_early_onset;
            idls = dls_early_onset;
            im12 = m1_early_onset2;
            idls2 = dls_early_onset2;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts(im12,[12:2:86])/sum(histcounts(im12,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls2,[12:2:86])/sum(histcounts(idls2,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts([im1 im12],[12:2:86])/sum(histcounts([im1 im12],[12:2:86])),'gaussian',5),'color',[0 0 0],'LineWidth',2)
            plot(smoothdata(histcounts([idls idls2],[12:2:86])/sum(histcounts([idls idls2],[12:2:86])),'gaussian',5),'color',[1 0 0],'LineWidth',2)
            xlim([1 38]); ylim([0 0.2]);
            [h, p, ci, stats] = ttest2([im1 im12],[idls idls2]);
            [p2, ~, ~] = ranksum([im1 im12],[idls idls2]);
            title(['onset | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        subplot(2,2,2); hold on;
            im1 = m1_late_onset;
            idls = dls_late_onset;
            im12 = m1_late_onset2;
            idls2 = dls_late_onset2;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts(im12,[12:2:86])/sum(histcounts(im12,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls2,[12:2:86])/sum(histcounts(idls2,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts([im1 im12],[12:2:86])/sum(histcounts([im1 im12],[12:2:86])),'gaussian',5),'color',[0 0 0],'LineWidth',2)
            plot(smoothdata(histcounts([idls idls2],[12:2:86])/sum(histcounts([idls idls2],[12:2:86])),'gaussian',5),'color',[1 0 0],'LineWidth',2)
            xlim([1 38]); ylim([0 0.2]);
            [h, p, ci, stats] = ttest2([im1 im12],[idls idls2]);
            [p2, ~, ~] = ranksum([im1 im12],[idls idls2]);
            title(['onset | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        subplot(2,2,3); hold on;
            im1 = m1_early_peak;
            idls = dls_early_peak;
            im12 = m1_early_peak2;
            idls2 = dls_early_peak2;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts(im12,[12:2:86])/sum(histcounts(im12,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls2,[12:2:86])/sum(histcounts(idls2,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts([im1 im12],[12:2:86])/sum(histcounts([im1 im12],[12:2:86])),'gaussian',5),'color',[0 0 0],'LineWidth',2)
            plot(smoothdata(histcounts([idls idls2],[12:2:86])/sum(histcounts([idls idls2],[12:2:86])),'gaussian',5),'color',[1 0 0],'LineWidth',2)
            xlim([1 38]); ylim([0 0.2]);
            [h, p, ci, stats] = ttest2([im1 im12],[idls idls2]);
            [p2, ~, ~] = ranksum([im1 im12],[idls idls2]);
            title(['peak | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        subplot(2,2,4); hold on;
            im1 = m1_late_peak;
            idls = dls_late_peak;
            im12 = m1_late_peak2;
            idls2 = dls_late_peak2;
            plot(smoothdata(histcounts(im1,[12:2:86])/sum(histcounts(im1,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls,[12:2:86])/sum(histcounts(idls,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts(im12,[12:2:86])/sum(histcounts(im12,[12:2:86])),'gaussian',5),'color',[0 0 0])
            plot(smoothdata(histcounts(idls2,[12:2:86])/sum(histcounts(idls2,[12:2:86])),'gaussian',5),'color',[1 0 0])
            plot(smoothdata(histcounts([im1 im12],[12:2:86])/sum(histcounts([im1 im12],[12:2:86])),'gaussian',5),'color',[0 0 0],'LineWidth',2)
            plot(smoothdata(histcounts([idls idls2],[12:2:86])/sum(histcounts([idls idls2],[12:2:86])),'gaussian',5),'color',[1 0 0],'LineWidth',2)
            xlim([1 38]); ylim([0 0.2]);
            [h, p, ci, stats] = ttest2([im1 im12],[idls idls2]);
            [p2, ~, ~] = ranksum([im1 im12],[idls idls2]);
            title(['peak | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

        if plot_params.save == 1
            print([plot_params.save_path '\lfp_information_timing.eps'],'-painters','-depsc');
        end

%% NAIVE VS SKILLED

   figure;
    subplot(2,2,1); hold on;
        cdfplot(max(m1_early_info'));
        cdfplot(max(m1_late_info'));
        xlim([0 0.25]);
        grid off;
    subplot(2,2,2); hold on;
        cdfplot(max(dls_early_info'));
        cdfplot(max(dls_late_info'));
        xlim([0 0.25]);
        grid off;
    if plot_params.save == 1
        print([plot_params.save_path '\lfp_info_niave_skilled_lfp.eps'],'-painters','-depsc');
    end   

end

