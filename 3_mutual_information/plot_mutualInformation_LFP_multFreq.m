function [] = plot_mutualInformation_LFP_multFreq(params,clusterParams,features_to_plot,save_params,broadband_path,filtered_path,newBinTimes); 

figure('units','normalized','outerposition',[0 0 1 1])

n_lfp = 1;

load('G:\CurrBiol_revision\data\processed_data\pooled\MI\MI_23-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_POOLED.mat');
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
    subplot(2,2,1); hold on;
        plot(nanmean(m1_early_info),'g')
    subplot(2,2,2); hold on;
        plot(nanmean(m1_late_info),'g')
    subplot(2,2,3); hold on;
        plot(nanmean(dls_early_info),'g')
    subplot(2,2,4); hold on;
        plot(nanmean(dls_late_info),'g')

color_count = 1;
for freq = [5 15 30]

    load(['G:\CurrBiol_revision\data\processed_data\pooled\MI\filtered\MI_04-Oct-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_1_lowpass_' num2str(freq) 'hz'])
    info_mean_sig.M1 = cell(2,2);
    info_mean_sig.DLS = cell(2,2);
    info_mean_all.M1 = cell(2,2);
    info_mean_all.DLS = cell(2,2);
    for animal = 1:numel(params.animals)
        for aidx = 1:length(params.areas)
            for early_late = 1:2
                for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{9}).info,1)
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
    subplot(2,2,1); hold on;
        plot(nanmean(m1_early_info),'color',[0 0 (color_count-1)/3])
    subplot(2,2,2); hold on;
        plot(nanmean(m1_late_info),'color',[0 0 (color_count-1)/3])
    subplot(2,2,3); hold on;
        plot(nanmean(dls_early_info),'color',[0 0 (color_count-1)/3])
    subplot(2,2,4); hold on;
        plot(nanmean(dls_late_info),'color',[0 0 (color_count-1)/3])

    color_count = color_count + 1;

end

color_count = 1;
for freq = [5 15 30]

    load(['G:\CurrBiol_revision\data\processed_data\pooled\MI\filtered\MI_07-Oct-2021_LFPearlylate_1_LFPday_0_Spikesearlylate_0_Spikesday_0_SpikesdayAve_0_MatchTrials_1_' num2str(freq) 'Hzhighpass'])
    info_mean_sig.M1 = cell(2,2);
    info_mean_sig.DLS = cell(2,2);
    info_mean_all.M1 = cell(2,2);
    info_mean_all.DLS = cell(2,2);
    for animal = 1:numel(params.animals)
        for aidx = 1:length(params.areas)
            for early_late = 1:2
                for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{9}).info,1)
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
    subplot(2,2,1); hold on;
        plot(nanmean(m1_early_info),'color',[(color_count-1)/3 0 0 ])
    subplot(2,2,2); hold on;
        plot(nanmean(m1_late_info),'color',[(color_count-1)/3 0 0 ])
    subplot(2,2,3); hold on;
        plot(nanmean(dls_early_info),'color',[(color_count-1)/3 0 0 ])
    subplot(2,2,4); hold on;
        plot(nanmean(dls_late_info),'color',[(color_count-1)/3 0 0 ])

    color_count = color_count + 1;

end

if save_params.save == 1
    print([save_params.save_path '\lfp_information_freq_decomp.eps'],'-painters','-depsc');
end

end