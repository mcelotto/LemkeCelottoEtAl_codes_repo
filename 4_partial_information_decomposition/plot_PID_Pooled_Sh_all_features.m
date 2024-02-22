function [] = plot_PID_Pooled_Sh_all_features(clusterParams,params,reach_bins,plot_params)

animal_colors = distinguishable_colors(numel(params.animals));
n_lfp = 1;

allFeatureInfo_PID_naive = cell(1,6);
allFeatureInfo_PID_skilled = cell(1,6);

%%% TRAJ LENGTH

    PIDpath = 'G:\CurrBiol_revision\data\processed_data\pooled\PID\PID_allFeatures\maxVel_distTrav_100Sh_recRef_070621_traj\';
    tmpFiles = ls([PIDpath 'PID*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'distTrav';

    early_info_all = cell(1,size(tmpFiles,1));
    late_info_all = cell(1,size(tmpFiles,1));
    for animal = 1:size(tmpFiles,1)
        load([PIDpath tmpFiles(animal,:)])
        %%% EARLY
        early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            early_info(pairs,:,:) = tmp_info;
        end
        early_info_all{animal} = early_info;
        late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            late_info(pairs,:,:) = tmp_info;
        end
        late_info_all{animal} = late_info;
    end
    allFeatureInfo_PID_naive{1} = early_info_all;
    allFeatureInfo_PID_skilled{1} = late_info_all;
    
%%% MAX VEL

    PIDpath = 'G:\CurrBiol_revision\data\processed_data\pooled\PID\PID_allFeatures\maxVel_distTrav_100Sh_recRef_070621_traj\';
    tmpFiles = ls([PIDpath 'PID*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxVel';

    early_info_all = cell(1,size(tmpFiles,1));
    late_info_all = cell(1,size(tmpFiles,1));
    for animal = 1:size(tmpFiles,1)
        load([PIDpath tmpFiles(animal,:)])
        %%% EARLY
        early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            early_info(pairs,:,:) = tmp_info;
        end
        early_info_all{animal} = early_info;
        late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            late_info(pairs,:,:) = tmp_info;
        end
        late_info_all{animal} = late_info;
    end
    allFeatureInfo_PID_naive{2} = early_info_all;
    allFeatureInfo_PID_skilled{2} = late_info_all;

%%% MAX ACC

    PIDpath = 'G:\CurrBiol_revision\data\processed_data\pooled\PID\PID_allFeatures\maxAcc_Pos_movDur_100Sh_recRef_070621_traj\';
    tmpFiles = ls([PIDpath 'PID*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxAcc';

    early_info_all = cell(1,size(tmpFiles,1));
    late_info_all = cell(1,size(tmpFiles,1));
    for animal = 1:size(tmpFiles,1)
        load([PIDpath tmpFiles(animal,:)])
        %%% EARLY
        early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            early_info(pairs,:,:) = tmp_info;
        end
        early_info_all{animal} = early_info;
        late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            late_info(pairs,:,:) = tmp_info;
        end
        late_info_all{animal} = late_info;
    end
    allFeatureInfo_PID_naive{3} = early_info_all;
    allFeatureInfo_PID_skilled{3} = late_info_all;

%%% MOV DUR

    PIDpath = 'G:\CurrBiol_revision\data\processed_data\pooled\PID\PID_allFeatures\maxAcc_Pos_movDur_100Sh_recRef_070621_traj\';
    tmpFiles = ls([PIDpath 'PID*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'movDur';

    early_info_all = cell(1,size(tmpFiles,1));
    late_info_all = cell(1,size(tmpFiles,1));
    for animal = 1:size(tmpFiles,1)
        load([PIDpath tmpFiles(animal,:)])
        %%% EARLY
        early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            early_info(pairs,:,:) = tmp_info;
        end
        early_info_all{animal} = early_info;
        late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            late_info(pairs,:,:) = tmp_info;
        end
        late_info_all{animal} = late_info;
    end
    allFeatureInfo_PID_naive{4} = early_info_all;
    allFeatureInfo_PID_skilled{4} = late_info_all;

%%% MAX Y POS

    PIDpath = 'G:\CurrBiol_revision\data\processed_data\pooled\PID\PID_allFeatures\maxAcc_Pos_movDur_100Sh_recRef_070621_traj\';
    tmpFiles = ls([PIDpath 'PID*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxPosY';

    early_info_all = cell(1,size(tmpFiles,1));
    late_info_all = cell(1,size(tmpFiles,1));
    for animal = 1:size(tmpFiles,1)
        load([PIDpath tmpFiles(animal,:)])
        %%% EARLY
        early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            early_info(pairs,:,:) = tmp_info;
        end
        early_info_all{animal} = early_info;
        late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            late_info(pairs,:,:) = tmp_info;
        end
        late_info_all{animal} = late_info;
    end
    allFeatureInfo_PID_naive{5} = early_info_all;
    allFeatureInfo_PID_skilled{5} = late_info_all;

%%% MAX X POS

    PIDpath = 'G:\CurrBiol_revision\data\processed_data\pooled\PID\PID_allFeatures\maxAcc_Pos_movDur_100Sh_recRef_070621_traj\';
    tmpFiles = ls([PIDpath 'PID*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxPosX';

    early_info_all = cell(1,size(tmpFiles,1));
    late_info_all = cell(1,size(tmpFiles,1));
    for animal = 1:size(tmpFiles,1)
        load([PIDpath tmpFiles(animal,:)])
        %%% EARLY
        early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            early_info(pairs,:,:) = tmp_info;
        end
        early_info_all{animal} = early_info;
        late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1),39,51);
        for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared,1)
            tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).shared(pairs,:,:));
            tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).sharedSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            late_info(pairs,:,:) = tmp_info;
        end
        late_info_all{animal} = late_info;
    end
    allFeatureInfo_PID_naive{6} = early_info_all;
    allFeatureInfo_PID_skilled{6} = late_info_all;

%% PLOT
    
    for feature = 1:6
    
        early_info_all = allFeatureInfo_PID_naive{feature};
        late_info_all = allFeatureInfo_PID_skilled{feature};
        
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
                y = [mean(tmp_infoe(:,:,27:51),3); mean(tmp_infoe(:,:,1:25),3)];
                c = [ones(1,size(mean(tmp_infoe(:,:,27:51),3),1)) 2*ones(1,size(mean(tmp_infoe(:,:,1:25),3),1))];
                g(2,1)=gramm('x',x,'y',y,'color',c);
                g(2,1).stat_summary('type','sem');
                g(2,1).axe_property('XLim',[6 36]);
                g(2,1).axe_property('YLim',[0 0.016]);
                g(2,1).set_title(['early']);
                y = [mean(tmp_infol(:,:,27:51),3); mean(tmp_infol(:,:,1:25),3)];
                c = [ones(1,size(mean(tmp_infol(:,:,27:51),3),1)) 2*ones(1,size(mean(tmp_infol(:,:,1:25),3),1))];
                g(2,2)=gramm('x',x,'y',y,'color',c);
                g(2,2).stat_summary('type','sem');
                g(2,2).axe_property('XLim',[6 36]);
                g(2,2).axe_property('YLim',[0 0.016]);
                g(2,2).set_title(['late']);
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\PID_imagesc_feature_' num2str(feature) '.eps'],'-painters','-depsc');
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
            print([plot_params.save_path '\PID_cdf_feature_' num2str(feature) '.eps'],'-painters','-depsc');
        end

    end

end