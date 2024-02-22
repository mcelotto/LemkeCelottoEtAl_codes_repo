function [] = plot_FIT_Pooled_Sh_all_features(clusterParams,params,reach_bins,plot_params)

animal_colors = distinguishable_colors(numel(params.animals));
n_lfp = 1;

allFeatureInfo_FIT_naive_M1receiver = cell(1,6);
allFeatureInfo_FIT_naive_DLSreceiver = cell(1,6);
allFeatureInfo_FIT_skilled_M1receiver = cell(1,6);
allFeatureInfo_FIT_skilled_DLSreceiver = cell(1,6);

%%% TRAJ LENGTH

    FITpath = 'G:\CurrBiol_revision\data\processed_data\pooled\FIT\FIT_allFeatures\distTrav_movDur_070621_traj_100Shuff\';
    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'distTrav';

    early_info_all_M1receiver = cell(1,size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(1,size(tmpFiles,1));
    late_info_all_M1receiver = cell(1,size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(1,size(tmpFiles,1));

    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        %%% EARLY
        early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        early_info_all_M1receiver{animal} = early_info_M1receiver;
        early_info_all_DLSreceiver{animal} = early_info_DLSreceiver;
        %%% LATE
        late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        late_info_all_M1receiver{animal} = late_info_M1receiver;
        late_info_all_DLSreceiver{animal} = late_info_DLSreceiver;
    end
    
    allFeatureInfo_FIT_naive_M1receiver{1} = early_info_all_M1receiver;
    allFeatureInfo_FIT_naive_DLSreceiver{1} = early_info_all_DLSreceiver;
    allFeatureInfo_FIT_skilled_M1receiver{1} = late_info_all_M1receiver;
    allFeatureInfo_FIT_skilled_DLSreceiver{1} = late_info_all_DLSreceiver;

%%% MAX VEL

    FITpath = 'G:\CurrBiol_revision\data\processed_data\pooled\FIT\FIT_allFeatures\maxVel_distTrav\';
    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxVel';

    early_info_all_M1receiver = cell(1,size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(1,size(tmpFiles,1));
    late_info_all_M1receiver = cell(1,size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(1,size(tmpFiles,1));

    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        %%% EARLY
        early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        early_info_all_M1receiver{animal} = early_info_M1receiver;
        early_info_all_DLSreceiver{animal} = early_info_DLSreceiver;
        %%% LATE
        late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        late_info_all_M1receiver{animal} = late_info_M1receiver;
        late_info_all_DLSreceiver{animal} = late_info_DLSreceiver;
    end
    
    allFeatureInfo_FIT_naive_M1receiver{2} = early_info_all_M1receiver;
    allFeatureInfo_FIT_naive_DLSreceiver{2} = early_info_all_DLSreceiver;
    allFeatureInfo_FIT_skilled_M1receiver{2} = late_info_all_M1receiver;
    allFeatureInfo_FIT_skilled_DLSreceiver{2} = late_info_all_DLSreceiver;

%%% MAX ACC

    FITpath = 'G:\CurrBiol_revision\data\processed_data\pooled\FIT\FIT_allFeatures\movDur_maxAcc\';
    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxAcc';

    early_info_all_M1receiver = cell(1,size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(1,size(tmpFiles,1));
    late_info_all_M1receiver = cell(1,size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(1,size(tmpFiles,1));

    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        %%% EARLY
        early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        early_info_all_M1receiver{animal} = early_info_M1receiver;
        early_info_all_DLSreceiver{animal} = early_info_DLSreceiver;
        %%% LATE
        late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        late_info_all_M1receiver{animal} = late_info_M1receiver;
        late_info_all_DLSreceiver{animal} = late_info_DLSreceiver;
    end
    
    allFeatureInfo_FIT_naive_M1receiver{3} = early_info_all_M1receiver;
    allFeatureInfo_FIT_naive_DLSreceiver{3} = early_info_all_DLSreceiver;
    allFeatureInfo_FIT_skilled_M1receiver{3} = late_info_all_M1receiver;
    allFeatureInfo_FIT_skilled_DLSreceiver{3} = late_info_all_DLSreceiver;

%%% MOV DUR

    FITpath = 'G:\CurrBiol_revision\data\processed_data\pooled\FIT\FIT_allFeatures\distTrav_movDur_070621_traj_100Shuff\';
    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'movDur';

    early_info_all_M1receiver = cell(1,size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(1,size(tmpFiles,1));
    late_info_all_M1receiver = cell(1,size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(1,size(tmpFiles,1));

    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        %%% EARLY
        early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        early_info_all_M1receiver{animal} = early_info_M1receiver;
        early_info_all_DLSreceiver{animal} = early_info_DLSreceiver;
        %%% LATE
        late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        late_info_all_M1receiver{animal} = late_info_M1receiver;
        late_info_all_DLSreceiver{animal} = late_info_DLSreceiver;
    end
    
    allFeatureInfo_FIT_naive_M1receiver{4} = early_info_all_M1receiver;
    allFeatureInfo_FIT_naive_DLSreceiver{4} = early_info_all_DLSreceiver;
    allFeatureInfo_FIT_skilled_M1receiver{4} = late_info_all_M1receiver;
    allFeatureInfo_FIT_skilled_DLSreceiver{4} = late_info_all_DLSreceiver;

%%% MAX X POS

    FITpath = 'G:\CurrBiol_revision\data\processed_data\pooled\FIT\FIT_allFeatures\movDur_maxAcc\';
    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxPosX';

    early_info_all_M1receiver = cell(1,size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(1,size(tmpFiles,1));
    late_info_all_M1receiver = cell(1,size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(1,size(tmpFiles,1));

    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        %%% EARLY
        early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        early_info_all_M1receiver{animal} = early_info_M1receiver;
        early_info_all_DLSreceiver{animal} = early_info_DLSreceiver;
        %%% LATE
        late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        late_info_all_M1receiver{animal} = late_info_M1receiver;
        late_info_all_DLSreceiver{animal} = late_info_DLSreceiver;
    end
    
    allFeatureInfo_FIT_naive_M1receiver{5} = early_info_all_M1receiver;
    allFeatureInfo_FIT_naive_DLSreceiver{5} = early_info_all_DLSreceiver;
    allFeatureInfo_FIT_skilled_M1receiver{5} = late_info_all_M1receiver;
    allFeatureInfo_FIT_skilled_DLSreceiver{5} = late_info_all_DLSreceiver;

%%% MAX Y POS

    FITpath = 'G:\CurrBiol_revision\data\processed_data\pooled\FIT\FIT_allFeatures\movDur_maxAcc\';
    tmpFiles = ls([FITpath 'FIT*']);
    tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
    reach_feature = 'maxPosY';

    early_info_all_M1receiver = cell(1,size(tmpFiles,1));
    early_info_all_DLSreceiver = cell(1,size(tmpFiles,1));
    late_info_all_M1receiver = cell(1,size(tmpFiles,1));
    late_info_all_DLSreceiver = cell(1,size(tmpFiles,1));

    for animal = 1:size(tmpFiles,1)
        load([FITpath tmpFiles(animal,:)])
        %%% EARLY
        early_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        early_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            early_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        early_info_all_M1receiver{animal} = early_info_M1receiver;
        early_info_all_DLSreceiver{animal} = early_info_DLSreceiver;
        %%% LATE
        late_info_M1receiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1),39,25);
        late_info_DLSreceiver = zeros(size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver,1),39,25);
        for pairs = 1:size(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver,1)
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_M1receiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_M1receiver(pairs,:,:) = tmp_info;
            tmp_info = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiver(pairs,:,:));
            tmp_infoSh = squeeze(FITout{animal}.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(reach_feature).FIT_negd_DLSreceiverSh(pairs,:,:,:));
            tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
            tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
            %                 tmp_info(~tmp_sig) = 0;
            late_info_DLSreceiver(pairs,:,:) = tmp_info;
        end
        late_info_all_M1receiver{animal} = late_info_M1receiver;
        late_info_all_DLSreceiver{animal} = late_info_DLSreceiver;
    end
    
    allFeatureInfo_FIT_naive_M1receiver{6} = early_info_all_M1receiver;
    allFeatureInfo_FIT_naive_DLSreceiver{6} = early_info_all_DLSreceiver;
    allFeatureInfo_FIT_skilled_M1receiver{6} = late_info_all_M1receiver;
    allFeatureInfo_FIT_skilled_DLSreceiver{6} = late_info_all_DLSreceiver;

%% PLOT
    
    close all

    for feature = 1:6
    
        early_info_all_M1receiver = allFeatureInfo_FIT_naive_M1receiver{feature};
        late_info_all_M1receiver = allFeatureInfo_FIT_skilled_M1receiver{feature};
        early_info_all_DLSreceiver = allFeatureInfo_FIT_naive_DLSreceiver{feature};
        late_info_all_DLSreceiver = allFeatureInfo_FIT_skilled_DLSreceiver{feature};
        
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
                if feature == 1 
                    imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 3e-3])
                elseif feature == 2 
                    imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 3e-3])
                elseif feature == 3 
                    imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 3e-3])
                elseif feature == 4 
                    imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 3e-3])
                elseif feature == 5 
                    imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 4e-3])
                elseif feature == 6 
                    imagesc(squeeze(mean(tmp_infoe_M1receiver,1))',[0 3e-3])
                end
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
            subplot(3,2,3); hold on;
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
                if feature == 1 
                    imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 3e-3])
                elseif feature == 2 
                    imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 3e-3])
                elseif feature == 3 
                    imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 3e-3])
                elseif feature == 4 
                    imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 3e-3])
                elseif feature == 5 
                    imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 4e-3])
                elseif feature == 6 
                    imagesc(squeeze(mean(tmp_infoe_DLSreceiver,1))',[0 3e-3])
                end
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
            subplot(3,2,2); hold on;
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
                if feature == 1 
                    imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 8e-3])
                elseif feature == 2 
                    imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 3e-3])
                elseif feature == 3 
                    imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 3e-3])
                elseif feature == 4 
                    imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 3e-3])
                elseif feature == 5 
                    imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 3e-3])
                elseif feature == 6 
                    imagesc(squeeze(mean(tmp_infol_M1receiver,1))',[0 3e-3])
                end                
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
                if feature == 1 
                    imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 8e-3])
                elseif feature == 2 
                    imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 3e-3])
                elseif feature == 3 
                    imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 3e-3])
                elseif feature == 4 
                    imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 3e-3])
                elseif feature == 5 
                    imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 7e-3])
                elseif feature == 6 
                    imagesc(squeeze(mean(tmp_infol_DLSreceiver,1))',[0 3e-3])
                end 
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
                y = [mean(tmp_infoe_M1receiver,3); mean(tmp_infoe_DLSreceiver,3)];
                c = [ones(1,size(mean(tmp_infoe_M1receiver,3),1)) 2*ones(1,size(mean(tmp_infoe_DLSreceiver,3),1))];
                g(3,1)=gramm('x',x,'y',y,'color',c);
                g(3,1).stat_summary('type','sem');
                g(3,1).axe_property('XLim',[6 36]);
                if feature == 1
                    g(3,1).axe_property('YLim',[0 0.0035]);
                elseif feature == 2
                    g(3,1).axe_property('YLim',[0 0.0035]);
                elseif feature == 3
                    g(3,1).axe_property('YLim',[0 0.0035]);
                elseif feature == 4
                    g(3,1).axe_property('YLim',[0 0.0035]);
                elseif feature == 5
                    g(3,1).axe_property('YLim',[0 0.0035]);
                elseif feature == 6
                    g(3,1).axe_property('YLim',[0 0.0035]);
                end
                g(3,1).set_title(['early']);
                y = [mean(tmp_infol_M1receiver,3); mean(tmp_infol_DLSreceiver,3)];
                c = [ones(1,size(mean(tmp_infol_M1receiver,3),1)) 2*ones(1,size(mean(tmp_infol_DLSreceiver,3),1))];
                g(3,2)=gramm('x',x,'y',y,'color',c);
                g(3,2).stat_summary('type','sem');
                g(3,2).axe_property('XLim',[6 36]);
                if feature == 1
                    g(3,2).axe_property('YLim',[0 0.0035]);
                elseif feature == 2
                    g(3,2).axe_property('YLim',[0 0.0035]);
                elseif feature == 3
                    g(3,2).axe_property('YLim',[0 0.0035]);
                elseif feature == 4
                    g(3,2).axe_property('YLim',[0 0.0035]);
                elseif feature == 5
                    g(3,2).axe_property('YLim',[0 0.0035]);
                elseif feature == 6
                    g(3,2).axe_property('YLim',[0 0.0035]);
                end                
                g(3,2).set_title(['late']);
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\FIT_imagesc_feature_' num2str(feature) '.eps'],'-painters','-depsc');
                end

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
            subplot(1,2,1); hold on;
                hold on;
                cdfplot(max_early_M1toDLS);
                cdfplot(max_late_M1toDLS);
                if feature == 1
                    xlim([0 0.08]);
                elseif feature == 2
                    xlim([0 0.05]);
                elseif feature == 3
                    xlim([0 0.05]);
                elseif feature == 4
                    xlim([0 0.05]);
                elseif feature == 5
                    xlim([0 0.08]);
                elseif feature == 6
                    xlim([0 0.05]);
                end
                [p,~,stats] = ranksum(max_early_M1toDLS,max_early_DLStoM1);
                title(p);
                grid off;
            subplot(1,2,2); hold on;
                cdfplot(max_early_DLStoM1);
                cdfplot(max_late_DLStoM1);
                if feature == 1
                    xlim([0 0.08]);
                elseif feature == 2
                    xlim([0 0.05]);
                elseif feature == 3
                    xlim([0 0.05]);
                elseif feature == 4
                    xlim([0 0.05]);
                elseif feature == 5
                    xlim([0 0.08]);
                elseif feature == 6
                    xlim([0 0.05]);
                end
                grid off;
                [p,~,stats] = ranksum(max_late_M1toDLS,max_late_DLStoM1);
                title(p);
                %
        if plot_params.save == 1
            print([plot_params.save_path '\FIT_cdf_feature_' num2str(feature) '.eps'],'-painters','-depsc');
        end

    end

end