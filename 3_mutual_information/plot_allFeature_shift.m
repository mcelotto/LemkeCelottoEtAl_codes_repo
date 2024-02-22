function [] = plot_allFeature_shift(MI_LFP,MI_Spikes,MImatched,MIsuccessFail,params,clusterParams,save_params,newBinTimes)

prctile_id = 70;
smooth_factor = 1;
n_lfp = 1;

for max_thresh = [.5 .6 .7 .8 .9]

%% ORIG LFP
    
    all_info_onsets = cell(4,14);
    
    for fidx = 1:14
    
        info_orig.M1 = cell(1,2);
        info_orig.DLS = cell(1,2);
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MI_LFP.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MI_LFP.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MI_LFP.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_orig.(params.areas{aidx}){early_late} = [info_orig.(params.areas{aidx}){early_late}; infQuant];
                    end
                end
            end
        end
    
        timeBins4LFP = round(newBinTimes/10)+150;
    
        m1_early_info = [];
        m1_early_info_all = [];
        for unit = 1:size(info_orig.M1{1},1)
            tmp_info = info_orig.M1{1}(unit,:);
            m1_early_info_all = [m1_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info];
            end
        end
    
        m1_late_info = [];
        m1_late_info_all = [];
        for unit = 1:size(info_orig.M1{2},1)
            tmp_info = info_orig.M1{2}(unit,:);
            m1_late_info_all = [m1_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info];
            end
        end
    
        dls_early_info = [];
        dls_early_info_all = [];
        for unit = 1:size(info_orig.DLS{1},1)
            tmp_info = info_orig.DLS{1}(unit,:);
            %         dls_early_info_all = [dls_early_info_all; tmp_info];
            dls_early_info_all = [dls_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info];
            end
        end
    
        dls_late_info = [];
        dls_late_info_all = [];
        for unit = 1:size(info_orig.DLS{2},1)
            tmp_info = info_orig.DLS{2}(unit,:);
            dls_late_info_all = [dls_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info];
            end
        end
    
        tmp = mean(m1_early_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{1,fidx} = i;
    
        tmp = mean(m1_late_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{2,fidx} = i;
    
        tmp = mean(dls_early_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{3,fidx} = i;
    
        tmp = mean(dls_late_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{4,fidx} = i;
    
    end
    
    early_diff = cellfun(@median,all_info_onsets(1,:))-cellfun(@median,all_info_onsets(3,:));
    late_diff = cellfun(@median,all_info_onsets(2,:))-cellfun(@median,all_info_onsets(4,:));
    
    figure;
    hold on;
    for n = 1:14
        plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
    end
    count = 1;
    for n = 15:28
        plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
        count = count+1;
    end
    sgtitle(['LFP | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
    pause(.1);

%% ORIG SPIKES
    
    all_info_onsets = cell(4,14);
    
    for fidx = 1:14
    
        info_orig.M1 = cell(1,2);
        info_orig.DLS = cell(1,2);
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                % early
                for day = 1:params.num_earlylate_days{animal}
                    for unit = 1:size(MI_Spikes.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                            infQuant = MI_Spikes.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                            infQuantSh = MI_Spikes.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_orig.(params.areas{aidx}){1} = [info_orig.(params.areas{aidx}){1}; infQuant];
                    end
                end
                % late
                for day = length(MI_Spikes.MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MI_Spikes.MIout{animal}.spikesDay)
                    for unit = 1:size(MI_Spikes.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                            infQuant = MI_Spikes.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                            infQuantSh = MI_Spikes.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_orig.(params.areas{aidx}){2} = [info_orig.(params.areas{aidx}){2}; infQuant];
                    end
                end
            end
        end
    
        timeBins4LFP = round(newBinTimes/10)+150;
    
        m1_early_info = [];
        m1_early_info_all = [];
        for unit = 1:size(info_orig.M1{1},1)
            tmp_info = info_orig.M1{1}(unit,:);
            m1_early_info_all = [m1_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info];
            end
        end
    
        m1_late_info = [];
        m1_late_info_all = [];
        for unit = 1:size(info_orig.M1{2},1)
            tmp_info = info_orig.M1{2}(unit,:);
            m1_late_info_all = [m1_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info];
            end
        end
    
        dls_early_info = [];
        dls_early_info_all = [];
        for unit = 1:size(info_orig.DLS{1},1)
            tmp_info = info_orig.DLS{1}(unit,:);
            dls_early_info_all = [dls_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info];
            end
        end
    
        dls_late_info = [];
        dls_late_info_all = [];
        for unit = 1:size(info_orig.DLS{2},1)
            tmp_info = info_orig.DLS{2}(unit,:);
            dls_late_info_all = [dls_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info];
            end
        end
    
        tmp = mean(m1_early_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{1,fidx} = i;
    
        tmp = mean(m1_late_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{2,fidx} = i;
    
        tmp = mean(dls_early_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{3,fidx} = i;
    
        tmp = mean(dls_late_info_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{4,fidx} = i;

%         figure;
%         subplot(2,2,1);
%         hold on;
%         plot(mean(m1_early_info))
%         plot(mean(dls_early_info))
%         subplot(2,2,2);
%         hold on;
%         plot(mean(m1_late_info_all))
%         plot(mean(dls_late_info_all))
%         subplot(2,2,3);
%         hold on;
%         plot(median(m1_early_info))
%         plot(median(dls_early_info))
%         subplot(2,2,4);
%         hold on;
%         plot(median(m1_late_info))
%         plot(median(dls_late_info))
%         sgtitle(num2str(fidx))
    
    end
    
    early_diff = cellfun(@median,all_info_onsets(1,:))-cellfun(@median,all_info_onsets(3,:));
    late_diff = cellfun(@median,all_info_onsets(2,:))-cellfun(@median,all_info_onsets(4,:));
    
    figure;
    hold on;
    for n = 1:14
        plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
    end
    count = 1;
    for n = 15:28
        plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
        count = count+1;
    end
    sgtitle(['Spikes | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
    pause(.1);

% %% MATCHED LFP
%     
%     all_info_onsets = cell(4,14);
%     
%     for fidx = 1:14
%     
%         info_match.M1 = cell(1,2);
%         info_match.DLS = cell(1,2);
%         for animal = 1:numel(params.animals)
%             for aidx = 1:length(params.areas)
%                 for early_late = 1:2
%                     for chan = 1:size(MImatched.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
% 
%                         infQuant = MImatched.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
%                         infQuantSh = MImatched.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
%                         tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
%                         infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
%                         infQuant = infQuant-mean(squeeze(infQuantSh),2)';
%                         info_match.(params.areas{aidx}){early_late} = [info_match.(params.areas{aidx}){early_late}; infQuant];
% 
%                     end
%                 end
%             end
%         end
% 
%         timeBins4LFP = round(newBinTimes/10)+150;
%     
%         m1_early_info = [];
%         m1_early_info_all = [];
%         for unit = 1:size(info_match.M1{1},1)
%             tmp_info = info_match.M1{1}(unit,:);
%             m1_early_info_all = [m1_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 m1_early_info = [m1_early_info; tmp_info];
%             end
%         end
%     
%         m1_late_info = [];
%         m1_late_info_all = [];
%         for unit = 1:size(info_match.M1{2},1)
%             tmp_info = info_match.M1{2}(unit,:);
%             m1_late_info_all = [m1_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 m1_late_info = [m1_late_info; tmp_info];
%             end
%         end
%     
%         dls_early_info = [];
%         dls_early_info_all = [];
%         for unit = 1:size(info_match.DLS{1},1)
%             tmp_info = info_match.DLS{1}(unit,:);
%             %         dls_early_info_all = [dls_early_info_all; tmp_info];
%             dls_early_info_all = [dls_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 dls_early_info = [dls_early_info; tmp_info];
%             end
%         end
%     
%         dls_late_info = [];
%         dls_late_info_all = [];
%         for unit = 1:size(info_match.DLS{2},1)
%             tmp_info = info_match.DLS{2}(unit,:);
%             dls_late_info_all = [dls_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 dls_late_info = [dls_late_info; tmp_info];
%             end
%         end
%     
%         tmp = mean(m1_early_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{1,fidx} = i;
%     
%         tmp = mean(m1_late_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{2,fidx} = i;
%     
%         tmp = mean(dls_early_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{3,fidx} = i;
%     
%         tmp = mean(dls_late_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{4,fidx} = i;
%     
%     end
%     
%     early_diff = cellfun(@median,all_info_onsets(1,:))-cellfun(@median,all_info_onsets(3,:));
%     late_diff = cellfun(@median,all_info_onsets(2,:))-cellfun(@median,all_info_onsets(4,:));
%     
%     figure;
%     hold on;
%     for n = 1:14
%         plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
%     end
%     count = 1;
%     for n = 15:28
%         plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
%         count = count+1;
%     end
%     sgtitle(['Beh Matched LFP | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
%     pause(.1);
% 
% %% MATCHED SPIKES
%     
%     all_info_onsets = cell(4,14);
%     
%     for fidx = 1:14
%     
%         info_match.M1 = cell(1,2);
%         info_match.DLS = cell(1,2);
%         for animal = 1:numel(params.animals)
%             for aidx = 1:length(params.areas)
%                 % early
%                 for day = 1:params.num_earlylate_days{animal}
%                     if isfield(MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}),(params.reachFeatures{fidx}))
%                         for unit = 1:size(MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                             infQuant = MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
%                             infQuantSh = MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
%                             tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
%                             infQuant = infQuant-mean(squeeze(infQuantSh),2)';
%                             infQuant(~tmp_sig) = 0;
%                             info_match.(params.areas{aidx}){1} = [info_match.(params.areas{aidx}){1}; infQuant];
%                         end
%                     end
%                 end
%                 % late
%                 for day = length(MImatched.MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MImatched.MIout{animal}.spikesDay)
%                     if isfield(MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}),(params.reachFeatures{fidx}))                    
%                         for unit = 1:size(MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                             infQuant = MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
%                             infQuantSh = MImatched.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
%                             tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
%                             infQuant = infQuant-mean(squeeze(infQuantSh),2)';
%                             infQuant(~tmp_sig) = 0;
%                             info_match.(params.areas{aidx}){2} = [info_match.(params.areas{aidx}){2}; infQuant];
%                         end
%                     end
%                 end
%             end
%         end
%     
%         timeBins4LFP = round(newBinTimes/10)+150;
%     
%         m1_early_info = [];
%         m1_early_info_all = [];
%         for unit = 1:size(info_match.M1{1},1)
%             tmp_info = info_match.M1{1}(unit,:);
%             m1_early_info_all = [m1_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 m1_early_info = [m1_early_info; tmp_info];
%             end
%         end
%     
%         m1_late_info = [];
%         m1_late_info_all = [];
%         for unit = 1:size(info_match.M1{2},1)
%             tmp_info = info_match.M1{2}(unit,:);
%             m1_late_info_all = [m1_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 m1_late_info = [m1_late_info; tmp_info];
%             end
%         end
%     
%         dls_early_info = [];
%         dls_early_info_all = [];
%         for unit = 1:size(info_match.DLS{1},1)
%             tmp_info = info_match.DLS{1}(unit,:);
%             %         dls_early_info_all = [dls_early_info_all; tmp_info];
%             dls_early_info_all = [dls_early_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 dls_early_info = [dls_early_info; tmp_info];
%             end
%         end
%     
%         dls_late_info = [];
%         dls_late_info_all = [];
%         for unit = 1:size(info_match.DLS{2},1)
%             tmp_info = info_match.DLS{2}(unit,:);
%             dls_late_info_all = [dls_late_info_all; smoothdata(tmp_info,'gaussian',smooth_factor)];
%             [v,i] = max(tmp_info);
%             if v>0 && i>=12 && i<=86
%                 dls_late_info = [dls_late_info; tmp_info];
%             end
%         end
%     
%         tmp = mean(m1_early_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{1,fidx} = i;
%     
%         tmp = mean(m1_late_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{2,fidx} = i;
%     
%         tmp = mean(dls_early_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{3,fidx} = i;
%     
%         tmp = mean(dls_late_info_all);
% %         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
%         tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
%         [~, i] = max(tmp);
%         all_info_onsets{4,fidx} = i;
%     
%     end
%     
%     early_diff = cellfun(@median,all_info_onsets(1,:))-cellfun(@median,all_info_onsets(3,:));
%     late_diff = cellfun(@median,all_info_onsets(2,:))-cellfun(@median,all_info_onsets(4,:));
%     
%     figure;
%     hold on;
%     for n = 1:14
%         plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
%     end
%     count = 1;
%     for n = 15:28
%         plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
%         count = count+1;
%     end
%     sgtitle(['Beh Matched Spikes | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
%     pause(.1);
% 
% 
%% SUCCESS FAIL LFP
    
    all_info_onsets = cell(8,14);
    
    for fidx = 1:14

        info_success.M1 = cell(1,2);
        info_success.DLS = cell(1,2);
        info_fail.M1 = cell(1,2);
        info_fail.DLS = cell(1,2);        
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIsuccessFail.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess,1)

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
    
        tmp = mean(m1_early_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{1,fidx} = i;
    
        tmp = mean(m1_late_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{2,fidx} = i;
    
        tmp = mean(dls_early_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{3,fidx} = i;
    
        tmp = mean(dls_late_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{4,fidx} = i;
    
        tmp = mean(m1_early_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{5,fidx} = i;
    
        tmp = mean(m1_late_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{6,fidx} = i;
    
        tmp = mean(dls_early_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{7,fidx} = i;
    
        tmp = mean(dls_late_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{8,fidx} = i;
    
    end
    
    early_diff = cellfun(@median,all_info_onsets(1,:))-cellfun(@median,all_info_onsets(3,:));
    late_diff = cellfun(@median,all_info_onsets(2,:))-cellfun(@median,all_info_onsets(4,:));
    
    figure;
    hold on;
    for n = 1:14
        plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
    end
    count = 1;
    for n = 15:28
        plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
        count = count+1;
    end
    sgtitle(['LFP success | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
    pause(.1);

    early_diff = cellfun(@median,all_info_onsets(5,:))-cellfun(@median,all_info_onsets(7,:));
    late_diff = cellfun(@median,all_info_onsets(6,:))-cellfun(@median,all_info_onsets(8,:));
    
    figure;
    hold on;
    for n = 1:14
        plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
    end
    count = 1;
    for n = 15:28
        plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
        count = count+1;
    end
    sgtitle(['LFP fail | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
    pause(.1);    

%% SUCCESS FAIL SPIKES
    
    all_info_onsets = cell(4,14);
    
    for fidx = 1:14
       
        info_success.M1 = cell(1,2);
        info_success.DLS = cell(1,2);
        info_fail.M1 = cell(1,2);
        info_fail.DLS = cell(1,2);
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)

                % early
                for day = 1:params.num_earlylate_days{animal}
                    if isfield(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}),'infoSuccess')
                        for unit = 1:size(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess,1)
                            infQuant = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess(unit,:);
                            infQuantSh = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoShSuccess(unit,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_success.(params.areas{aidx}){1} = [info_success.(params.areas{aidx}){1}; infQuant];
                        end
                    end
                    if isfield(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}),'infoFail')
                        for unit = 1:size(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail,1)
                            infQuant = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail(unit,:);
                            infQuantSh = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoShFail(unit,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_fail.(params.areas{aidx}){1} = [info_fail.(params.areas{aidx}){1}; infQuant];
                        end
                    end
                end

                % late
                for day = length(MIsuccessFail.MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MIsuccessFail.MIout{animal}.spikesDay)
                    if isfield(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}),'infoSuccess')
                        for unit = 1:size(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess,1)
                            infQuant = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSuccess(unit,:);
                            infQuantSh = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoShSuccess(unit,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_success.(params.areas{aidx}){2} = [info_success.(params.areas{aidx}){2}; infQuant];
                        end
                    end
                    if isfield(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}),'infoFail')
                        for unit = 1:size(MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail,1)
                            infQuant = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoFail(unit,:);
                            infQuantSh = MIsuccessFail.MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoShFail(unit,:,:);
                            tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                            infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                            infQuant(~tmp_sig) = 0;
                            info_fail.(params.areas{aidx}){2} = [info_fail.(params.areas{aidx}){2}; infQuant];
                        end
                    end
                end
            end
        end
    
        timeBins4LFP = round(newBinTimes/10)+150;

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

        tmp = mean(m1_early_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{1,fidx} = i;
    
        tmp = mean(m1_late_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{2,fidx} = i;
    
        tmp = mean(dls_early_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{3,fidx} = i;
    
        tmp = mean(dls_late_info_success_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{4,fidx} = i;

        tmp = mean(m1_early_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{5,fidx} = i;
    
        tmp = mean(m1_late_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{6,fidx} = i;
    
        tmp = mean(dls_early_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{7,fidx} = i;
    
        tmp = mean(dls_late_info_fail_all);
%         tmp(tmp>prctile(tmp,prctile_id)) = prctile(tmp,prctile_id);
        tmp(tmp>(max(tmp)*max_thresh)) = (max(tmp)*max_thresh);
        [~, i] = max(tmp);
        all_info_onsets{8,fidx} = i;

    end
    
    early_diff = cellfun(@median,all_info_onsets(1,:))-cellfun(@median,all_info_onsets(3,:));
    late_diff = cellfun(@median,all_info_onsets(2,:))-cellfun(@median,all_info_onsets(4,:));
    
    figure;
    hold on;
    for n = 1:14
        plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
    end
    count = 1;
    for n = 15:28
        plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
        count = count+1;
    end
    sgtitle(['LFP success | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
    pause(.1);

    early_diff = cellfun(@median,all_info_onsets(5,:))-cellfun(@median,all_info_onsets(7,:));
    late_diff = cellfun(@median,all_info_onsets(6,:))-cellfun(@median,all_info_onsets(8,:));
    
    figure;
    hold on;
    for n = 1:14
        plot([0 early_diff(n)],[n n],'color',[.5 .5 .5],'LineWidth',2)
    end
    count = 1;
    for n = 15:28
        plot([0 late_diff(count)],[n n],'color',[.5 0 0],'LineWidth',2)
        count = count+1;
    end
    sgtitle(['LFP fail | smooth ' num2str(smooth_factor) ' max_thresh ' num2str(max_thresh)])
    pause(.1);   

end

end
