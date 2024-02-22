function [] = plot_PID_examples(PIDpath,clusterParams,params,reach_bins,plot_params,PID_params,LFP_early_late,reachFeatures_early_late)

n_lfp = 1;
fidx = 1;
MI_info = load('E:\final_IIT\data\processed_data\pooled\miSh\MI_24-Jun-2021_LFPearlylate_1_LFPday_0_Spikesday_0_MatchTrials_0_POOLED.mat');
sig_channels = cell(numel(params.animals),length(params.lfpFeatures),length(params.reachFeatures),length(params.areas),2);
for animal = 1:numel(params.animals)
    for n_lfp = 1:length(params.lfpFeatures)
        for fidx = 1:length(params.reachFeatures)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    sig_chan = [];
                    for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,PID_params.MIclusterparams(1),PID_params.MIclusterparams(2));
                        infQuant(~tmp_sig) = 0;
                        [~,i] = max(infQuant);
                        sigTime_start = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(1));
                        sigTime_end = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(2));
                        sig_chan = [sig_chan i>=sigTime_start & i<=sigTime_end];
                    end
                    sig_channels{animal,n_lfp,fidx,aidx,early_late} = sig_chan;
                end
            end
        end
    end
end

load('E:\final_IIT\data\processed_data\realBinTimes.mat');
timeBins4LFP = round(newBinTimes/10)+150;

tmpFiles = ls([PIDpath 'PID*']);
tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);

n_lfp = 1;
fidx = 1;

% for animal = 1:numel(params.animals)
%     
%     info_mean_sig.M1 = cell(1,2);
%     info_mean_sig.DLS = cell(1,2);
%     raw_LFP.M1 = cell(1,2);
%     raw_LFP.M1{1}=cell(1,1);
%     raw_LFP.M1{2}=cell(1,1);
%     raw_LFP.DLS = cell(1,2);
%     raw_LFP.DLS{1}=cell(1,1);
%     raw_LFP.DLS{2}=cell(1,1);    
% 
%     for aidx = 1:length(params.areas)
%         for early_late = 1:2
%             for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                 if sig_channels{animal,n_lfp,fidx,aidx,early_late}(chan)==1
%                     infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
%                     infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
%                     tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
%                     infQuant = infQuant-mean(squeeze(infQuantSh),2)';
%                     infQuant(~tmp_sig) = 0;
%                     info_mean_sig.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
%                     raw_LFP.(params.areas{aidx}){early_late}{end+1} =  zscore(squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:)))';
%                 end
%             end
%         end
%     end
% 
%     early_info_all = [];
%     early_info_Orig_all = [];
%     early_info_Sh_all = [];
%     early_pairs_per_animal = 0;
%     late_info_all = [];
%     late_info_Orig_all = [];
%     late_info_Sh_all = [];
%     late_pairs_per_animal = 0;
%     
%     load([PIDpath tmpFiles(animal,:)])
%     
%      %%% EARLY
%      early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
%      early_info_Orig = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
%      early_info_Sh = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
%      for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
%          tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
%          tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
%          tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
%          early_info_Orig(pairs,:,:) = tmp_info;
%          early_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
%          tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%          tmp_info(~tmp_sig) = 0;
%          early_info(pairs,:,:) = tmp_info;
%      end
%      early_info_all = cat(1,early_info_all,early_info);
%      early_info_Orig_all = cat(1,early_info_Orig_all,early_info_Orig);
%      early_info_Sh_all = cat(1,early_info_Sh_all,early_info_Sh);
%      early_pairs_per_animal = [early_pairs_per_animal size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)];
%     
%      %%% LATE
%      late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
%      late_info_Orig = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
%      late_info_Sh = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
%      for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
%          tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
%          tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
%          tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
%          late_info_Orig(pairs,:,:) = tmp_info;
%          late_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
%          tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%          tmp_info(~tmp_sig) = 0;
%          late_info(pairs,:,:) = tmp_info;
%      end
%      late_info_all = cat(1,late_info_all,late_info);
%      late_info_Orig_all = cat(1,late_info_Orig_all,late_info_Orig);
%      late_info_Sh_all = cat(1,late_info_Sh_all,late_info_Sh);
%      late_pairs_per_animal = [late_pairs_per_animal size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)];
% 
%      %%% EARLY
%      for pid_pair = 1:size(early_info_all,1)
%          figure; hold on;
% %          subplot(4,2,[1 3]); hold on;
%              imagesc(squeeze(early_info_all(pid_pair,:,:))')
%              plot([0.5 39.5],[25.5 25.5],'color','r')
%              plot([0.5 39.5],[26.5 26.5],'color','r')
%              plot([26 26],[.5 51.5],'color','k')
%              colorbar;
%              xlim([6 36])
%              xticks([6 16 26 36])
%              xticklabels([-1 -.5 0 .5])
%              xlabel('time from pellet touch (ms)')
%              ylim([.5 51.5])
%              yticks(1:2:51)
%              yticklabels(-250:20:250)
%              ylabel('DLS time delay (ms)')
%              sgtitle(['animal ' num2str(animal) ' early pid pair ' num2str(pid_pair)])
%              
%      end
% 
%      %%% LATE
%      for pid_pair = 1:size(late_info_all,1)
%          figure; hold on;
% %          subplot(4,2,[1 3]); hold on;
%              imagesc(squeeze(late_info_all(pid_pair,:,:))')
%              plot([0.5 39.5],[25.5 25.5],'color','r')
%              plot([0.5 39.5],[26.5 26.5],'color','r')
%              plot([26 26],[.5 51.5],'color','k')
%              colorbar;
%              xlim([6 36])
%              xticks([6 16 26 36])
%              xticklabels([-1 -.5 0 .5])
%              xlabel('time from pellet touch (ms)')
%              ylim([.5 51.5])
%              yticks(1:2:51)
%              yticklabels(-250:20:250)
%              ylabel('DLS time delay (ms)')
%              sgtitle(['animal ' num2str(animal) ' late pid pair ' num2str(pid_pair)])
%      end
% 
%      figure;
%      subplot(1,2,1); hold on;
%      imagesc(squeeze(nanmean(early_info_Orig_all,1))')
%      plot([0.5 39.5],[25.5 25.5],'color','r')
%      plot([0.5 39.5],[26.5 26.5],'color','r')
%      plot([26 26],[.5 51.5],'color','k')
%      colorbar;
%      xlim([6 36])
%      xticks([6 16 26 36])
%      xticklabels([-1 -.5 0 .5])
%      xlabel('time from pellet touch (ms)')
%      ylim([.5 51.5])
%      yticks(1:2:51)
%      yticklabels(-250:20:250)
%      ylabel('DLS time delay (ms)')
%      subplot(1,2,2); hold on;
%      imagesc(squeeze(nanmean(late_info_Orig_all,1))')
%      plot([0.5 39.5],[25.5 25.5],'color','r')
%      plot([0.5 39.5],[26.5 26.5],'color','r')
%      plot([26 26],[.5 51.5],'color','k')
%      colorbar;
%      xlim([6 36])
%      xticks([6 16 26 36])
%      xticklabels([-1 -.5 0 .5])
%      xlabel('time from pellet touch (ms)')
%      ylim([.5 51.5])
%      yticks(1:2:51)
%      yticklabels(-250:20:250)
%      ylabel('DLS time delay (ms)')
% 
%      pause;
%      close all;
% 
% end

animal = 8;
info_mean_sig.M1 = cell(1,2);
info_mean_sig.DLS = cell(1,2);
raw_LFP.M1 = cell(1,2);
raw_LFP.M1{1}=cell(1,1);
raw_LFP.M1{2}=cell(1,1);
raw_LFP.DLS = cell(1,2);
raw_LFP.DLS{1}=cell(1,1);
raw_LFP.DLS{2}=cell(1,1);    
for aidx = 1:length(params.areas)
    for early_late = 1:2
        for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
            if sig_channels{animal,n_lfp,fidx,aidx,early_late}(chan)==1
                infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                tmp_sig = clusterStat_v2(infQuant,infQuantSh,.05,0.95);%clusterParams(1),clusterParams(2));
                infQuant = infQuant-mean(squeeze(infQuantSh),2)';
%                 infQuant(~tmp_sig) = 0;
                info_mean_sig.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
%                 raw_LFP.(params.areas{aidx}){early_late}{end+1} =  zscore(squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:)))';
                raw_LFP.(params.areas{aidx}){early_late}{end+1} =  squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:));
            end
        end
    end
end

early_info_all = [];
early_info_Orig_all = [];
early_info_Sh_all = [];
early_pairs_per_animal = 0;
late_info_all = [];
late_info_Orig_all = [];
late_info_Sh_all = [];
late_pairs_per_animal = 0;

load([PIDpath tmpFiles(animal,:)])

%%% EARLY
early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
early_info_Orig = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
early_info_Sh = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
    tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
    tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
    early_info_Orig(pairs,:,:) = tmp_info;
    early_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%     tmp_info(~tmp_sig) = 0;
    early_info(pairs,:,:) = tmp_info;
end
early_info_all = cat(1,early_info_all,early_info);
early_info_Orig_all = cat(1,early_info_Orig_all,early_info_Orig);
early_info_Sh_all = cat(1,early_info_Sh_all,early_info_Sh);
early_pairs_per_animal = [early_pairs_per_animal size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)];

%%% LATE
late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
late_info_Orig = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
late_info_Sh = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1),39,51);
for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)
    tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared(pairs,:,:));
    tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).sharedSh(pairs,:,:,:));
    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
    late_info_Orig(pairs,:,:) = tmp_info;
    late_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%     tmp_info(~tmp_sig) = 0;
    late_info(pairs,:,:) = tmp_info;
end
late_info_all = cat(1,late_info_all,late_info);
late_info_Orig_all = cat(1,late_info_Orig_all,late_info_Orig);
late_info_Sh_all = cat(1,late_info_Sh_all,late_info_Sh);
late_pairs_per_animal = [late_pairs_per_animal size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).shared,1)];

early_M1_MI = [];
early_DLS_MI = [];
early_M1_LFP = cell(1,1);
early_DLS_LFP = cell(1,1);
for n1 = 1:size(info_mean_sig.M1{1},1)
    for n2 = 1:size(info_mean_sig.DLS{1},1)
        early_M1_MI = [early_M1_MI; info_mean_sig.M1{1}(n1,:)];
        early_DLS_MI = [early_DLS_MI; info_mean_sig.DLS{1}(n2,:)];

        tmp_LFP = [];
        for trial = 1:size(raw_LFP.M1{1}{n1+1},2)
            tmp_trial = [];
            for n = 1:size(timeBins4LFP,1)
                tmp_trial = [tmp_trial mean(raw_LFP.M1{1}{n1+1}(timeBins4LFP(n,1):timeBins4LFP(n,2),trial))];
            end
            tmp_LFP = [tmp_LFP; tmp_trial];
        end
        early_M1_LFP{end+1} = tmp_LFP;

        tmp_LFP = [];
        for trial = 1:size(raw_LFP.DLS{1}{n1+1},2)
            tmp_trial = [];
            for n = 1:size(timeBins4LFP,1)
                tmp_trial = [tmp_trial mean(raw_LFP.DLS{1}{n1+1}(timeBins4LFP(n,1):timeBins4LFP(n,2),trial))];
            end
            tmp_LFP = [tmp_LFP; tmp_trial];
        end
        early_DLS_LFP{end+1} = tmp_LFP;

    end
end

late_M1_MI = [];
late_DLS_MI = [];
late_M1_LFP = cell(1,1);
late_DLS_LFP = cell(1,1);
for n1 = 1:size(info_mean_sig.M1{2},1)
    for n2 = 1:size(info_mean_sig.DLS{2},1)
        late_M1_MI = [late_M1_MI; info_mean_sig.M1{2}(n1,:)];
        late_DLS_MI = [late_DLS_MI; info_mean_sig.DLS{2}(n2,:)];

        tmp_LFP = [];
        for trial = 1:size(raw_LFP.M1{2}{n1+1},2)
            tmp_trial = [];
            for n = 1:size(timeBins4LFP,1)
                tmp_trial = [tmp_trial mean(raw_LFP.M1{2}{n1+1}(timeBins4LFP(n,1):timeBins4LFP(n,2),trial))];
            end
            tmp_LFP = [tmp_LFP; tmp_trial];
        end
        late_M1_LFP{end+1} = tmp_LFP;

        tmp_LFP = [];
        for trial = 1:size(raw_LFP.DLS{2}{n1+1},2)
            tmp_trial = [];
            for n = 1:size(timeBins4LFP,1)
                tmp_trial = [tmp_trial mean(raw_LFP.DLS{2}{n1+1}(timeBins4LFP(n,1):timeBins4LFP(n,2),trial))];
            end
            tmp_LFP = [tmp_LFP; tmp_trial];
        end
        late_DLS_LFP{end+1} = tmp_LFP;

    end
end
 
%%% MEAN VERSION

    %%% EARLY 
    
        M1_MI_bins = {[49:52]};
        DLS_MI_bins = {[59:62]};
        units_plot = [46];
        
%         close all;
        count = 1;
        for early_to_plot = units_plot
        
             figure; hold on;
            
                 subplot(4,2,[1 3]); hold on;
                     imagesc(squeeze(early_info_all(early_to_plot,:,:))',[.005 .02])
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
                     sgtitle(['animal ' num2str(animal) ' early pid pair ' num2str(early_to_plot)])
            
                 subplot(4,2,5); hold on;
                    plot(early_M1_MI(early_to_plot,:))
                    xlim([12 87])
        
                 subplot(4,2,7); hold on;
                    plot(early_DLS_MI(early_to_plot,:))
                    xlim([12 87])
        
                 subplot(4,2,[2 4]); hold on; ylim([2.2 3.5])
                    tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
                    errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
                    errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==3,M1_MI_bins{count}),2);
                    errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
                    tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
                    errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
                    errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==3,M1_MI_bins{count}),2);
                    errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
                    title('MI M1 bins')
        
                 subplot(4,2,[6 8]); hold on; ylim([2.7 3.8])
                    tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
                    errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
                    errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==3,DLS_MI_bins{count}),2);
                    errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
                    tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
                    errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
                    errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==3,DLS_MI_bins{count}),2);
                    errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
                    title('DLS MI DLS bins')            
        
                    print([plot_params.save_path '\examplePIDearly.eps'],'-painters','-depsc');
                    
                count = count + 1;
        end
    
    %%% LATE
        
        M1_MI_bins = {[56:60]};
        DLS_MI_bins = {[42:50]};
        units_plot = [174];
        
        close all;
        count = 1;
        for late_to_plot = units_plot
        
             figure; hold on;
            
                 subplot(4,2,[1 3]); hold on;
                     imagesc(squeeze(late_info_all(late_to_plot,:,:))',[.0 .02])
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
                     sgtitle(['animal ' num2str(animal) ' late pid pair ' num2str(late_to_plot)])
            
                 subplot(4,2,5); hold on;
                    plot(late_M1_MI(late_to_plot,:))
                    xlim([12 87])
        
                 subplot(4,2,7); hold on;
                    plot(late_DLS_MI(late_to_plot,:))
                    xlim([12 87])
        
                 subplot(4,2,[2 4]); hold on; ylim([2.2 3.3])
                    tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
                    errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
                    errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==3,M1_MI_bins{count}),2);
                    errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
                    tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
                    errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
                    errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==3,M1_MI_bins{count}),2);
                    errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
                    title('MI M1 bins')
        
                 subplot(4,2,[6 8]); hold on; ylim([2.6 3.3])
                    tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
                    errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
                    errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
                    tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==3,DLS_MI_bins{count}),2);
                    errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
                    tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
                    errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
                    errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
                    tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==3,DLS_MI_bins{count}),2);
                    errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
                    title('DLS MI DLS bins')             
        
                    print([plot_params.save_path '\examplePIDlate.eps'],'-painters','-depsc');
                    
                count = count + 1;
        end

end

%% OLD
% % 
% % 
% % %1
% % M1_MI_bins = {[56:62],[56:62],[56:62],[56:62],[56:62],[56:62]};
% % DLS_MI_bins = {[49:55],[49:55],[46:55],[46:55],[49:55],[46:55]};
% % units_plot = [39 96 98 117 172 174];
% % 
% % close all;
% % count = 1;
% % for late_to_plot = units_plot
% % 
% %      figure; hold on;
% %     
% %          subplot(4,2,[1 3]); hold on;
% %              imagesc(squeeze(late_info_all(late_to_plot,:,:))',[.01 .02])
% %              plot([0.5 39.5],[25.5 25.5],'color','r')
% %              plot([0.5 39.5],[26.5 26.5],'color','r')
% %              plot([26 26],[.5 51.5],'color','k')
% %              colorbar;
% %              xlim([6 36])
% %              xticks([6 16 26 36])
% %              xticklabels([-1 -.5 0 .5])
% %              xlabel('time from pellet touch (ms)')
% %              ylim([.5 51.5])
% %              yticks(1:2:51)
% %              yticklabels(-250:20:250)
% %              ylabel('DLS time delay (ms)')
% %              sgtitle(['animal ' num2str(animal) ' late pid pair ' num2str(late_to_plot)])
% %     
% %          subplot(4,2,5); hold on;
% %             plot(late_M1_MI(late_to_plot,:))
% %             xlim([12 87])
% % 
% %          subplot(4,2,7); hold on;
% %             plot(late_DLS_MI(late_to_plot,:))
% %             xlim([12 87])
% % 
% %          subplot(4,2,[2 4]); hold on; %ylim([2.8 3.4])
% %             tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
% %             errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %             tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %             errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %             tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %             errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
% %             tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
% %             errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %             tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %             errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %             tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %             errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
% %             title('MI M1 bins')
% % 
% %          subplot(4,2,[6 8]); hold on; %ylim([2.8 3.8])
% %             tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
% %             errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %             tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %             errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %             tmp = mean(late_M1_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %             errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
% %             tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
% %             errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %             tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %             errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %             tmp = mean(late_DLS_LFP{late_to_plot+1}(reachFeatures_early_late{animal}{2}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %             errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
% %             title('DLS MI DLS bins')             
% % 
% % %             print([plot_params.save_path '\examplePIDlate.eps'],'-painters','-depsc');
% %             
% %         count = count + 1;
% % end
% 
% % 
% % for late_to_plot = [39 91 94 95 96 97 98 99 112 113 116 117 118 119 120 121 167 168 171 172 173 174 175 176]
% % 
% %      figure; hold on;
% %     
% %          subplot(4,2,[1 3]); hold on;
% %              imagesc(squeeze(late_info_all(late_to_plot,:,:))')
% %              plot([0.5 39.5],[25.5 25.5],'color','r')
% %              plot([0.5 39.5],[26.5 26.5],'color','r')
% %              plot([26 26],[.5 51.5],'color','k')
% %              colorbar;
% %              xlim([6 36])
% %              xticks([6 16 26 36])
% %              xticklabels([-1 -.5 0 .5])
% %              xlabel('time from pellet touch (ms)')
% %              ylim([.5 51.5])
% %              yticks(1:2:51)
% %              yticklabels(-250:20:250)
% %              ylabel('DLS time delay (ms)')
% %              sgtitle(['animal ' num2str(animal) ' late pid pair ' num2str(late_to_plot)])
% %     
% %          subplot(4,2,5); hold on;
% %             plot(late_M1_MI(late_to_plot,:))
% %             xlim([12 87])
% % 
% %          subplot(4,2,7); hold on;
% %             plot(late_DLS_MI(late_to_plot,:))
% %             xlim([12 87])
% % 
% % end
% 
% % %%% EARLY MIN VERSION
% %     
% %     %1
% %     M1_MI_bins = {[49:52]};
% %     DLS_MI_bins = {[59:62]};
% %     units_plot = [46];
% %     % 
% %     % %2
% %     % M1_MI_bins = {[40:50],[48:55],[48:49]};
% %     % DLS_MI_bins = {[56:63],[56:63],[62:69]};
% %     % units_plot = [41 42 45];
% %     
% %     close all;
% %     count = 1;
% %     for early_to_plot = units_plot
% %     
% %          figure; hold on;
% %         
% %              subplot(4,2,[1 3]); hold on;
% %                  imagesc(squeeze(early_info_all(early_to_plot,:,:))',[.01 .02])
% %                  plot([0.5 39.5],[25.5 25.5],'color','r')
% %                  plot([0.5 39.5],[26.5 26.5],'color','r')
% %                  plot([26 26],[.5 51.5],'color','k')
% %                  colorbar;
% %                  xlim([6 36])
% %                  xticks([6 16 26 36])
% %                  xticklabels([-1 -.5 0 .5])
% %                  xlabel('time from pellet touch (ms)')
% %                  ylim([.5 51.5])
% %                  yticks(1:2:51)
% %                  yticklabels(-250:20:250)
% %                  ylabel('DLS time delay (ms)')
% %                  sgtitle(['animal ' num2str(animal) ' early pid pair ' num2str(early_to_plot)])
% %         
% %              subplot(4,2,5); hold on;
% %                 plot(early_M1_MI(early_to_plot,:))
% %                 xlim([12 87])
% %     
% %              subplot(4,2,7); hold on;
% %                 plot(early_DLS_MI(early_to_plot,:))
% %                 xlim([12 87])
% %     
% %              subplot(4,2,[2 4]); hold on; ylim([2.2 3])
% %     %             tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
% %                 tmp = min(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,M1_MI_bins{count})');
% %                 errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %     %             tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %                 tmp = min(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count})');
% %                 errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %     %             tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %                 tmp = min(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count})');
% %                 errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
% %     %             tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,M1_MI_bins{count}),2);
% %                 tmp = min(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,M1_MI_bins{count})');
% %                 errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %     %             tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %                 tmp = min(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count})');
% %                 errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %     %             tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count}),2);
% %                 tmp = min(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,M1_MI_bins{count})');
% %                 errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
% %                 title('MI M1 bins')
% %     
% %              subplot(4,2,[6 8]); hold on; ylim([2.4 3.4])
% %     %             tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
% %                 tmp = min(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count})');
% %                 errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %     %             tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %                 tmp = min(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count})');
% %                 errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);
% %     %             tmp = mean(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %                 tmp = min(early_M1_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count})');
% %                 errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0]);   
% %                 title('M1 MI DLS bins')
% %     %             tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count}),2);
% %                 tmp = min(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==1,DLS_MI_bins{count})');
% %                 errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %     %             tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %                 tmp = min(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count})');
% %                 errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);
% %     %             tmp = mean(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count}),2);
% %                 tmp = min(early_DLS_LFP{early_to_plot+1}(reachFeatures_early_late{animal}{1}.(params.reachFeatures{fidx})==2,DLS_MI_bins{count})');
% %                 errorbar(7,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0]);   
% %                 title('DLS MI DLS bins')            
% %     
% %                 print([plot_params.save_path '\examplePIDearly.eps'],'-painters','-depsc');
% %                 
% %             count = count + 1;
% %     end
% 
% % 
% 
% % 
% % 
% %         
% 
% % 
% 
% %         for animal = 1:size(tmpFiles,1)
% %             
% % 
% %         end
% % 
% % 
% %         for animal = 1:size(tmpFiles,1)
% % 
% %             % determine channels with significant mutual information
% %             sig_channels = cell(length(params.lfpFeatures),length(params.reachFeatures),length(params.areas),2);
% %             for n_lfp = 1:length(params.lfpFeatures)
% %                 for fidx = 1:length(params.reachFeatures)
% %                     for aidx = 1:length(params.areas)
% %                         for early_late = 1:2
% %                             sig_chan = [];
% %                             for chan = 1:size(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
% %                                 infQuant = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
% %                                 infQuantSh = MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
% %                                 tmp_sig = clusterStat_v2(infQuant,infQuantSh,PID_params.MIclusterparams(1),PID_params.MIclusterparams(2));
% %                                 infQuant(~tmp_sig) = 0;
% %                                 [~,i] = max(infQuant);
% %                                 sigTime_start = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(1));
% %                                 sigTime_end = find(MI_info.MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).binTimes==PID_params.MIsigtiming(2));
% %                                 sig_chan = [sig_chan i>=sigTime_start & i<=sigTime_end];
% %                             end
% %                             sig_channels{animal,n_lfp,fidx,aidx,early_late} = sig_chan;
% %                         end
% %                     end
% %                 end
% %             end
% %         
% %             info_mean_sig.M1 = cell(1,2);
% %             info_mean_sig.DLS = cell(1,2);
% %             raw_Spikes.M1 = cell(1,2);
% %             raw_Spikes.M1{1}=cell(1,1);
% %             raw_Spikes.M1{2}=cell(1,1);
% %             raw_Spikes.DLS = cell(1,2);
% %             raw_Spikes.DLS{1}=cell(1,1);
% %             raw_Spikes.DLS{2}=cell(1,1);
% %             reach_features.M1 = cell(1,2);
% %             reach_features.DLS = cell(1,2);
% %             
% %             for animal = 1:numel(params.animals)
% %                 for aidx = 1:length(params.areas)
% %                     % early
% %                     for day = 1:params.num_earlylate_days{animal}
% %                         for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
% %                             infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
% %                             infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
% %                             tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
% %                             infQuant = infQuant-mean(squeeze(infQuantSh),2)';
% %                             infQuant(~tmp_sig) = 0;
% %                             info_mean_sig.(params.areas{aidx}){1} = [info_mean_sig.(params.areas{aidx}){1}; infQuant];
% %                             raw_Spikes.(params.areas{aidx}){1}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
% %                             reach_features.(params.areas{aidx}){1}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
% %                         end
% %                     end
% %                     % late
% %                     for day = length(MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MIout{animal}.spikesDay)
% %                         for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
% %                             infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
% %                             infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
% %                             tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
% %                             infQuant = infQuant-mean(squeeze(infQuantSh),2)';
% %                             infQuant(~tmp_sig) = 0;
% %                             info_mean_sig.(params.areas{aidx}){2} = [info_mean_sig.(params.areas{aidx}){2}; infQuant];
% %                             raw_Spikes.(params.areas{aidx}){2}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
% %                             reach_features.(params.areas{aidx}){2}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});                    
% %                         end
% %                     end
% %                 end
% %             end
% %             
% %             load('F:\final_IIT\data\processed_data\realBinTimes.mat');
% %             timeBins4spikes = round(newBinTimes/10)+150;
% % 
% % 
% % 
% %         end
% % 
% % 
% %         if fidx==1
% %             ylim_info = [0 0.01];
% %             ylim_info_summary = [0 0.1];
% %             ylim_info_diff = [-.006 .006];
% %         elseif fidx==2
% %             ylim_info = [0 0.02];
% %             ylim_info_summary = [0 0.15];            
% %             ylim_info_diff = [-.015 .015];
% %         end
% %         
% %         %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)
% % 
% %             figure('units','normalized','outerposition',[0 0 1 1]);
% % 
% %                 subplot(2,2,1); hold on;
% %                     imagesc(squeeze(mean(early_info_all,1))',ylim_info)
% %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% %                     plot([26 26],[.5 51.5],'color','k')
% %                     colorbar;
% %                     xlim([6 36])
% %                     xticks([6 16 26 36])
% %                     xticklabels([-1 -.5 0 .5])
% %                     xlabel('time from pellet touch (ms)')
% %                     ylim([.5 51.5])
% %                     yticks(1:2:51)
% %                     yticklabels(-250:20:250)
% %                     ylabel('DLS time delay (ms)')
% %                     title(num2str(size(early_info_all,1)));
% %                 subplot(2,2,2); hold on;
% %                     imagesc(squeeze(mean(late_info_all,1))',ylim_info)
% %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% %                     plot([26 26],[.5 51.5],'color','k')
% %                     colorbar;
% %                     xlim([6 36])
% %                     xticks([6 16 26 36])
% %                     xticklabels([-1 -.5 0 .5])
% %                     xlabel('time from pellet touch (ms)')
% %                     ylim([.5 51.5])
% %                     yticks(1:2:51)
% %                     yticklabels(-250:20:250)
% %                     ylabel('DLS time delay (ms)')
% %                     title(num2str(size(late_info_all,1)));                    
% %                 clearvars g
% %                     x = [1:39];
% %                     y = [mean(early_info_all(:,:,27:51),3); mean(early_info_all(:,:,1:25),3)];
% %                     c = [ones(1,size(mean(early_info_all(:,:,27:51),3),1)) 2*ones(1,size(mean(early_info_all(:,:,1:25),3),1))];
% %                     g(2,1)=gramm('x',x,'y',y,'color',c);
% %                     g(2,1).stat_summary('type','sem');
% %                     g(2,1).axe_property('XLim',[6 36]);
% %                     g(2,1).axe_property('YLim',ylim_info);    
% %                     g(2,1).set_title(['early']);        
% %                     y = [mean(late_info_all(:,:,27:51),3); mean(late_info_all(:,:,1:25),3)];
% %                     c = [ones(1,size(mean(late_info_all(:,:,27:51),3),1)) 2*ones(1,size(mean(late_info_all(:,:,1:25),3),1))];
% %                     g(2,2)=gramm('x',x,'y',y,'color',c);
% %                     g(2,2).stat_summary('type','sem');
% %                     g(2,2).axe_property('XLim',[6 36]);
% %                     g(2,2).axe_property('YLim',ylim_info);    
% %                     g(2,2).set_title(['late']);        
% %                     g.draw();
% %                 if plot_params.save == 1
% %                     print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
% %                 end
% %                     
% %             figure('units','normalized','outerposition',[0 0 1 1]);
% %                 subplot(2,2,1); hold on;
% %                     imagesc(squeeze(mean(late_info_all,1))'-squeeze(mean(early_info_all,1))',ylim_info_diff)
% %                     colormap(redblue)
% %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% %                     plot([26 26],[.5 51.5],'color','k')
% %                     colorbar;
% %                     xlim([6 36])
% %                     xticks([6 16 26 36])
% %                     xticklabels([-1 -.5 0 .5])
% %                     xlabel('time from pellet touch (ms)')
% %                     ylim([.5 51.5])
% %                     yticks(1:2:51)
% %                     yticklabels(-250:20:250)
% %                     ylabel('DLS time delay (ms)')
% %                 subplot(2,2,3); hold on;
% %                     plot(mean(mean(late_info_all(:,:,27:51),3),1)-mean(mean(early_info_all(:,:,27:51),3),1),'color',[.5 1 .25],'LineWidth',2);
% %                     plot(mean(mean(late_info_all(:,:,1:25),3),1)-mean(mean(early_info_all(:,:,1:25),3),1),'color',[0 .25 0],'LineWidth',2);
% %                     xlim([6 36])
% %                     ylim(ylim_info_diff)
% %                     xticks([6 16 26 36])
% %                     xticklabels([-1 -.5 0 .5])
% %                     xlabel('time from pellet touch (ms)')
% %                 if plot_params.save == 1
% %                     print([plot_params.save_path '\PID_sharedInfo_difference_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
% %                 end
% %             
% %         %% HISTOGRAM
% % 
% %             max_early_DLStoM1 = [];
% %             for pair = 1:size(early_info_all,1)
% %                 max_early_DLStoM1 = [max_early_DLStoM1 max(squeeze(early_info_all(pair,reach_bins,1:25)),[],'all')];
% %             end
% %             max_early_M1toDLS = [];
% %             for pair = 1:size(early_info_all,1)
% %                 max_early_M1toDLS = [max_early_M1toDLS max(squeeze(early_info_all(pair,reach_bins,27:51)),[],'all')];
% %             end
% %             max_late_DLStoM1 = [];
% %             for pair = 1:size(late_info_all,1)
% %                 max_late_DLStoM1 = [max_late_DLStoM1 max(squeeze(late_info_all(pair,reach_bins,1:25)),[],'all')];
% %             end
% %             max_late_M1toDLS = [];
% %             for pair = 1:size(late_info_all,1)
% %                 max_late_M1toDLS = [max_late_M1toDLS max(squeeze(late_info_all(pair,reach_bins,27:51)),[],'all')];
% %             end
% % 
% %             figure;
% %                 subplot(2,2,1); hold on;
% %                     histogram(max_early_M1toDLS,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
% %                     histogram(max_early_DLStoM1,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
% %                     [p,~,stats] = signrank(max_early_M1toDLS,max_early_DLStoM1);
% %                     title(p);
% %                     xlim([0 0.1]);
% %                 subplot(2,2,3);
% %                     hold on;
% %                     cdfplot(max_early_M1toDLS);
% %                     cdfplot(max_early_DLStoM1);
% %                     xlim([0 0.1]);
% %                     grid off;
% %                 subplot(2,2,2); hold on;
% %                     histogram(max_late_M1toDLS,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
% %                     histogram(max_late_DLStoM1,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
% %                     [p,~,stats] = signrank(max_late_M1toDLS,max_late_DLStoM1);
% %                     title(p);
% %                     xlim([0 0.1]);
% %                 subplot(2,2,4);
% %                     hold on;
% %                     cdfplot(max_late_M1toDLS);
% %                     cdfplot(max_late_DLStoM1);
% %                     xlim([0 0.1]);
% %                     grid off;
% % 
% %                 if plot_params.save == 1
% %                     print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '_early_late_Histogram.eps'],'-painters','-depsc');
% %                 end   
% % 
% %         %% PLOT DELAYS
% % 
% %             early_peakDelay = cell(1,length(early_pairs_per_animal)-1);
% %             late_peakDelay = cell(1,length(early_pairs_per_animal)-1);
% %             
% %             for animal = 1:length(early_pairs_per_animal)-1
% %                                 
% %                 tmp_peakDelay = [];
% %                 for pair=1+sum(early_pairs_per_animal(1:animal)):sum(early_pairs_per_animal(1:animal+1))
% %                     [~,time_i] = max(squeeze(max(early_info_all(pair,:,:),[],3)));                    
% %                     if ismember(time_i,reach_bins)
% %                         [~, i] = max(squeeze(max(early_info_all(pair,reach_bins,:),[],2)));
% %                         tmp_peakDelay = [tmp_peakDelay i];
% %                     else
% %                         tmp_peakDelay = [tmp_peakDelay nan];
% %                     end
% %                 end
% %                 early_peakDelay{animal} = tmp_peakDelay;
% %                 
% %                 tmp_peakDelay = [];
% %                 for pair = 1+sum(late_pairs_per_animal(1:animal)):sum(late_pairs_per_animal(1:animal+1))
% %                     [~,time_i] = max(squeeze(max(late_info_all(pair,:,:),[],3)));                    
% %                     if ismember(time_i,reach_bins)
% %                         [~, i] = max(squeeze(max(late_info_all(pair,reach_bins,:),[],2)));                        
% %                         tmp_peakDelay = [tmp_peakDelay i];
% %                     else
% %                         tmp_peakDelay = [tmp_peakDelay nan];
% %                     end
% %                 end
% %                 late_peakDelay{animal} = tmp_peakDelay;
% %                                 
% %             end
% % 
% %             figure;
% %                 subplot(1,2,1); hold on;
% %                     histogram([early_peakDelay{:}],[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
% %                     histogram([late_peakDelay{:}],[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
% %                 subplot(1,2,2); hold on;
% %                     tmp_e = histcounts([early_peakDelay{:}],[1:5:51],'normalization','probability');
% %                     tmp_l = histcounts([late_peakDelay{:}],[1:5:51],'normalization','probability');
% %                     scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
% %                     plot(tmp_l-tmp_e,'color',[0 0 0])
% %                     plot([0 11],[0 0],'color','k')
% %                     plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
% %                     xticks([1:1:10])     
% %                     xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
% %                 if plot_params.save == 1
% %                     print([plot_params.save_path '\PID_sharedInfo_delays_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
% %                 end          
% % 
% %             change_pos = [];
% %             change_neg = [];
% %             for animal = 1:length(early_pairs_per_animal)-1
% %                 change_pos = [change_pos sum(late_peakDelay{animal}>=27 & late_peakDelay{animal}<=41)/length(late_peakDelay{animal})-sum(early_peakDelay{animal}>=27 & early_peakDelay{animal}<=41)/length(early_peakDelay{animal})];
% %                 change_neg = [change_neg sum(late_peakDelay{animal}>=6 & late_peakDelay{animal}<=21)/length(late_peakDelay{animal})-sum(early_peakDelay{animal}>=6 & early_peakDelay{animal}<=21)/length(early_peakDelay{animal})];
% %             end         
% %             
% %             figure;
% %                 hold on;
% %                 for animal = 1:length(early_pairs_per_animal)-1
% %                     plot([1 2],[change_neg(animal) change_pos(animal)],'color',animal_colors(animal,:),'LineWidth',1)
% %                 end
% %                 errorbar(1,nanmean(change_neg),nanstd(change_neg)/sqrt(sum(~isnan(change_neg))),'LineWidth',2,'color','k');            
% %                 errorbar(2,nanmean(change_pos),nanstd(change_pos)/sqrt(sum(~isnan(change_pos))),'LineWidth',2,'color','k');
% %                 plot([.5 2.5],[0 0],'color','k','LineStyle','--')
% %                 xlim([.5 2.5])
% %                 [h p1, ~, stats] = ttest(change_pos,change_neg);
% %                 disp(['all animal channel pair change: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                                
% % 
% %                 title([params.reachFeatures{fidx} ' | ' num2str(p1)])
% %             if plot_params.save == 1
% %                 print([plot_params.save_path '\PID_sharedInfo_delays_by_animal_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
% %             end    
% % 
% % 
% % %% OTHER PLOTS
% % % 
% % %         %%% SUMMARY STATS
% % % 
% % %             max_early_neg = [];
% % %             for pair = 1:size(early_info_Orig_all,1)
% % %                 max_early_neg = [max_early_neg max(squeeze(early_info_Orig_all(pair,reach_bins,1:25)),[],'all')];
% % %             end
% % %             max_early_neg_Sh = [];
% % %             for pair = 1:size(early_info_Sh_all,1)
% % %                 max_early_neg_Sh = [max_early_neg_Sh max(squeeze(early_info_Sh_all(pair,reach_bins,1:25)),[],'all')];
% % %             end       
% % %             max_early_pos = [];
% % %             for pair = 1:size(early_info_Orig_all,1)
% % %                 max_early_pos = [max_early_pos max(squeeze(early_info_Orig_all(pair,reach_bins,27:51)),[],'all')];
% % %             end
% % %             max_early_pos_Sh = [];
% % %             for pair = 1:size(early_info_Sh_all,1)
% % %                 max_early_pos_Sh = [max_early_pos_Sh max(squeeze(early_info_Sh_all(pair,reach_bins,27:51)),[],'all')];
% % %             end
% % % 
% % %             max_late_neg = [];
% % %             for pair = 1:size(late_info_Orig_all,1)
% % %                 max_late_neg = [max_late_neg max(squeeze(late_info_Orig_all(pair,reach_bins,1:25)),[],'all')];
% % %             end
% % %             max_late_neg_Sh = [];
% % %             for pair = 1:size(late_info_Orig_all,1)
% % %                 max_late_neg_Sh = [max_late_neg_Sh max(squeeze(late_info_Sh_all(pair,reach_bins,1:25)),[],'all')];
% % %             end       
% % %             max_late_pos = [];
% % %             for pair = 1:size(late_info_Orig_all,1)
% % %                 max_late_pos = [max_late_pos max(squeeze(late_info_Orig_all(pair,reach_bins,27:51)),[],'all')];
% % %             end
% % %             max_late_pos_Sh = [];
% % %             for pair = 1:size(late_info_Orig_all,1)
% % %                 max_late_pos_Sh = [max_late_pos_Sh max(squeeze(late_info_Sh_all(pair,reach_bins,27:51)),[],'all')];
% % %             end
% % % 
% % %             figure('units','normalized','outerposition',[0 0 1 1]);
% % %                 clearvars g    
% % %                 
% % %                 Y = [max_early_neg max_early_neg_Sh max_early_pos max_early_pos_Sh];   
% % %                 X = cell(1,length(Y));
% % %                 X(1:(length(max_early_neg)+length(max_early_neg_Sh))) = {'DLS -> M1'};
% % %                 X(1+(length(max_early_neg)+length(max_early_neg_Sh)):end) = {'M1 -> DLS'};
% % %                 ShReal = [ones(1,length(max_early_neg)) 2*ones(1,length(max_early_neg_Sh)) ones(1,length(max_early_pos)) 2*ones(1,length(max_early_pos_Sh))];        
% % %                 g(1,1)=gramm('x',X,'y',Y,'color',ShReal); 
% % %                 g(1,1).stat_boxplot();
% % %                 [h p1,~,stats] = ttest(max_early_neg,max_early_neg_Sh);
% % %                 disp(['DLS->M1 naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                
% % %                 [h p2,~,stats] = ttest(max_early_pos,max_early_pos_Sh);
% % %                 disp(['M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p2)]);                                
% % %                 neg_tmp = max_early_neg-max_early_neg_Sh;
% % %                 pos_tmp = max_early_pos-max_early_pos_Sh;
% % %                 [h p3,~,stats] = ttest2(neg_tmp,pos_tmp);
% % %                 disp(['DLS->M1 vs. M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p3)]);                                                
% % %                 g(1,1).set_title(['M1->DLS p=' num2str(p2) ' | DLS->M1 p=' num2str(p1) ' | Real M1->DLS vs DLS->M1 p=' num2str(p3)]);
% % %                 g(1,1).axe_property('YLim',ylim_info_summary);
% % % 
% % %                 Y = [max_late_neg max_late_neg_Sh max_late_pos max_late_pos_Sh];   
% % %                 X = cell(1,length(Y));
% % %                 X(1:(length(max_late_neg)+length(max_late_neg_Sh))) = {'DLS -> M1'};
% % %                 X(1+(length(max_late_neg)+length(max_late_neg_Sh)):end) = {'M1 -> DLS'};
% % %                 ShReal = [ones(1,length(max_late_neg)) 2*ones(1,length(max_late_neg_Sh)) ones(1,length(max_late_pos)) 2*ones(1,length(max_late_pos_Sh))];        
% % %                 g(1,2)=gramm('x',X,'y',Y,'color',ShReal); 
% % %                 g(1,2).stat_boxplot();
% % %                 [h p1,~,stats] = ttest(max_late_neg,max_late_neg_Sh);
% % %                 disp(['DLS->M1 naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                
% % %                 [h p2,~,stats] = ttest(max_late_pos,max_late_pos_Sh);
% % %                 disp(['M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p2)]);                                                
% % %                 neg_tmp = max_late_neg-max_late_neg_Sh;
% % %                 pos_tmp = max_late_pos-max_late_pos_Sh;
% % %                 [h p3,~,stats] = ttest2(neg_tmp,pos_tmp);
% % %                 disp(['DLS->M1 vs. M1->DLS naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p3)]);                                                                
% % %                 g(1,2).set_title(['M1->DLS p=' num2str(p2) ' | DLS->M1 p=' num2str(p1) ' | Real M1->DLS vs DLS->M1 p=' num2str(p3)]);
% % %                 g(1,2).axe_property('YLim',ylim_info_summary);    
% % % 
% % %                 g.set_title(['EARLY (left) and LATE (right) | ' params.reachFeatures{fidx}]);                    
% % %                 g.draw();
% % %                 
% % %                 if plot_params.save == 1
% % %                     print([plot_params.save_path '\PID_sharedInfo_summary_stats_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
% % %                 end    
% % %                 
% % %             figure('units','normalized','outerposition',[0 0 1 1]);
% % %                 subplot(3,2,1); hold on;
% % %                     imagesc(squeeze(mean(early_info_all,1))',[0 0.01])
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(3,2,2); hold on;
% % %                     plot(squeeze(mean(mean(early_info_all,1),2))-std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2);
% % %                     plot(squeeze(mean(mean(early_info_all,1),2))+std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2);
% % %                     for ny = 5:5:45
% % %                         if ny == 25
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',2)
% % %                         else
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
% % %                         end
% % %                     end
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(3,2,3); hold on;
% % %                     imagesc(squeeze(mean(late_info_all,1))',[0 0.01])
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(3,2,4); hold on;
% % %                     plot(squeeze(mean(mean(late_info_all,1),2))-std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
% % %                     plot(squeeze(mean(mean(late_info_all,1),2))+std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
% % %                     for ny = 5:5:45
% % %                         if ny == 25
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',2)
% % %                         else
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
% % %                         end
% % %                     end
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(3,2,5); hold on;
% % %                     plot(squeeze(mean(mean(early_info_all,1),2))-std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2);
% % %                     plot(squeeze(mean(mean(early_info_all,1),2))+std(squeeze(mean(early_info_all,1)))'/sqrt(size(early_info_all,1)),[1:51],'color',[.5 1 .25],'LineWidth',2); 
% % %                     plot(squeeze(mean(mean(late_info_all,1),2))-std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
% % %                     plot(squeeze(mean(mean(late_info_all,1),2))+std(squeeze(mean(late_info_all,1)))'/sqrt(size(late_info_all,1)),[1:51],'color',[0 .25 0],'LineWidth',2);
% % %                     for ny = 5:5:45
% % %                         if ny == 25
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',2)
% % %                         else
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
% % %                         end
% % %                     end
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')            
% % %                 subplot(3,2,6); hold on;
% % %                     plot(squeeze(mean(mean(late_info_all,1),2))-squeeze(mean(mean(early_info_all,1),2)),[1:51],'color','k','LineWidth',2);
% % %                     for ny = 5:5:45
% % %                         if ny == 25
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',2)
% % %                         else
% % %                             plot([xlim],[ny ny],'color','k','LineWidth',.5)
% % %                         end
% % %                     end
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')                   
% % 
% %         %%% ORIG VS. SHUFFLED. VS. SIG
% %         
% % %             figure; 
% % %                 subplot(2,3,1); hold on;
% % %                     imagesc(squeeze(mean(early_info_Orig_all,1))'/size(early_info_Orig_all,1));
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(2,3,2); hold on;
% % %                     imagesc(squeeze(mean(early_info_Sh_all,1))'/size(early_info_Sh_all,1));
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(2,3,3); hold on;
% % %                     imagesc(squeeze(mean(early_info_all,1))'/length(find(squeeze(sum(sum(early_info_all,2),3)))))
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(2,3,4); hold on;
% % %                     imagesc(squeeze(mean(late_info_Orig_all,1))'/size(late_info_Orig_all,1));
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(2,3,5); hold on;
% % %                     imagesc(squeeze(mean(late_info_Sh_all,1))'/size(late_info_Sh_all,1));
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % %                 subplot(2,3,6); hold on;
% % %                     imagesc(squeeze(mean(late_info_all,1))'/length(find(squeeze(sum(sum(late_info_all,2),3)))))
% % %                     plot([0.5 39.5],[25.5 25.5],'color','r')
% % %                     plot([0.5 39.5],[26.5 26.5],'color','r')
% % %                     plot([26 26],[.5 51.5],'color','k')
% % %                     colorbar;
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                     ylim([.5 51.5])
% % %                     yticks(1:2:51)
% % %                     yticklabels(-250:20:250)
% % %                     ylabel('DLS time delay (ms)')
% % 
% %         %%% COMPARE EARLY/LATE INFO SIG
% %             
% % %             early_sig_chans = find(squeeze(sum(sum(early_info_all,2),3)));
% % %             late_sig_chans = find(squeeze(sum(sum(late_info_all,2),3)));
% % %             
% % %             figure
% % %                 subplot(2,3,[1 2]); hold on;
% % %                     plot(mean(mean(early_info_all(early_sig_chans,:,27:51),3),1)+(std(mean(early_info_all(early_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(early_info_all(early_sig_chans,:,27:51),3),1)-(std(mean(early_info_all(early_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(late_info_all(late_sig_chans,:,27:51),3),1)+(std(mean(late_info_all(late_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
% % %                     plot(mean(mean(late_info_all(late_sig_chans,:,27:51),3),1)-(std(mean(late_info_all(late_sig_chans,:,27:51),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
% % %                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
% % %                     title(['early/late sig M1 -> DLS shared information | ' params.reachFeatures{fidx}])
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                 subplot(2,3,[4 5]); hold on;
% % %                     plot(mean(mean(early_info_all(early_sig_chans,:,1:25),3),1)+(std(mean(early_info_all(early_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(early_info_all(early_sig_chans,:,1:25),3),1)-(std(mean(early_info_all(early_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(late_info_all(late_sig_chans,:,1:25),3),1)+(std(mean(late_info_all(late_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
% % %                     plot(mean(mean(late_info_all(late_sig_chans,:,1:25),3),1)-(std(mean(late_info_all(late_sig_chans,:,1:25),3),1)/sqrt(length(early_sig_chans))),'b','LineWidth',2);
% % %                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
% % %                     title(['early/late sig DLS -> M1 shared information | ' params.reachFeatures{fidx}])
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                 subplot(2,3,3); hold on;
% % %                     errorbar(1,mean(max(mean(early_info_all(early_sig_chans,6:36,27:51),3)')),std(max(mean(early_info_all(early_sig_chans,6:36,27:51),3)'))/sqrt(length(early_sig_chans)),'LineWidth',2,'color',[.5 .5 .5])
% % %                     errorbar(2,mean(max(mean(late_info_all(late_sig_chans,6:36,27:51),3)')),std(max(mean(late_info_all(late_sig_chans,6:36,27:51),3)'))/sqrt(length(late_sig_chans)),'LineWidth',2,'color','b')
% % %                     [h p] = ttest2(max(mean(early_info_all(early_sig_chans,6:36,27:51),3)'),max(mean(late_info_all(late_sig_chans,6:36,27:51),3)'));
% % %                     title(num2str(p))
% % %                     xlim([.5 2.5])
% % %                 subplot(2,3,6); hold on;
% % %                     errorbar(1,mean(max(mean(early_info_all(early_sig_chans,6:36,1:25),3)')),std(max(mean(early_info_all(early_sig_chans,6:36,1:25),3)'))/sqrt(length(early_sig_chans)),'LineWidth',2,'color',[.5 .5 .5])
% % %                     errorbar(2,mean(max(mean(late_info_all(late_sig_chans,6:36,1:25),3)')),std(max(mean(late_info_all(late_sig_chans,6:36,1:25),3)'))/sqrt(length(late_sig_chans)),'LineWidth',2,'color','b')
% % %                     [h p] = ttest2(max(mean(early_info_all(early_sig_chans,6:36,1:25),3)'),max(mean(late_info_all(late_sig_chans,6:36,1:25),3)'));
% % %                     title(num2str(p))
% % %                     xlim([.5 2.5])
% %                
% %         %%% COMPARE EARLY/LATE INFO ALL
% % 
% % %             figure
% % %                 subplot(2,3,[1 2]); hold on;
% % %                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)+(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)-(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)+(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)-(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
% % %                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
% % %                     title(['early/late sig M1 -> DLS shared information | ' params.reachFeatures{fidx}])
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                 subplot(2,3,[4 5]); hold on;
% % %                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)+(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)-(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 .5 .5],'LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)+(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)-(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'b','LineWidth',2);
% % %                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
% % %                     title(['early/late sig DLS -> M1 shared information | ' params.reachFeatures{fidx}])
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                 subplot(2,3,3); hold on;
% % %                     errorbar(1,mean(max(mean(early_info_Orig_all(:,:,27:51),3)')),std(max(mean(early_info_Orig_all(:,:,27:51),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[.5 .5 .5])
% % %                     errorbar(2,mean(max(mean(late_info_Orig_all(:,:,27:51),3)')),std(max(mean(late_info_Orig_all(:,:,27:51),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color','b')
% % %                     [h p] = ttest2(max(mean(early_info_Orig_all(:,:,27:51),3)'),max(mean(late_info_Orig_all(:,:,27:51),3)'));
% % %                     title(num2str(p))
% % %                     xlim([.5 2.5])
% % %                 subplot(2,3,6); hold on;
% % %                     errorbar(1,mean(max(mean(early_info_Orig_all(:,:,1:25),3)')),std(max(mean(early_info_Orig_all(:,:,1:25),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[.5 .5 .5])
% % %                     errorbar(2,mean(max(mean(late_info_Orig_all(:,:,1:25),3)')),std(max(mean(late_info_Orig_all(:,:,1:25),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color','b')
% % %                     [h p] = ttest2(max(mean(early_info_Orig_all(:,:,1:25),3)'),max(mean(late_info_Orig_all(:,:,1:25),3)'));
% % %                     title(num2str(p))
% % %                     xlim([.5 2.5])
% % 
% %         %%% COMPARE M1->DLS/DLS->M1 INFO ALL
% %     
% % %             figure
% % %                 subplot(2,3,[1 2]); hold on;
% % %                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)+(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
% % %                     plot(mean(mean(early_info_Orig_all(:,:,27:51),3),1)-(std(mean(early_info_Orig_all(:,:,27:51),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
% % %                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)+(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
% % %                     plot(mean(mean(early_info_Orig_all(:,:,1:25),3),1)-(std(mean(early_info_Orig_all(:,:,1:25),3),1)/sqrt(size(early_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
% % %                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
% % %                     title(['early sig shared information | ' params.reachFeatures{fidx}])
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                 subplot(2,3,[4 5]); hold on;
% % %                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)+(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,27:51),3),1)-(std(mean(late_info_Orig_all(:,:,27:51),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[.5 1 .25],'LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)+(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
% % %                     plot(mean(mean(late_info_Orig_all(:,:,1:25),3),1)-(std(mean(late_info_Orig_all(:,:,1:25),3),1)/sqrt(size(late_info_Orig_all,1))),'color',[0 .25 0],'LineWidth',2);
% % %                     plot([26 26],[ylim],'color','k','LineWidth',2,'LineStyle','--')
% % %                     title(['late sig shared information | ' params.reachFeatures{fidx}])
% % %                     xlim([6 36])
% % %                     xticks([6 16 26 36])
% % %                     xticklabels([-1 -.5 0 .5])
% % %                     xlabel('time from pellet touch (ms)')
% % %                 subplot(2,3,3); hold on;
% % %                     errorbar(1,mean(max(mean(early_info_Orig_all(:,:,27:51),3)')),std(max(mean(early_info_Orig_all(:,:,27:51),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[.5 1 .25])
% % %                     errorbar(2,mean(max(mean(early_info_Orig_all(:,:,1:25),3)')),std(max(mean(early_info_Orig_all(:,:,1:25),3)'))/sqrt(size(early_info_Orig_all,1)),'LineWidth',2,'color',[0 .25 0])
% % %                     [h p] = ttest2(max(mean(early_info_Orig_all(:,:,27:51),3)'),max(mean(early_info_Orig_all(:,:,1:25),3)'));
% % %                     title(num2str(p))
% % %                     xlim([.5 2.5])
% % %                 subplot(2,3,6); hold on;
% % %                     errorbar(1,mean(max(mean(late_info_Orig_all(:,:,27:51),3)')),std(max(mean(late_info_Orig_all(:,:,27:51),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color',[.5 1 .25])
% % %                     errorbar(2,mean(max(mean(late_info_Orig_all(:,:,1:25),3)')),std(max(mean(late_info_Orig_all(:,:,1:25),3)'))/sqrt(size(late_info_Orig_all,1)),'LineWidth',2,'color',[0 .25 0])
% % %                     [h p] = ttest2(max(mean(late_info_Orig_all(:,:,27:51),3)'),max(mean(late_info_Orig_all(:,:,1:25),3)'));
% % %                     title(num2str(p))
% % %                     xlim([.5 2.5])
% %        
% %     end
% % end
% % 
% % end
% % 
% % function c = redblue(m)
% % %REDBLUE    Shades of red and blue color map
% % %   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
% % %   The colors begin with bright blue, range through shades of
% % %   blue to white, and then through shades of red to bright red.
% % %   REDBLUE, by itself, is the same length as the current figure's
% % %   colormap. If no figure exists, MATLAB creates one.
% % %
% % %   For example, to reset the colormap of the current figure:
% % %
% % %             colormap(redblue)
% % %
% % %   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
% % %   COLORMAP, RGBPLOT.
% % 
% % %   Adam Auton, 9th October 2009
% % 
% % if nargin < 1, m = size(get(gcf,'colormap'),1); end
% % 
% % if (mod(m,2) == 0)
% %     % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
% %     m1 = m*0.5;
% %     r = (0:m1-1)'/max(m1-1,1);
% %     g = r;
% %     r = [r; ones(m1,1)];
% %     g = [g; flipud(g)];
% %     b = flipud(r);
% % else
% %     % From [0 0 1] to [1 1 1] to [1 0 0];
% %     m1 = floor(m*0.5);
% %     r = (0:m1-1)'/max(m1,1);
% %     g = r;
% %     r = [r; ones(m1+1,1)];
% %     g = [g; 1; flipud(g)];
% %     b = flipud(r);
% % end
% % 
% % c = [r g b]; 
% % 
% % end
