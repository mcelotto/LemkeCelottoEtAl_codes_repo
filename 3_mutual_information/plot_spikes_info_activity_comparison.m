function [] = plot_spikes_info_activity_comparison(MIout,spikes,reachFeatures,params,clusterParams,save_params,newBinTimes)
%% COMPUTE MAX VEL

    fidx = 9;%f2plot
    
    info_mean_sig.M1 = cell(1,2);
    info_mean_sig.DLS = cell(1,2);
    raw_Spikes.M1 = cell(1,2);
    raw_Spikes.M1{1}=cell(1,1);
    raw_Spikes.M1{2}=cell(1,1);
    raw_Spikes.DLS = cell(1,2);
    raw_Spikes.DLS{1}=cell(1,1);
    raw_Spikes.DLS{2}=cell(1,1);
    reach_features.M1 = cell(1,2);
    reach_features.DLS = cell(1,2);
    
    for animal = 1:numel(params.animals)
        for aidx = 1:length(params.areas)
            % early
            for day = 1:params.num_earlylate_days{animal}
                for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                    infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                    infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                    tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                    infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                    info_mean_sig.(params.areas{aidx}){1} = [info_mean_sig.(params.areas{aidx}){1}; infQuant];
                    raw_Spikes.(params.areas{aidx}){1}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
                    reach_features.(params.areas{aidx}){1}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
                end
            end
            % late
            for day = length(MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MIout{animal}.spikesDay)
                for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                    infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                    infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                    tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                    infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                    info_mean_sig.(params.areas{aidx}){2} = [info_mean_sig.(params.areas{aidx}){2}; infQuant];
                    raw_Spikes.(params.areas{aidx}){2}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
                    reach_features.(params.areas{aidx}){2}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});                    
                end
            end
        end
    end
    
    timeBins4spikes = round(newBinTimes/10)+150;
      
%% PLOT M1 LATE MAX VELOCITY EXAMPLES
        
        %%% PLOT ALL

%             unit_idx = find(sum(info_mean_sig.M1{2}')>0);
%             for unit = unit_idx
%                 all_binned = [];
%                 for trial = 1:size(raw_Spikes.M1{2}{unit+1},1)
%                     tmp_trial = [];
%                     for n = 1:size(timeBins4spikes,1)
%                         tmp_trial = [tmp_trial mean(raw_Spikes.M1{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
%                     end
%                     all_binned = [all_binned; tmp_trial];
%                 end
%                 figure;
%                     subplot(2,1,1);
%                     plot(mean(all_binned))            
%                     subplot(2,1,2);
%                     plot(info_mean_sig.M1{2}(unit,:))
%                     sgtitle(['M1 | late | ' num2str(fidx) ' | unit | ' num2str(unit)])
%             end

        %%% PLOT SPECIFIC EXAMPLES

            units_to_plot = [163 219 227 249 256 266 317 319 322]; % 289 290 328 330 331
            high_info_bins = {[46:53],[40:51],[38:50],[40:52],[44:54],[36:51],[44:53],[33:43],[44:58]};
            high_activity_bins = {[61:68],[54:65],[59:71],[56:68],[63:73],[62:77],[55:64],[53:63],[60:74]};
            count = 1;
            for unit = units_to_plot
                if count ~= 4 
                    count = count + 1;
                    continue
                end
                all_binned = [];
                for trial = 1:size(raw_Spikes.M1{2}{unit+1},1)
                    tmp_trial = [];
                    for n = 1:size(timeBins4spikes,1)
                        tmp_trial = [tmp_trial mean(raw_Spikes.M1{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
                    end
                    all_binned = [all_binned; smoothdata(tmp_trial,'gaussian',10)];
                end
                figure;
                    subplot(3,1,1); hold on;
%                         plot(mean(all_binned)+std(all_binned)/sqrt(size(all_binned,1)))
%                         plot(mean(all_binned)-std(all_binned)/sqrt(size(all_binned,1)))
                        plot(all_binned(1:50,:)','color',[.5 .5 .5]);                        
                        title(['number of trials: ' num2str(size(all_binned,1))])
%                         ylim([0 0.5])
                    subplot(3,1,2); hold on;
                        plot(info_mean_sig.M1{2}(unit,:))
                    subplot(3,1,3); hold on;
                        tmp = mean(all_binned(reach_features.M1{2}{unit}==1,high_activity_bins{count}),2);
                        errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0])
                        tmp = mean(all_binned(reach_features.M1{2}{unit}==2,high_activity_bins{count}),2);
                        errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0])
                        tmp = mean(all_binned(reach_features.M1{2}{unit}==3,high_activity_bins{count}),2);
                        errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0])
                        tmp = mean(all_binned(reach_features.M1{2}{unit}==1,high_info_bins{count}),2);
                        errorbar(4,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0])
                        tmp = mean(all_binned(reach_features.M1{2}{unit}==2,high_info_bins{count}),2);
                        errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0])
                        tmp = mean(all_binned(reach_features.M1{2}{unit}==3,high_info_bins{count}),2);
                        errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0])
%                         ylim([0 0.4])                        
                    sgtitle(['M1 | late | feature ' num2str(fidx) ' | unit | ' num2str(unit)])
                    if save_params.save==1
                        print([save_params.save_path '\exampleMI_M1_feature_' num2str(fidx) '_unit_' num2str(unit) '.eps'],'-painters','-depsc');
                    end
                count = count + 1;
            end

%% PLOT DLS LATE MAX VELOCITY EXAMPLES
        
        %%% PLOT ALL

%             unit_idx = find(sum(info_mean_sig.DLS{2}')>0);
%             for unit = unit_idx
%                 all_binned = [];
%                 for trial = 1:size(raw_Spikes.DLS{2}{unit+1},1)
%                     tmp_trial = [];
%                     for n = 1:size(timeBins4spikes,1)
%                         tmp_trial = [tmp_trial mean(raw_Spikes.DLS{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
%                     end
%                     all_binned = [all_binned; tmp_trial];
%                 end
%                 figure;
%                     subplot(2,1,1);
%                     plot(mean(all_binned))            
%                     subplot(2,1,2);
%                     plot(info_mean_sig.DLS{2}(unit,:))
%                     sgtitle(['DLS | late | ' num2str(fidx) ' | unit | ' num2str(unit)])
%             end

        %%% PLOT SPECIFIC EXAMPLES

            units_to_plot = [131]; % 95 163
            high_info_bins = {[35:45]}; % [52:66] [56:66]
            high_activity_bins = {[50:60]}; % [26:40] [43:53]
            count = 1;
            for unit = units_to_plot
                all_binned = [];
                for trial = 1:size(raw_Spikes.DLS{2}{unit+1},1)
                    tmp_trial = [];
                    for n = 1:size(timeBins4spikes,1)
                        tmp_trial = [tmp_trial mean(raw_Spikes.DLS{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
                    end
                    all_binned = [all_binned; smoothdata(tmp_trial,'gaussian',10)];
                end
                figure;
                    subplot(3,1,1); hold on;
%                         plot(mean(all_binned)+std(all_binned)/sqrt(size(all_binned,1)))
%                         plot(mean(all_binned)-std(all_binned)/sqrt(size(all_binned,1)))
                        plot(all_binned(1:50,:)','color',[.5 .5 .5]);    
                        title(['number of trials: ' num2str(size(all_binned,1))])
%                         ylim([0 0.5])                        
                    subplot(3,1,2); hold on;
                        plot(info_mean_sig.DLS{2}(unit,:))
                    subplot(3,1,3); hold on;
                        tmp = mean(all_binned(reach_features.DLS{2}{unit}==1,high_activity_bins{count}),2);
                        errorbar(1,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0])
                        tmp = mean(all_binned(reach_features.DLS{2}{unit}==2,high_activity_bins{count}),2);
                        errorbar(2,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0])
                        tmp = mean(all_binned(reach_features.DLS{2}{unit}==3,high_activity_bins{count}),2);
                        errorbar(3,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[0 0 0])
                        tmp = mean(all_binned(reach_features.DLS{2}{unit}==1,high_info_bins{count}),2);
                        errorbar(4,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0])
                        tmp = mean(all_binned(reach_features.DLS{2}{unit}==2,high_info_bins{count}),2);
                        errorbar(5,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0])
                        tmp = mean(all_binned(reach_features.DLS{2}{unit}==3,high_info_bins{count}),2);
                        errorbar(6,mean(tmp),std(tmp)/sqrt(length(tmp)),'color',[1 0 0])
                        ylim([0 0.4])                                                
                    sgtitle(['DLS | late | feature ' num2str(fidx) ' | unit | ' num2str(unit)])
            if save_params.save==1
                print([save_params.save_path '\exampleMI_DLS_feature_' num2str(fidx) '_unit_' num2str(unit) '.eps'],'-painters','-depsc');
            end
            count = count + 1;
        end

%% PLOT HISTOGRAMS
    
    both_feature_tmp_max_activity_m1 = [];
    both_feature_tmp_max_variance_m1 = [];
    both_feature_tmp_max_info_m1 = [];
    
    both_feature_tmp_max_activity_dls = [];
    both_feature_tmp_max_variance_dls = [];
    both_feature_tmp_max_info_dls = [];

    all_M1_hist_diff = [];
    all_DLS_hist_diff = [];
    all_M1_hist_diff_var = [];
    all_DLS_hist_diff_var = [];

    %%% MAX VEL
        
        timeBins4spikes = round(newBinTimes/10)+150;
        
        fidx = 9;%f2plot
        info_mean_sig.M1 = cell(1,2);
        info_mean_sig.DLS = cell(1,2);
        info_mean_all.M1 = cell(1,2);
        info_mean_all.DLS = cell(1,2);
        raw_Spikes.M1 = cell(1,2);
        raw_Spikes.M1{1}=cell(1,1);
        raw_Spikes.M1{2}=cell(1,1);
        raw_Spikes.DLS = cell(1,2);
        raw_Spikes.DLS{1}=cell(1,1);
        raw_Spikes.DLS{2}=cell(1,1);
        reach_features.M1 = cell(1,2);
        reach_features.DLS = cell(1,2);
        
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                % early
                for day = 1:params.num_earlylate_days{animal}
                    for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                        infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){1} = [info_mean_all.(params.areas{aidx}){1}; infQuant];
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        info_mean_sig.(params.areas{aidx}){1} = [info_mean_sig.(params.areas{aidx}){1}; infQuant];
                        raw_Spikes.(params.areas{aidx}){1}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
                        reach_features.(params.areas{aidx}){1}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
                    end
                end
                % late
                for day = length(MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MIout{animal}.spikesDay)
                    for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                        infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){2} = [info_mean_all.(params.areas{aidx}){2}; infQuant];
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        info_mean_sig.(params.areas{aidx}){2} = [info_mean_sig.(params.areas{aidx}){2}; infQuant];
                        raw_Spikes.(params.areas{aidx}){2}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
                        reach_features.(params.areas{aidx}){2}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
                    end
                end
            end
        end
    
        tmp_hist_diff = [];
        tmp_hist_diff_var = [];
        tmp_max_info_m1 = [];
        tmp_max_activity_m1 = [];
        tmp_maxvar_activity_m1 = [];
        for unit = 1:size(info_mean_sig.M1{2},1)
            tmp_info = info_mean_sig.M1{2}(unit,:);
            [v,i] = max(tmp_info);  
            if v>0.001 && i>=12 && i<=86
                all_binned = [];
                for trial = 1:size(raw_Spikes.M1{2}{unit+1},1)
                    tmp_trial = [];
                    for n = 1:size(timeBins4spikes,1)
                        tmp_trial = [tmp_trial mean(raw_Spikes.M1{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
                    end
                    all_binned = [all_binned; tmp_trial];
                end
                [~,i_mag] = max(abs(mean(all_binned)));
                [~,i_var] = max(var(all_binned));
                tmp_max_activity_m1 = [tmp_max_activity_m1 i_mag];
                tmp_maxvar_activity_m1 = [tmp_maxvar_activity_m1 i_var];
                tmp_max_info_m1 = [tmp_max_info_m1 i];
                tmp_hist_diff = [tmp_hist_diff i_mag-i];
                tmp_hist_diff_var = [tmp_hist_diff_var i_var-i];
            end
        end
        both_feature_tmp_max_activity_m1 = [both_feature_tmp_max_activity_m1 tmp_max_activity_m1];
        both_feature_tmp_max_variance_m1 = [both_feature_tmp_max_variance_m1 tmp_maxvar_activity_m1];
        both_feature_tmp_max_info_m1 = [both_feature_tmp_max_info_m1 tmp_max_info_m1];
        all_M1_hist_diff = [all_M1_hist_diff tmp_hist_diff];
        all_M1_hist_diff_var = [all_M1_hist_diff_var tmp_hist_diff_var];

        figure;
            subplot(1,2,1);
                histogram(tmp_hist_diff,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)
            subplot(1,2,2);
                histogram(tmp_hist_diff_var,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)            
            sgtitle('M1 max vel late')

        tmp_hist_diff = [];
        tmp_hist_diff_var = [];
        tmp_max_info_dls = [];
        tmp_max_activity_dls = [];
        tmp_maxvar_activity_dls = [];
        for unit = 1:size(info_mean_sig.DLS{2},1)
            tmp_info = info_mean_sig.DLS{2}(unit,:);
            [v,i] = max(tmp_info);  
            if v>0.001 && i>=12 && i<=86
                all_binned = [];
                for trial = 1:size(raw_Spikes.DLS{2}{unit+1},1)
                    tmp_trial = [];
                    for n = 1:size(timeBins4spikes,1)
                        tmp_trial = [tmp_trial mean(raw_Spikes.DLS{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
                    end
                    all_binned = [all_binned; tmp_trial];
                end
                [~,i_mag] = max(abs(mean(all_binned)));
                [~,i_var] = max(var(all_binned));

                tmp_max_activity_dls = [tmp_max_activity_dls i_mag];
                tmp_maxvar_activity_dls = [tmp_maxvar_activity_dls i_var];
                tmp_max_info_dls = [tmp_max_info_dls i];
                tmp_hist_diff = [tmp_hist_diff i_mag-i];
                tmp_hist_diff_var = [tmp_hist_diff_var i_var-i];
            end
        end
        both_feature_tmp_max_activity_dls = [both_feature_tmp_max_activity_dls tmp_max_activity_dls];
        both_feature_tmp_max_variance_dls = [both_feature_tmp_max_variance_dls tmp_maxvar_activity_dls];
        both_feature_tmp_max_info_dls = [both_feature_tmp_max_info_dls tmp_max_info_dls];
        all_DLS_hist_diff = [all_DLS_hist_diff tmp_hist_diff];
        all_DLS_hist_diff_var = [all_DLS_hist_diff_var tmp_hist_diff_var];

        figure;
            subplot(1,2,1);
                histogram(tmp_hist_diff,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)
            subplot(1,2,2);
                histogram(tmp_hist_diff_var,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)            
            sgtitle('DLS max vel late')

    %%% DIST TRAV
    
        fidx = 13;%f2plot
        info_mean_sig.M1 = cell(1,2);
        info_mean_sig.DLS = cell(1,2);
        info_mean_all.M1 = cell(1,2);
        info_mean_all.DLS = cell(1,2);
        raw_Spikes.M1 = cell(1,2);
        raw_Spikes.M1{1}=cell(1,1);
        raw_Spikes.M1{2}=cell(1,1);
        raw_Spikes.DLS = cell(1,2);
        raw_Spikes.DLS{1}=cell(1,1);
        raw_Spikes.DLS{2}=cell(1,1);
        reach_features.M1 = cell(1,2);
        reach_features.DLS = cell(1,2);
        
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                % early
                for day = 1:params.num_earlylate_days{animal}
                    for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                        infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){1} = [info_mean_all.(params.areas{aidx}){1}; infQuant];
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        info_mean_sig.(params.areas{aidx}){1} = [info_mean_sig.(params.areas{aidx}){1}; infQuant];
                        raw_Spikes.(params.areas{aidx}){1}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
                        reach_features.(params.areas{aidx}){1}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
                    end
                end
                % late
                for day = length(MIout{animal}.spikesDay)-params.num_earlylate_days{animal}+1:length(MIout{animal}.spikesDay)
                    for unit = 1:size(MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).info(unit,:);
                        infQuantSh = MIout{animal}.spikesDay{day}.(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(unit,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){2} = [info_mean_all.(params.areas{aidx}){2}; infQuant];
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        info_mean_sig.(params.areas{aidx}){2} = [info_mean_sig.(params.areas{aidx}){2}; infQuant];
                        raw_Spikes.(params.areas{aidx}){2}{end+1} = squeeze(spikes{animal}{day}.(params.areas{aidx})(unit,:,:))';
                        reach_features.(params.areas{aidx}){2}{end+1} = reachFeatures{animal}{day}.(params.reachFeatures{fidx});
                    end
                end
            end
        end
    
        tmp_hist_diff = [];
        tmp_hist_diff_var = [];
        tmp_max_info_m1 = [];
        tmp_max_activity_m1 = [];
        tmp_maxvar_activity_m1 = [];
        for unit = 1:size(info_mean_sig.M1{2},1)
            tmp_info = info_mean_sig.M1{2}(unit,:);
            [v,i] = max(tmp_info);  
            if v>0.001 && i>=12 && i<=86
                all_binned = [];
                for trial = 1:size(raw_Spikes.M1{2}{unit+1},1)
                    tmp_trial = [];
                    for n = 1:size(timeBins4spikes,1)
                        tmp_trial = [tmp_trial mean(raw_Spikes.M1{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
                    end
                    all_binned = [all_binned; tmp_trial];
                end
                [~,i_mag] = max(abs(mean(all_binned)));
                [~,i_var] = max(var(all_binned));
                tmp_max_activity_m1 = [tmp_max_activity_m1 i_mag];
                tmp_maxvar_activity_m1 = [tmp_maxvar_activity_m1 i_var];
                tmp_max_info_m1 = [tmp_max_info_m1 i];
                tmp_hist_diff = [tmp_hist_diff i_mag-i];
                tmp_hist_diff_var = [tmp_hist_diff_var i_var-i];
            end
        end
        both_feature_tmp_max_activity_m1 = [both_feature_tmp_max_activity_m1 tmp_max_activity_m1];
        both_feature_tmp_max_variance_m1 = [both_feature_tmp_max_variance_m1 tmp_maxvar_activity_m1];
        both_feature_tmp_max_info_m1 = [both_feature_tmp_max_info_m1 tmp_max_info_m1];
        all_M1_hist_diff = [all_M1_hist_diff tmp_hist_diff];
        all_M1_hist_diff_var = [all_M1_hist_diff_var tmp_hist_diff_var];

        figure;
            subplot(1,2,1);
                histogram(tmp_hist_diff,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)
            subplot(1,2,2);
                histogram(tmp_hist_diff_var,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)            
            sgtitle('M1 dist trav late')

        tmp_hist_diff = [];
        tmp_hist_diff_var = [];
        tmp_max_info_dls = [];
        tmp_max_activity_dls = [];
        tmp_maxvar_activity_dls = [];
        for unit = 1:size(info_mean_sig.DLS{2},1)
            tmp_info = info_mean_sig.DLS{2}(unit,:);
            [v,i] = max(tmp_info);  
            if v>0.001 && i>=12 && i<=86
                all_binned = [];
                for trial = 1:size(raw_Spikes.DLS{2}{unit+1},1)
                    tmp_trial = [];
                    for n = 1:size(timeBins4spikes,1)
                        tmp_trial = [tmp_trial mean(raw_Spikes.DLS{2}{unit+1}(trial,timeBins4spikes(n,1):timeBins4spikes(n,2)))];
                    end
                    all_binned = [all_binned; tmp_trial];
                end
                [~,i_mag] = max(abs(mean(all_binned)));
                [~,i_var] = max(var(all_binned));

                tmp_max_activity_dls = [tmp_max_activity_dls i_mag];
                tmp_maxvar_activity_dls = [tmp_maxvar_activity_dls i_var];
                tmp_max_info_dls = [tmp_max_info_dls i];
                tmp_hist_diff = [tmp_hist_diff i_mag-i];
                tmp_hist_diff_var = [tmp_hist_diff_var i_var-i];
            end
        end
        both_feature_tmp_max_activity_dls = [both_feature_tmp_max_activity_dls tmp_max_activity_dls];
        both_feature_tmp_max_variance_dls = [both_feature_tmp_max_variance_dls tmp_maxvar_activity_dls];
        both_feature_tmp_max_info_dls = [both_feature_tmp_max_info_dls tmp_max_info_dls];
        all_DLS_hist_diff = [all_DLS_hist_diff tmp_hist_diff];
        all_DLS_hist_diff_var = [all_DLS_hist_diff_var tmp_hist_diff_var];

        figure;
            subplot(1,2,1);
                histogram(tmp_hist_diff,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)
            subplot(1,2,2);
                histogram(tmp_hist_diff_var,[-100:4:100],'DisplayStyle','Stairs')
                xline(0)            
            sgtitle('DLS dist trav late')
     
        tmp_max = [both_feature_tmp_max_activity_m1 both_feature_tmp_max_info_m1];
        all_M1_hist_diff_sh_MAG = [];
        for shuffle = 1:1000
            tmp_sh = randperm(length(tmp_max));
            max_sh = tmp_max(tmp_sh);
            max_sh = reshape(max_sh,2,length(max_sh)/2)';
            all_M1_hist_diff_sh_MAG = [all_M1_hist_diff_sh_MAG; max_sh(:,1)-max_sh(:,2)];
        end
            
        tmp_max = [both_feature_tmp_max_variance_m1 both_feature_tmp_max_info_m1];
        all_M1_hist_diff_sh_VAR = [];
        for shuffle = 1:1000
            tmp_sh = randperm(length(tmp_max));
            max_sh = tmp_max(tmp_sh);
            max_sh = reshape(max_sh,2,length(max_sh)/2)';
            all_M1_hist_diff_sh_VAR = [all_M1_hist_diff_sh_VAR; max_sh(:,1)-max_sh(:,2)];
        end

        tmp_max = [both_feature_tmp_max_activity_dls both_feature_tmp_max_info_dls];
        all_DLS_hist_diff_sh_MAG = [];
        for shuffle = 1:1000
            tmp_sh = randperm(length(tmp_max));
            max_sh = tmp_max(tmp_sh);
            max_sh = reshape(max_sh,2,length(max_sh)/2)';
            all_DLS_hist_diff_sh_MAG = [all_DLS_hist_diff_sh_MAG; max_sh(:,1)-max_sh(:,2)];
        end
            
        tmp_max = [both_feature_tmp_max_variance_dls both_feature_tmp_max_info_dls];
        all_DLS_hist_diff_sh_VAR = [];
        for shuffle = 1:1000
            tmp_sh = randperm(length(tmp_max));
            max_sh = tmp_max(tmp_sh);
            max_sh = reshape(max_sh,2,length(max_sh)/2)';
            all_DLS_hist_diff_sh_VAR = [all_DLS_hist_diff_sh_VAR; max_sh(:,1)-max_sh(:,2)];
        end

        bin_size = 5;

        figure;
            hold on;
            histogram(all_M1_hist_diff,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            histogram(all_M1_hist_diff_sh_MAG,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            [h p] = kstest2(all_M1_hist_diff,all_M1_hist_diff_sh_MAG);
    %         [h p] = ttest2(all_M1_hist_diff,all_M1_hist_diff_sh_MAG)
            xline(0)
            ylim([0 0.15])                
            title(['m1 mag | number of M1 units ' num2str(length(all_M1_hist_diff)) ' |  p= ' num2str(p)])
            if save_params.save==1
                print([save_params.save_path '\M1_Spikes_histogram.eps'],'-painters','-depsc');
            end

        figure;
            hold on;
            histogram(all_M1_hist_diff_var,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            histogram(all_M1_hist_diff_sh_VAR,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            [h p] = kstest2(all_M1_hist_diff_var,all_M1_hist_diff_sh_VAR);
    %         [h p] = ttest2(all_M1_hist_diff_var,all_M1_hist_diff_sh_VAR);
            xline(0)
            ylim([0 0.15])
            title(['m1 var | number of M1 chans ' num2str(length(all_M1_hist_diff_var)) ' |  p= ' num2str(p)])
            if save_params.save==1
                print([save_params.save_path '\M1_Spikes_histogram_var.eps'],'-painters','-depsc');
            end

        figure;
            hold on;
            histogram(all_DLS_hist_diff,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            histogram(all_DLS_hist_diff_sh_MAG,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            [h p] = kstest2(all_DLS_hist_diff,all_DLS_hist_diff_sh_MAG);
            %         [h p] = ttest2(all_DLS_hist_diff,all_DLS_hist_diff_sh_MAG)
            xline(0)
            ylim([0 0.15])
            title(['number of DLS units ' num2str(length(all_DLS_hist_diff)) ' |  p= ' num2str(p)])
            if save_params.save==1
                print([save_params.save_path '\DLS_Spikes_histogram.eps'],'-painters','-depsc');
            end

        figure;
            hold on;
            histogram(all_DLS_hist_diff_var,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            histogram(all_DLS_hist_diff_sh_VAR,[-100:bin_size:100],'DisplayStyle','Stairs','normalization','probability')
            [h p] = kstest2(all_DLS_hist_diff_var,all_DLS_hist_diff_sh_VAR);
    %         [h p] = ttest2(all_DLS_hist_diff_var,all_DLS_hist_diff_sh_VAR);
            xline(0)
            ylim([0 0.15])
            title(['dls var | number of DLS chans ' num2str(length(all_DLS_hist_diff_var)) ' |  p= ' num2str(p)])
            if save_params.save==1
                print([save_params.save_path '\DLS_Spikes_histogram_var.eps'],'-painters','-depsc');
            end

end
