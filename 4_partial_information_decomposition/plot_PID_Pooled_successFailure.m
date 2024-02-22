function [] = plot_PID_Pooled_successFailure(PIDpath,clusterParams,features_to_plot,params,reach_bins,plot_params)
%% Add pooled features plot

% Plot statistics for different trials subgroups 
% Adding MI reconstruction from PID
animal_colors = distinguishable_colors(numel(params.animals));
trial_group_names = {'success','failure'};
epsilon = 10^(-9);

for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};

    early_info_all_pooled.(tgLabel) = [];
    early_info_Orig_all_pooled.(tgLabel) = [];
    early_info_Sh_all_pooled.(tgLabel) = [];
    early_MI_all_pooled.(tgLabel).M1 = [];
    early_MI_all_pooled.(tgLabel).DLS = [];
    early_pairs_per_animal_pooled.(tgLabel) = zeros(1,9);
    late_info_all_pooled.(tgLabel) = [];
    late_info_Orig_all_pooled.(tgLabel) = [];
    late_info_Sh_all_pooled.(tgLabel) = [];
    late_MI_all_pooled.(tgLabel).M1 = [];
    late_MI_all_pooled.(tgLabel).DLS = [];
    late_pairs_per_animal_pooled.(tgLabel) = zeros(1,9);
end


for n_lfp = 1:length(params.lfpFeatures)
    for fidx = features_to_plot
        %%
        disp(params.reachFeatures{fidx})
        
        % Initialize information structures
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            
            early_info_all.(tgLabel) = [];
            early_info_Orig_all.(tgLabel) = [];
            early_info_Sh_all.(tgLabel) = [];
            early_MI_all.(tgLabel).M1 = [];
            early_MI_all.(tgLabel).DLS = [];
            early_pairs_per_animal.(tgLabel) = 0;
            late_info_all.(tgLabel) = [];
            late_info_Orig_all.(tgLabel) = [];
            late_info_Sh_all.(tgLabel) = [];
            late_MI_all.(tgLabel).M1 = [];
            late_MI_all.(tgLabel).DLS = [];
            late_pairs_per_animal.(tgLabel) = 0;
        end
        
        tmpFiles = dir([PIDpath 'PID*']);        
    %%
    aIdx = 0;
        for animal = [3 4 5 6 7 1 2 8]
            aIdx = aIdx+1;
            disp(['Animal ',num2str(animal)])
            load([PIDpath tmpFiles(animal).name],'PIDout')
            for trialGroupIdx = 1:length(trial_group_names)
                tgLabel = trial_group_names{trialGroupIdx};
                %%% EARLY
                tmp_PID = PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).(trial_group_names{trialGroupIdx});

                early_info = zeros(size(tmp_PID.shared,1),39,51);
                early_info_Orig = zeros(size(tmp_PID.shared,1),39,51);
                early_info_Sh = zeros(size(tmp_PID.shared,1),39,51);
                early_MI_M1 = early_info; early_MI_DLS = early_info;
                for pairs = 1:size(tmp_PID.shared,1)
                    tmp_info = squeeze(tmp_PID.shared(pairs,:,:));
                    tmp_infoSh = squeeze(tmp_PID.sharedSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    early_info_Orig(pairs,:,:) = tmp_info;
                    early_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
%                     tmp_info(~tmp_sig) = 0;
                    early_info(pairs,:,:) = tmp_info;
                    
                    tmpMI_m1_shuff = mean(tmp_PID.sharedSh(pairs,:,:,:) + tmp_PID.unique_m1Sh(pairs,:,:,:),4);
                    early_MI_M1(pairs,:,:) = (tmp_PID.shared(pairs,:,:) + tmp_PID.unique_m1(pairs,:,:))-tmpMI_m1_shuff;
                    tmpMI_dls_shuff = mean(tmp_PID.sharedSh(pairs,:,:,:) + tmp_PID.unique_dlsSh(pairs,:,:,:),4);
                    early_MI_DLS(pairs,:,:) = tmp_PID.shared(pairs,:,:) + tmp_PID.unique_dls(pairs,:,:)-tmpMI_dls_shuff;
                end
                early_info_all.(tgLabel) = cat(1,early_info_all.(tgLabel),early_info);
                early_info_Orig_all.(tgLabel) = cat(1,early_info_Orig_all.(tgLabel),early_info_Orig);
                early_info_Sh_all.(tgLabel) = cat(1,early_info_Sh_all.(tgLabel),early_info_Sh);
                early_pairs_per_animal.(tgLabel) = [early_pairs_per_animal.(tgLabel) size(tmp_PID.shared,1)];
                early_MI_all.(tgLabel).M1 = cat(1,early_MI_all.(tgLabel).M1,early_MI_M1);
                early_MI_all.(tgLabel).DLS = cat(1,early_MI_all.(tgLabel).DLS,early_MI_DLS);

                %%% LATE
                tmp_PID = PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.reachFeatures{fidx}).(trial_group_names{trialGroupIdx});

                late_info = zeros(size(tmp_PID.shared,1),39,51);
                late_info_Orig = zeros(size(tmp_PID.shared,1),39,51);
                late_info_Sh = zeros(size(tmp_PID.shared,1),39,51);
                late_MI_M1 = late_info; late_MI_DLS = late_info;

                for pairs = 1:size(tmp_PID.shared,1)
                    tmp_info = squeeze(tmp_PID.shared(pairs,:,:));
                    tmp_infoSh = squeeze(tmp_PID.sharedSh(pairs,:,:,:));
                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                    late_info_Orig(pairs,:,:) = tmp_info;
                    late_info_Sh(pairs,:,:) = squeeze(mean(tmp_infoSh,3));
                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));                
%                     tmp_info(~tmp_sig) = 0;
                    late_info(pairs,:,:) = tmp_info;
                    
                    tmpMI_m1_shuff = mean(tmp_PID.sharedSh(pairs,:,:,:) + tmp_PID.unique_m1Sh(pairs,:,:,:),4);
                    late_MI_M1(pairs,:,:) = tmp_PID.shared(pairs,:,:) + tmp_PID.unique_m1(pairs,:,:) - tmpMI_m1_shuff;
                    tmpMI_dls_shuff = mean(tmp_PID.sharedSh(pairs,:,:,:) + tmp_PID.unique_dlsSh(pairs,:,:,:),4);
                    late_MI_DLS(pairs,:,:) = tmp_PID.shared(pairs,:,:) + tmp_PID.unique_dls(pairs,:,:) - tmpMI_dls_shuff;
                end
                late_info_all.(tgLabel) = cat(1,late_info_all.(tgLabel),late_info);
                late_info_Orig_all.(tgLabel) = cat(1,late_info_Orig_all.(tgLabel),late_info_Orig);
                late_info_Sh_all.(tgLabel) = cat(1,late_info_Sh_all.(tgLabel),late_info_Sh);
                late_pairs_per_animal.(tgLabel) = [late_pairs_per_animal.(tgLabel) size(tmp_PID.shared,1)];
                late_MI_all.(tgLabel).M1 = cat(1,late_MI_all.(tgLabel).M1,late_MI_M1);
                late_MI_all.(tgLabel).DLS = cat(1,late_MI_all.(tgLabel).DLS,late_MI_DLS);
            
               
            end
            
            if fidx==9
                ylim_info = [0 0.015];
                ylim_info_summary = [0 0.1];
                ylim_info_diff = [-.006 .006];
            elseif fidx==13
                ylim_info = [0 0.02];
                ylim_info_summary = [0 0.15];            
                ylim_info_diff = [-.015 .015];
            else
                ylim_info = [0 0.015];
                ylim_info_summary = [0 0.15];            
                ylim_info_diff = [-.015 .015];
            end
        end
        %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)
        disp('Plot spatiotemporal maps')
        figure('units','normalized','outerposition',[0 0 1 1]);
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
                
            % Remove pairs for animal T107 early (only 2 success trials
            % --> the shuffled info is equal to the original one and it
            % introduces zeros)
            nanpairs = find(isnan(squeeze(mean(mean(early_info_all.(tgLabel),2),3))) | (abs(squeeze(mean(mean(early_info_all.(tgLabel),2),3)))<epsilon));
            early_info_all.(tgLabel)(nanpairs,:,:) = [];
            early_MI_all.(tgLabel).M1(nanpairs,:,:) = [];
            early_MI_all.(tgLabel).DLS(nanpairs,:,:) = [];
            tmp = early_pairs_per_animal.(tgLabel);
            rem_animal = min(find(cumsum(tmp)>nanpairs(1)));
            tmp(rem_animal) = 0;
            early_pairs_per_animal.(tgLabel) = tmp;
            
            subplot(2,2,(trialGroupIdx*2)-1); hold on;
            imagesc(squeeze(mean(early_info_all.(tgLabel),1))',ylim_info)
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
            title([tgLabel,' early; ',num2str(size(early_info_all.(tgLabel),1)),' pairs']);
                
            subplot(2,2,trialGroupIdx*2); hold on;
            imagesc(squeeze(mean(late_info_all.(tgLabel),1))',ylim_info)
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
            title([tgLabel,' late; ',num2str(size(late_info_all.(tgLabel),1)),' pairs']); 
            
            % Pool features (doing it here, after removing early pairs)
            early_info_all_pooled.(tgLabel) = cat(1,early_info_all_pooled.(tgLabel),early_info_all.(tgLabel));
            early_MI_all_pooled.(tgLabel).M1 = cat(1,early_MI_all_pooled.(tgLabel).M1,early_MI_all.(tgLabel).M1);
            early_MI_all_pooled.(tgLabel).DLS = cat(1,early_MI_all_pooled.(tgLabel).DLS,early_MI_all.(tgLabel).DLS);
            early_pairs_per_animal_pooled.(tgLabel) = early_pairs_per_animal_pooled.(tgLabel) + early_pairs_per_animal.(tgLabel);

            late_info_all_pooled.(tgLabel) = cat(1,late_info_all_pooled.(tgLabel),late_info_all.(tgLabel));
            late_MI_all_pooled.(tgLabel).M1 = cat(1,late_MI_all_pooled.(tgLabel).M1,late_MI_all.(tgLabel).M1);
            late_MI_all_pooled.(tgLabel).DLS = cat(1,late_MI_all_pooled.(tgLabel).DLS,late_MI_all.(tgLabel).DLS);
            late_pairs_per_animal_pooled.(tgLabel) = late_pairs_per_animal_pooled.(tgLabel) + late_pairs_per_animal.(tgLabel);
        end
        sgtitle([params.reachFeatures{fidx}, ' time-delay maps'])

        % Info time profiles
        figure('units','normalized','outerposition',[0 0 1 1]);
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};

            subplot(2,2,1+2*(trialGroupIdx-1)); hold on;
            plot(squeeze(mean(early_MI_all.(tgLabel).M1(:,:,26),1))')
            plot(squeeze(mean(early_MI_all.(tgLabel).DLS(:,:,26),1))')
            xlim([6 36])
            xticks([6 16 26 36])
            xticklabels([-1 -.5 0 .5])
            xlabel('time from pellet touch (ms)')
            ylabel('bits')
            title([tgLabel,' early; ',num2str(size(early_info_all.(tgLabel),1)),' pairs']);

            subplot(2,2,2+2*(trialGroupIdx-1)); hold on;
            plot(squeeze(mean(late_MI_all.(tgLabel).M1(:,:,26),1))')
            plot(squeeze(mean(late_MI_all.(tgLabel).DLS(:,:,26),1))')
            xlim([6 36])
            xticks([6 16 26 36])
            xticklabels([-1 -.5 0 .5])
            xlabel('time from pellet touch (ms)')
            ylabel('bits')
            title([tgLabel,' late; ',num2str(size(late_info_all.(tgLabel),1)),' pairs']); 

        end
        sgtitle([params.reachFeatures{fidx}, ' info profiles'])
            
        % Mean info
        figure('units','normalized','outerposition',[0 0 1 1]);
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            subplot(1,2,trialGroupIdx); hold on;
            bar(1,squeeze(mean(mean(early_MI_all.(tgLabel).M1(:,6:36,26),1),2))')
            errorbar(1,squeeze(std(mean(early_MI_all.(tgLabel).M1(:,6:36,26),2),[],1))','ko')
            bar(2,squeeze(mean(mean(early_MI_all.(tgLabel).DLS(:,6:36,26),1),2))')
            errorbar(2,squeeze(std(mean(early_MI_all.(tgLabel).DLS(:,6:36,26),2),[],1))','ko')
            bar(3,squeeze(mean(mean(late_MI_all.(tgLabel).M1(:,6:36,26),1),2))')
            errorbar(3,squeeze(std(mean(late_MI_all.(tgLabel).M1(:,6:36,26),2),[],1))','ko')
            bar(4,squeeze(mean(mean(late_MI_all.(tgLabel).DLS(:,6:36,26),1),2))')
            errorbar(4,squeeze(std(mean(late_MI_all.(tgLabel).DLS(:,6:36,26),2),[],1))','ko')
            ylim([0,0.01])

            xticks([1 2 3 4])
            xticklabels({'M1 naive','DLS naive','M1 skilled','DLS skilled'})
            title(tgLabel)
        end
        sgtitle(['Mean info ', params.reachFeatures{fidx}])
            
        figure('units','normalized','outerposition',[0 0 1 1])
        clearvars g
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};

            x = [1:39];
            y = [mean(early_info_all.(tgLabel)(:,:,27:51),3); mean(early_info_all.(tgLabel)(:,:,1:25),3)];
            c = [ones(1,size(mean(early_info_all.(tgLabel)(:,:,27:51),3),1)) 2*ones(1,size(mean(early_info_all.(tgLabel)(:,:,1:25),3),1))];
            g(trialGroupIdx,1)=gramm('x',x,'y',y,'color',c);
            g(trialGroupIdx,1).stat_summary('type','sem');
            g(trialGroupIdx,1).axe_property('XLim',[6 36]);
            g(trialGroupIdx,1).axe_property('YLim',ylim_info);    
            g(trialGroupIdx,1).set_title([tgLabel ' early']);        

            y = [mean(late_info_all.(tgLabel)(:,:,27:51),3); mean(late_info_all.(tgLabel)(:,:,1:25),3)];
            c = [ones(1,size(mean(late_info_all.(tgLabel)(:,:,27:51),3),1)) 2*ones(1,size(mean(late_info_all.(tgLabel)(:,:,1:25),3),1))];
            g(trialGroupIdx,2)=gramm('x',x,'y',y,'color',c);
            g(trialGroupIdx,2).stat_summary('type','sem');
            g(trialGroupIdx,2).axe_property('XLim',[6 36]);
            g(trialGroupIdx,2).axe_property('YLim',ylim_info);    
            g(trialGroupIdx,2).set_title([tgLabel ' late']);        

        end
        g.draw();

        if plot_params.save == 1
            print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
        end

        %%% PLOT Time-Delay max stat (MC July '23)
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            info_M1_DLS_early = squeeze(max(early_info_all.(tgLabel)(:,:,27:51),[],3));
            info_DLS_M1_early = squeeze(max(early_info_all.(tgLabel)(:,:,1:25),[],3));
            info_M1_DLS_late = squeeze(max(late_info_all.(tgLabel)(:,:,27:51),[],3));
            info_DLS_M1_late = squeeze(max(late_info_all.(tgLabel)(:,:,1:25),[],3));

            % take maximum over time per pair
%                     peaks_M1_DLS_early.(tgLabel) = max(info_M1_DLS_early(:,reach_bins),[],2);
            peaks_M1_DLS_early.(tgLabel) = max(info_M1_DLS_early(:,reach_bins),[],2);
            peaks_DLS_M1_early.(tgLabel) = max(info_DLS_M1_early(:,reach_bins),[],2);
            peaks_M1_DLS_late.(tgLabel) = max(info_M1_DLS_late(:,reach_bins),[],2);
            peaks_DLS_M1_late.(tgLabel) = max(info_DLS_M1_late(:,reach_bins),[],2);
            %sgtitle([params.reachFeatures{fidx},' '])
        end

        figure()
        hold on
        diff_early_succ = peaks_M1_DLS_early.success-peaks_DLS_M1_early.success;
        diff_early_fail = peaks_M1_DLS_early.failure-peaks_DLS_M1_early.failure;
        diff_late_succ = peaks_M1_DLS_late.success-peaks_DLS_M1_late.success;
        diff_late_fail = peaks_M1_DLS_late.failure-peaks_DLS_M1_late.failure;

        bar([1],[mean(diff_early_succ)],0.3,'b')
        bar([2],[mean(diff_late_succ)],0.3,'r')
        bar([4],[mean(diff_early_fail)],0.3,'b')
        bar([5],[mean(diff_late_fail)],0.3,'r')
        errorbar([1,2,4,5],[mean(diff_early_succ), mean(diff_late_succ), mean(diff_early_fail), mean(diff_late_fail)],...
            [std(diff_early_succ)/numel(diff_early_succ), std(diff_late_succ)/numel(diff_late_succ), std(diff_early_fail)/numel(diff_early_fail), std(diff_late_fail)/numel(diff_late_fail)],'o')

        [~,p_diff_succ] = ttest2(diff_early_succ,diff_late_succ);
        [~,p_naive_succ] = ttest(diff_early_succ);
        [~,p_skilled_succ] = ttest(diff_late_succ);
        [~,p_diff_fail] = ttest2(diff_early_fail,diff_late_fail);
        [~,p_naive_fail] = ttest(diff_early_fail);
        [~,p_skilled_fail] = ttest(diff_late_fail);
        ylabel('M1->DLS - DLS->M1 [bits]')
        title(['shared info | p_{diff}=' num2str(p_diff_succ,2) ' | p_{naive}=' num2str(p_naive_succ,2) ' | p_{skill}=' num2str(p_skilled_succ,2) ' \n' ...
            ' | p_{diff}=' num2str(p_diff_fail,2) ' | p_{naive}=' num2str(p_naive_fail,2) ' | p_{skill}=' num2str(p_skilled_fail,2)])
        ylim([-0.007,0.007])
%                 g.draw()
        yline(0,'k--')

        %% PLOT DELAYS
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            early_peakDelay.(tgLabel) = cell(1,length(early_pairs_per_animal.(tgLabel))-1);
            late_peakDelay.(tgLabel) = cell(1,length(early_pairs_per_animal.(tgLabel))-1);

            for animal = 1:length(early_pairs_per_animal.(tgLabel))-1

                tmp_peakDelay = [];
                for pair=1+sum(early_pairs_per_animal.(tgLabel)(1:animal)):sum(early_pairs_per_animal.(tgLabel)(1:animal+1))
                    [~,time_i] = max(squeeze(max(early_info_all.(tgLabel)(pair,:,:),[],3)));                    
                    if ismember(time_i,reach_bins)
                        [~, i] = max(squeeze(max(early_info_all.(tgLabel)(pair,reach_bins,:),[],2)));
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                early_peakDelay.(tgLabel){animal} = tmp_peakDelay;

                tmp_peakDelay = [];
                for pair = 1+sum(late_pairs_per_animal.(tgLabel)(1:animal)):sum(late_pairs_per_animal.(tgLabel)(1:animal+1))
                    [~,time_i] = max(squeeze(max(late_info_all.(tgLabel)(pair,:,:),[],3)));                    
                    if ismember(time_i,reach_bins)
                        [~, i] = max(squeeze(max(late_info_all.(tgLabel)(pair,reach_bins,:),[],2)));                        
                        tmp_peakDelay = [tmp_peakDelay i];
                    else
                        tmp_peakDelay = [tmp_peakDelay nan];
                    end
                end
                late_peakDelay.(tgLabel){animal} = tmp_peakDelay;

            end
        end
% 
        figure;
        subplot(2,2,1); hold on;
        delays_early = [early_peakDelay.success{:}];
        delays_late = [late_peakDelay.success{:}];
        histogram(delays_early(~isnan(delays_early)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
        histogram(delays_late(~isnan(delays_late)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
        [a,p]=ttest2(delays_early,delays_late);
        title(['success, p=',num2str(p,2)])

        subplot(2,2,2); hold on;
        tmp_e = histcounts([delays_early(~isnan(delays_early))],[1:5:51],'normalization','probability');
        tmp_l = histcounts([delays_late(~isnan(delays_late))],[1:5:51],'normalization','probability');
        scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
        plot(tmp_l-tmp_e,'color',[0 0 0])
        plot([0 11],[0 0],'color','k')
        plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
        xticks([1:1:10])     
        xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
        title('success')

        subplot(2,2,3); hold on;
        delays_early = [early_peakDelay.failure{:}];
        delays_late = [late_peakDelay.failure{:}];
        histogram(delays_early(~isnan(delays_early)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
        histogram(delays_late(~isnan(delays_late)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
        [a,p]=ttest2(delays_early,delays_late);
        title(['failure, p=',num2str(p,2)])

        subplot(2,2,4); hold on;
        tmp_e = histcounts([delays_early(~isnan(delays_early))],[1:5:51],'normalization','probability');
        tmp_l = histcounts([delays_late(~isnan(delays_late))],[1:5:51],'normalization','probability');
        scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
        plot(tmp_l-tmp_e,'color',[0 0 0])
        plot([0 11],[0 0],'color','k')
        plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
        xticks([1:1:10])     
        xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
        title('failure')

        if plot_params.save == 1
            print([plot_params.save_path '\PID_sharedInfo_delays_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
        end          

        figure;
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            change_pos.(tgLabel) = [];
            change_neg.(tgLabel) = [];
            for animal = 1:length(early_pairs_per_animal.(tgLabel))-1
                change_pos.(tgLabel) = [change_pos.(tgLabel) sum(late_peakDelay.(tgLabel){animal}>=27 & late_peakDelay.(tgLabel){animal}<=51)/length(late_peakDelay.(tgLabel){animal})-sum(early_peakDelay.(tgLabel){animal}>=27 & early_peakDelay.(tgLabel){animal}<=51)/length(early_peakDelay.(tgLabel){animal})];
                change_neg.(tgLabel) = [change_neg.(tgLabel) sum(late_peakDelay.(tgLabel){animal}>=1 & late_peakDelay.(tgLabel){animal}<=25)/length(late_peakDelay.(tgLabel){animal})-sum(early_peakDelay.(tgLabel){animal}>=1 & early_peakDelay.(tgLabel){animal}<=25)/length(early_peakDelay.(tgLabel){animal})];
            end         
            change_pos_all_animals.(tgLabel) = sum([late_peakDelay.(tgLabel){:}]>=27 & [late_peakDelay.(tgLabel){:}]<=51)/length(late_peakDelay.(tgLabel){animal})-sum([early_peakDelay.(tgLabel){:}]>=27 & [early_peakDelay.(tgLabel){:}]<=51)/length([early_peakDelay.(tgLabel){:}]);
            change_neg_all_animals.(tgLabel) = sum([late_peakDelay.(tgLabel){:}]>=1 & [late_peakDelay.(tgLabel){:}]<=25)/length(late_peakDelay.(tgLabel){animal})-sum([early_peakDelay.(tgLabel){:}]>=1 & [early_peakDelay.(tgLabel){:}]<=25)/length([early_peakDelay.(tgLabel){:}]);

            subplot(1,2,trialGroupIdx)
            hold on;
            for animal = 1:length(early_pairs_per_animal.(tgLabel))-1
                plot([1 2],[change_neg.(tgLabel)(animal) change_pos.(tgLabel)(animal)],'color',animal_colors(animal,:),'LineWidth',1)
            end
            errorbar(1,nanmean(change_neg.(tgLabel)),nanstd(change_neg.(tgLabel))/sqrt(sum(~isnan(change_neg.(tgLabel)))),'LineWidth',2,'color','k');            
            errorbar(2,nanmean(change_pos.(tgLabel)),nanstd(change_pos.(tgLabel))/sqrt(sum(~isnan(change_pos.(tgLabel)))),'LineWidth',2,'color','k');
            plot([.5 2.5],[0 0],'color','k','LineStyle','--')
            xlim([.5 2.5])
            [h p1, ~, stats] = ttest(change_pos.(tgLabel),change_neg.(tgLabel));
            disp(['all animal channel pair change: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                                

            title([tgLabel ' | ' num2str(p1)])

            subplot(2,2,trialGroupIdx+2)
            hold on;
            for animal = 1:length(early_pairs_per_animal.(tgLabel))-1
                plot([1 2],[change_neg.(tgLabel)(animal) change_pos.(tgLabel)(animal)],'color',animal_colors(animal,:),'LineWidth',1)
            end
            errorbar(1,nanmean(change_neg.(tgLabel)),nanstd(change_neg.(tgLabel))/sqrt(sum(~isnan(change_neg.(tgLabel)))),'LineWidth',2,'color','k');            
            errorbar(2,nanmean(change_pos.(tgLabel)),nanstd(change_pos.(tgLabel))/sqrt(sum(~isnan(change_pos.(tgLabel)))),'LineWidth',2,'color','k');
            plot([.5 2.5],[0 0],'color','k','LineStyle','--')
            xlim([.5 2.5])
            [h p1, ~, stats] = ttest(change_pos.(tgLabel),change_neg.(tgLabel));
            disp(['all animal channel pair change: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                                

            title([tgLabel ' | ' num2str(p1)])
        end
        sgtitle([params.reachFeatures{fidx}])
        if plot_params.save == 1
            print([plot_params.save_path '\PID_sharedInfo_delays_by_animal_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
        end    

        %%% PLOT Time-Delay max stat (MC July '23)
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            info_M1_DLS_early = squeeze(mean(early_info_all.(tgLabel)(:,:,27:51),3));
            info_DLS_M1_early = squeeze(mean(early_info_all.(tgLabel)(:,:,1:25),3));
            info_M1_DLS_late = squeeze(mean(late_info_all.(tgLabel)(:,:,27:51),3));
            info_DLS_M1_late = squeeze(mean(late_info_all.(tgLabel)(:,:,1:25),3));

            % take maximum over time per pair
            peaks_M1_DLS_early.(tgLabel) = mean(info_M1_DLS_early,2);
            peaks_DLS_M1_early.(tgLabel) = mean(info_DLS_M1_early,2);
            peaks_M1_DLS_late.(tgLabel) = mean(info_M1_DLS_late,2);
            peaks_DLS_M1_late.(tgLabel) = mean(info_DLS_M1_late,2);
            %sgtitle([params.reachFeatures{fidx},' '])
        end
        
         %% HISTOGRAM
         
        figure();
         
        for trialGroupIdx = 1:length(trial_group_names)
            tgLabel = trial_group_names{trialGroupIdx};
            
            max_early_DLStoM1 = [];
            for pair = 1:size(early_info_all.(tgLabel),1)
                max_early_DLStoM1 = [max_early_DLStoM1 max(squeeze(early_info_all.(tgLabel)(pair,reach_bins,1:25)),[],'all')];
            end
            max_early_M1toDLS = [];
            for pair = 1:size(early_info_all.(tgLabel),1)
                max_early_M1toDLS = [max_early_M1toDLS max(squeeze(early_info_all.(tgLabel)(pair,reach_bins,27:51)),[],'all')];
            end
            max_late_DLStoM1 = [];
            for pair = 1:size(late_info_all.(tgLabel),1)
                max_late_DLStoM1 = [max_late_DLStoM1 max(squeeze(late_info_all.(tgLabel)(pair,reach_bins,1:25)),[],'all')];
            end
            max_late_M1toDLS = [];
            for pair = 1:size(late_info_all.(tgLabel),1)
                max_late_M1toDLS = [max_late_M1toDLS max(squeeze(late_info_all.(tgLabel)(pair,reach_bins,27:51)),[],'all')];
            end

            subplot(2,4,1+4*(trialGroupIdx-1)); hold on;
            histogram(max_early_M1toDLS,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
            histogram(max_early_DLStoM1,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
            [p,~,stats] = signrank(max_early_M1toDLS,max_early_DLStoM1);
            title([num2str(p) '' tgLabel]);
            xlim([0 0.1]);
            
            subplot(2,4,3+4*(trialGroupIdx-1));
            hold on;
            cdfplot(max_early_M1toDLS);
            cdfplot(max_early_DLStoM1);
            title(['Naive ' tgLabel])
            xlim([0 0.1]);
            grid off;
            
            subplot(2,4,2+4*(trialGroupIdx-1)); hold on;
            histogram(max_late_M1toDLS,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
            histogram(max_late_DLStoM1,[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
            [p,~,stats] = signrank(max_late_M1toDLS,max_late_DLStoM1);
            title([num2str(p) '' tgLabel]);
            xlim([0 0.1]);
            
            subplot(2,4,4+4*(trialGroupIdx-1));
            hold on;
            cdfplot(max_late_M1toDLS);
            cdfplot(max_late_DLStoM1);
            xlim([0 0.1]);
            title(['Skilled ' tgLabel])
            grid off;
        end
        sgtitle(params.reachFeatures{fidx})
                
        if plot_params.save == 1
            print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '_early_late_Histogram.svg'],'-painters','-dsvg');
            print([plot_params.save_path '\PID_sharedInfo_' params.reachFeatures{fidx} '_early_late_Histogram.png'],'-dpng');
        end  
    end
end

%%% Plot pooled features
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};
    early_info_all.(tgLabel) = early_info_all_pooled.(tgLabel);
    early_MI_all.(tgLabel).M1 = early_MI_all_pooled.(tgLabel).M1;
    early_MI_all.(tgLabel).DLS = early_MI_all_pooled.(tgLabel).DLS;
    early_pairs_per_animal.(tgLabel)(animal) = early_pairs_per_animal_pooled.(tgLabel)(animal);

    late_info_all.(tgLabel) = late_info_all_pooled.(tgLabel);
    late_MI_all.(tgLabel).M1 = late_MI_all_pooled.(tgLabel).M1;
    late_MI_all.(tgLabel).DLS = late_MI_all_pooled.(tgLabel).DLS;
    late_pairs_per_animal.(tgLabel)(animal) = late_pairs_per_animal_pooled.(tgLabel)(animal);
end
%% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)
disp('Plot spatiotemporal maps')
figure('units','normalized','outerposition',[0 0 1 1]);
y_lim_info_early = [0, 0.015];
y_lim_info_late = [0, 0.025];

for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};

    % Early pairs of animal T107 already removed
    
    subplot(2,2,(trialGroupIdx*2)-1); hold on;
    imagesc(squeeze(mean(early_info_all.(tgLabel),1))',y_lim_info_early)
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
    title([tgLabel,' early; ',num2str(size(early_info_all.(tgLabel),1)),' pairs']);

    subplot(2,2,trialGroupIdx*2); hold on;
    imagesc(squeeze(mean(late_info_all.(tgLabel),1))',y_lim_info_late)
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
    title([tgLabel,' late; ',num2str(size(late_info_all.(tgLabel),1)),' pairs']); 
end
sgtitle(['pooled time-delay maps'])

% Info time profiles
figure('units','normalized','outerposition',[0 0 1 1]);
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};

    subplot(2,2,1+2*(trialGroupIdx-1)); hold on;
    plot(squeeze(mean(early_MI_all.(tgLabel).M1(:,:,26),1))')
    plot(squeeze(mean(early_MI_all.(tgLabel).DLS(:,:,26),1))')
    xlim([6 36])
    xticks([6 16 26 36])
    xticklabels([-1 -.5 0 .5])
    xlabel('time from pellet touch (ms)')
    ylabel('bits')
    title([tgLabel,' early; ',num2str(size(early_info_all.(tgLabel),1)),' pairs']);

    subplot(2,2,2+2*(trialGroupIdx-1)); hold on;
    plot(squeeze(mean(late_MI_all.(tgLabel).M1(:,:,26),1))')
    plot(squeeze(mean(late_MI_all.(tgLabel).DLS(:,:,26),1))')
    xlim([6 36])
    xticks([6 16 26 36])
    xticklabels([-1 -.5 0 .5])
    xlabel('time from pellet touch (ms)')
    ylabel('bits')
    title([tgLabel,' late; ',num2str(size(late_info_all.(tgLabel),1)),' pairs']); 

end
sgtitle(['pooled info profiles'])

% Mean info
figure('units','normalized','outerposition',[0 0 1 1]);
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};
    subplot(1,2,trialGroupIdx); hold on;
    bar(1,squeeze(mean(mean(early_MI_all.(tgLabel).M1(:,6:36,26),1),2))')
    errorbar(1,squeeze(std(mean(early_MI_all.(tgLabel).M1(:,6:36,26),2),[],1))','ko')
    bar(2,squeeze(mean(mean(early_MI_all.(tgLabel).DLS(:,6:36,26),1),2))')
    errorbar(2,squeeze(std(mean(early_MI_all.(tgLabel).DLS(:,6:36,26),2),[],1))','ko')
    bar(3,squeeze(mean(mean(late_MI_all.(tgLabel).M1(:,6:36,26),1),2))')
    errorbar(3,squeeze(std(mean(late_MI_all.(tgLabel).M1(:,6:36,26),2),[],1))','ko')
    bar(4,squeeze(mean(mean(late_MI_all.(tgLabel).DLS(:,6:36,26),1),2))')
    errorbar(4,squeeze(std(mean(late_MI_all.(tgLabel).DLS(:,6:36,26),2),[],1))','ko')
    ylim([0,0.01])

    xticks([1 2 3 4])
    xticklabels({'M1 naive','DLS naive','M1 skilled','DLS skilled'})
    title(tgLabel)
end
sgtitle(['Mean info pooled'])

figure('units','normalized','outerposition',[0 0 1 1])
clearvars g
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};

    x = [1:39];
    y = [mean(early_info_all.(tgLabel)(:,:,27:51),3); mean(early_info_all.(tgLabel)(:,:,1:25),3)];
    c = [ones(1,size(mean(early_info_all.(tgLabel)(:,:,27:51),3),1)) 2*ones(1,size(mean(early_info_all.(tgLabel)(:,:,1:25),3),1))];
    g(trialGroupIdx,1)=gramm('x',x,'y',y,'color',c);
    g(trialGroupIdx,1).stat_summary('type','sem');
    g(trialGroupIdx,1).axe_property('XLim',[6 36]);
    g(trialGroupIdx,1).axe_property('YLim',ylim_info);    
    g(trialGroupIdx,1).set_title([tgLabel ' early']);        

    y = [mean(late_info_all.(tgLabel)(:,:,27:51),3); mean(late_info_all.(tgLabel)(:,:,1:25),3)];
    c = [ones(1,size(mean(late_info_all.(tgLabel)(:,:,27:51),3),1)) 2*ones(1,size(mean(late_info_all.(tgLabel)(:,:,1:25),3),1))];
    g(trialGroupIdx,2)=gramm('x',x,'y',y,'color',c);
    g(trialGroupIdx,2).stat_summary('type','sem');
    g(trialGroupIdx,2).axe_property('XLim',[6 36]);
    g(trialGroupIdx,2).axe_property('YLim',ylim_info);    
    g(trialGroupIdx,2).set_title([tgLabel ' late']);        

end
g.draw();

if plot_params.save == 1
    print([plot_params.save_path '\PID_sharedInfo_pooled.eps'],'-painters','-depsc');
end

%%% PLOT Time-Delay max stat (MC July '23)
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};
%     info_M1_DLS_early = squeeze(max(early_info_all.(tgLabel)(:,:,27:51),[],3));
%     info_DLS_M1_early = squeeze(max(early_info_all.(tgLabel)(:,:,1:25),[],3));
%     info_M1_DLS_late = squeeze(max(late_info_all.(tgLabel)(:,:,27:51),[],3));
%     info_DLS_M1_late = squeeze(max(late_info_all.(tgLabel)(:,:,1:25),[],3));
% 
%     % take maximum over time per pair
% %                     peaks_M1_DLS_early.(tgLabel) = max(info_M1_DLS_early(:,reach_bins),[],2);
%     peaks_M1_DLS_early.(tgLabel) = max(info_M1_DLS_early(:,reach_bins),[],2);
%     peaks_DLS_M1_early.(tgLabel) = max(info_DLS_M1_early(:,reach_bins),[],2);
%     peaks_M1_DLS_late.(tgLabel) = max(info_M1_DLS_late(:,reach_bins),[],2);
%     peaks_DLS_M1_late.(tgLabel) = max(info_DLS_M1_late(:,reach_bins),[],2);

    peaks_M1_DLS_early.(tgLabel) = squeeze(max(early_info_all.(tgLabel)(:,reach_bins,27:51),[],[2 3]));
    peaks_DLS_M1_early.(tgLabel) = squeeze(max(early_info_all.(tgLabel)(:,reach_bins,1:25),[],[2 3]));
    peaks_M1_DLS_late.(tgLabel) = squeeze(max(late_info_all.(tgLabel)(:,reach_bins,27:51),[],[2 3]));
    peaks_DLS_M1_late.(tgLabel) = squeeze(max(late_info_all.(tgLabel)(:,reach_bins,1:25),[],[2 3]));
    %sgtitle([params.reachFeatures{fidx},' '])
end

figure()
hold on
diff_early_succ = peaks_M1_DLS_early.success-peaks_DLS_M1_early.success;
diff_early_fail = peaks_M1_DLS_early.failure-peaks_DLS_M1_early.failure;
diff_late_succ = peaks_M1_DLS_late.success-peaks_DLS_M1_late.success;
diff_late_fail = peaks_M1_DLS_late.failure-peaks_DLS_M1_late.failure;

bar([1],[mean(diff_early_succ)],0.3,'b')
bar([2],[mean(diff_late_succ)],0.3,'r')
bar([4],[mean(diff_early_fail)],0.3,'b')
bar([5],[mean(diff_late_fail)],0.3,'r')
errorbar([1,2,4,5],[mean(diff_early_succ), mean(diff_late_succ), mean(diff_early_fail), mean(diff_late_fail)],...
    [std(diff_early_succ)/numel(diff_early_succ), std(diff_late_succ)/numel(diff_late_succ), std(diff_early_fail)/numel(diff_early_fail), std(diff_late_fail)/numel(diff_late_fail)],'o')

[~,p_diff_succ] = ttest2(diff_early_succ,diff_late_succ);
[~,p_naive_succ] = ttest(diff_early_succ);
[~,p_skilled_succ] = ttest(diff_late_succ);
[~,p_diff_fail] = ttest2(diff_early_fail,diff_late_fail);
[~,p_naive_fail] = ttest(diff_early_fail);
[~,p_skilled_fail] = ttest(diff_late_fail);
ylabel('M1->DLS - DLS->M1 [bits]')
title(['shared info | p_{diff}=' num2str(p_diff_succ,2) ' | p_{naive}=' num2str(p_naive_succ,2) ' | p_{skill}=' num2str(p_skilled_succ,2) ' \n' ...
    ' | p_{diff}=' num2str(p_diff_fail,2) ' | p_{naive}=' num2str(p_naive_fail,2) ' | p_{skill}=' num2str(p_skilled_fail,2)])
ylim([-0.007,0.007])
%                 g.draw()
yline(0,'k--')

%% PLOT DELAYS
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};
    early_peakDelay.(tgLabel) = cell(1,length(early_pairs_per_animal_pooled.(tgLabel))-1);
    late_peakDelay.(tgLabel) = cell(1,length(early_pairs_per_animal_pooled.(tgLabel))-1);

    for animal = 1:length(early_pairs_per_animal_pooled.(tgLabel))-1

        tmp_peakDelay = [];
        for pair=1+sum(early_pairs_per_animal_pooled.(tgLabel)(1:animal)):sum(early_pairs_per_animal_pooled.(tgLabel)(1:animal+1))
            [~,time_i] = max(squeeze(max(early_info_all.(tgLabel)(pair,:,:),[],3)));                    
            if ismember(time_i,reach_bins)
                [~, i] = max(squeeze(max(early_info_all.(tgLabel)(pair,reach_bins,:),[],2)));
                tmp_peakDelay = [tmp_peakDelay i];
            else
                tmp_peakDelay = [tmp_peakDelay nan];
%                 [~, i] = max(squeeze(max(early_info_all.(tgLabel)(pair,reach_bins,:),[],2)));
%                 tmp_peakDelay = [tmp_peakDelay i];
            end
        end
        early_peakDelay.(tgLabel){animal} = tmp_peakDelay;

        tmp_peakDelay = [];
        for pair = 1+sum(late_pairs_per_animal_pooled.(tgLabel)(1:animal)):sum(late_pairs_per_animal_pooled.(tgLabel)(1:animal+1))
            [~,time_i] = max(squeeze(max(late_info_all.(tgLabel)(pair,:,:),[],3)));                    
            if ismember(time_i,reach_bins)
                [~, i] = max(squeeze(max(late_info_all.(tgLabel)(pair,reach_bins,:),[],2)));
                tmp_peakDelay = [tmp_peakDelay i];
            else
                tmp_peakDelay = [tmp_peakDelay nan];
% %                 [~, i] = max(squeeze(max(late_info_all.(tgLabel)(pair,reach_bins,:),[],2)));
% %                 tmp_peakDelay = [tmp_peakDelay i];
            end
        end
        late_peakDelay.(tgLabel){animal} = tmp_peakDelay;

    end
end
% 
figure;
subplot(3,2,1); hold on;
delays_early = [cell2mat(early_peakDelay.success)];
delays_late = [cell2mat(late_peakDelay.success)];
histogram(delays_early(~isnan(delays_early)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
histogram(delays_late(~isnan(delays_late)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
[a,p]=ttest2(delays_early,delays_late);
title(['success, p=',num2str(p,2)])

subplot(3,2,2); hold on;
tmp_e = histcounts([delays_early(~isnan(delays_early))],[1:5:51],'normalization','probability');
tmp_l = histcounts([delays_late(~isnan(delays_late))],[1:5:51],'normalization','probability');
scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
plot(tmp_l-tmp_e,'color',[0 0 0])
plot([0 11],[0 0],'color','k')
plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
xticks([1:1:10])     
xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
title('success')

subplot(3,2,3); hold on;
delays_early = [cell2mat(early_peakDelay.failure)];
delays_late = [cell2mat(late_peakDelay.failure)];
histogram(delays_early(~isnan(delays_early)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
histogram(delays_late(~isnan(delays_late)),[1:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
[a,p]=ttest2(delays_early,delays_late);
title(['failure, p=',num2str(p,2)])

subplot(3,2,4); hold on;
tmp_e = histcounts([delays_early(~isnan(delays_early))],[1:5:51],'normalization','probability');
tmp_l = histcounts([delays_late(~isnan(delays_late))],[1:5:51],'normalization','probability');
scatter(1:10,tmp_l-tmp_e,30,[0 0 0],'filled')
plot(tmp_l-tmp_e,'color',[0 0 0])
plot([0 11],[0 0],'color','k')
plot([5.5 5.5],[-.2 .2],'color','r','LineStyle','--')
xticks([1:1:10])     
xticklabels({'-250 to -200','-200 to -150','-150 to -100','-100 to -50','-50 to 0','0 to 50','50 to 100','100 to 150','150 to 200','200 to 250'})
title('failure')

subplot(3,2,5); hold on;
delays_early = ([cell2mat(early_peakDelay.success)]-[cell2mat(early_peakDelay.failure)]);
delays_late = ([cell2mat(late_peakDelay.success)]-[cell2mat(late_peakDelay.failure)]);
histogram(delays_early(~isnan(delays_early)),[-51:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','b')
histogram(delays_late(~isnan(delays_late)),[-51:5:51],'normalization','probability','DisplayStyle','Stairs','EdgeColor','r')
[a,p]=ttest2(delays_early,delays_late);
title(['succ-fail, p=',num2str(p,2)])

subplot(3,2,6); hold on;
tmp_e = histcounts([delays_early(~isnan(delays_early))],[-51:5:51],'normalization','probability');
tmp_l = histcounts([delays_late(~isnan(delays_late))],[-51:5:51],'normalization','probability');
scatter(1:20,tmp_l-tmp_e,30,[0 0 0],'filled')
plot(tmp_l-tmp_e,'color',[0 0 0])
plot([0 21],[0 0],'color','k')
plot([10.5 10.5],[-.2 .2],'color','r','LineStyle','--')
xticks([1:2:20])     
xticklabels({'-500 to -400','-400 to -300','-300 to -200','-200 to -100','-100 to 0','0 to 100','100 to 200','200 to 300','300 to 400','400 to 500'})
title('succ-fail')

if plot_params.save == 1
    print([plot_params.save_path '\PID_sharedInfo_delays_pooled.eps'],'-painters','-depsc');
end          


early_delays_succ = ([early_peakDelay.success{:}]-26)*10;
late_delays_succ = ([late_peakDelay.success{:}]-26)*10;
early_delays_fail = ([early_peakDelay.failure{:}]-26)*10;
late_delays_fail = ([late_peakDelay.failure{:}]-26)*10;
tmp_e_succ = histcounts([early_delays_succ(~isnan(early_delays_succ))],[-25:5:25],'normalization','probability');
tmp_l_succ = histcounts([late_delays_succ(~isnan(late_delays_succ))],[-25:5:25],'normalization','probability');
tmp_e_fail = histcounts([early_delays_fail(~isnan(early_delays_fail))],[-25:5:25],'normalization','probability');
tmp_l_fail = histcounts([late_delays_fail(~isnan(late_delays_fail))],[-25:5:25],'normalization','probability');

figure;
lowPrc = 25;
highPrc = 75;
dpt = 10;
hold on
meanPlot1 = mean(early_delays_succ,'omitnan');
h(1)=bar([1],meanPlot1,'b');
meanPlot2 = mean(late_delays_succ,'omitnan');
h(2)=bar([2],meanPlot2,'r');
error1H = std(early_delays_succ,'omitnan')/sqrt(sum(~isnan(early_delays_succ)));
error1L = std(early_delays_succ,'omitnan')/sqrt(sum(~isnan(early_delays_succ)));
error2H = std(late_delays_succ,'omitnan')/sqrt(sum(~isnan(late_delays_succ)));
error2L = std(late_delays_succ,'omitnan')/sqrt(sum(~isnan(late_delays_succ)));
% error1H = prctile(early_delays_succ,highPrc)-meanPlot1;
% error1L = meanPlot1-prctile(early_delays_succ,lowPrc);
% error2H = prctile(late_delays_succ,highPrc)-meanPlot2;
% error2L = meanPlot2-prctile(late_delays_succ,lowPrc);

errorbar([1,2], [meanPlot1,meanPlot2],[error1L,error2L],[error1H,error2H],'ko');
[a,pval(1)]=ttest2(early_delays_succ,late_delays_succ);
maxY = max([nanmean(early_delays_succ),nanmean(late_delays_succ)]);
pvalues_plot(pval(1),1.5,1.4*maxY,dpt,1,14,0.5,'k',0)
[a,p1]=ttest(early_delays_succ);
[a,p2]=ttest(late_delays_succ);
pvalues_plot(p1,1,1.2*maxY,dpt,1,11,0,'k',0)
pvalues_plot(p2,2,1.2*maxY,dpt,1,11,0,'k',0)
legend([h(1) h(2)],'naive','skilled','AutoUpdate','off')

meanPlot1 = mean(early_delays_fail,'omitnan');
h(1)=bar([4],meanPlot1,'b');
meanPlot2 = mean(late_delays_fail,'omitnan');
h(2)=bar([5],meanPlot2,'r');
error1H = std(early_delays_fail,'omitnan')/sqrt(sum(~isnan(early_delays_fail)));
error1L = std(early_delays_fail,'omitnan')/sqrt(sum(~isnan(early_delays_fail)));
error2H = std(late_delays_fail,'omitnan')/sqrt(sum(~isnan(late_delays_fail)));
error2L = std(late_delays_fail,'omitnan')/sqrt(sum(~isnan(late_delays_fail)));
% error1H = prctile(early_delays_fail,highPrc)-meanPlot1;
% error1L = meanPlot1-prctile(early_delays_fail,lowPrc);
% error2H = prctile(late_delays_fail,highPrc)-meanPlot2;
% error2L = meanPlot2-prctile(late_delays_fail,lowPrc);
errorbar([4,5], [meanPlot1,meanPlot2],[error1L,error2L],[error1H,error2H],'ko');

[a,pval(2)]=ttest2(early_delays_fail,late_delays_fail);
maxY = max([nanmean(early_delays_fail),nanmean(late_delays_fail)]);
pvalues_plot(pval(2),4.5,1.5*maxY,dpt,1,14,0.5,'k',0)
[a,p1]=ttest(early_delays_fail);
[a,p2]=ttest(late_delays_fail);
pvalues_plot(p1,4,1.2*maxY,dpt,1,11,0,'k',0)
pvalues_plot(p2,5,1.2*maxY,dpt,1,11,0,'k',0)

early_delays_diff = early_delays_succ-early_delays_fail;
late_delays_diff = late_delays_succ-late_delays_fail;

meanPlot1 = mean(early_delays_diff,'omitnan');
h(1)=bar([7],meanPlot1,'b');
meanPlot2 = mean(late_delays_diff,'omitnan');
h(2)=bar([8],meanPlot2,'r');
error1H = std(early_delays_diff,'omitnan')/sqrt(sum(~isnan(early_delays_diff)));
error1L = std(early_delays_diff,'omitnan')/sqrt(sum(~isnan(early_delays_diff)));
error2H = std(late_delays_diff,'omitnan')/sqrt(sum(~isnan(late_delays_diff)));
error2L = std(late_delays_diff,'omitnan')/sqrt(sum(~isnan(late_delays_diff)));
% error1H = prctile(early_delays_diff,highPrc)-meanPlot1;
% error1L = meanPlot1-prctile(early_delays_diff,lowPrc);
% error2H = prctile(late_delays_diff,highPrc)-meanPlot2;
% error2L = meanPlot2-prctile(late_delays_diff,lowPrc);
errorbar([7,8], [meanPlot1,meanPlot2],[error1L,error2L],[error1H,error2H],'ko');

% bar([7], ([nanmean(early_delays_diff)]),'b')
% bar([8], ([nanmean(late_delays_diff)]),'r')
% errorbar([7,8], [nanmean(early_delays_diff),nanmean(late_delays_diff)],[std(early_delays_diff,'omitnan')/sqrt(sum(~isnan(early_delays_diff))),std(late_delays_diff,'omitnan')/sqrt(sum(~isnan(late_delays_diff)))],'ko');
[a,pval(3)]=ttest2(early_delays_diff,late_delays_diff);
maxY = max([nanmean(early_delays_diff),nanmean(late_delays_diff)]);
pvalues_plot(pval(3),7.5,1.6*maxY,dpt,1,14,0.5,'k',0)
[a,p1]=ttest(early_delays_diff);
[a,p2]=ttest(late_delays_diff);
pvalues_plot(p1,7,1.2*maxY,dpt,1,11,0,'k',0)
pvalues_plot(p2,8,1.2*maxY,dpt,1,11,0,'k',0)

ylabel('delay [ms]')
set(gca, 'XTick',[1.5,4.5,7.5], 'XTickLabel',{'Success','Failure','Succ-Fail'})
title('Delays statistics within success and failure trials')
ylim([-70,70])

figure;
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};
    change_pos.(tgLabel) = [];
    change_neg.(tgLabel) = [];
    for animal = 1:length(early_pairs_per_animal.(tgLabel))-1
        change_pos.(tgLabel) = [change_pos.(tgLabel) sum(late_peakDelay.(tgLabel){animal}>=27 & late_peakDelay.(tgLabel){animal}<=51)/length(late_peakDelay.(tgLabel){animal})-sum(early_peakDelay.(tgLabel){animal}>=27 & early_peakDelay.(tgLabel){animal}<=51)/length(early_peakDelay.(tgLabel){animal})];
        change_neg.(tgLabel) = [change_neg.(tgLabel) sum(late_peakDelay.(tgLabel){animal}>=1 & late_peakDelay.(tgLabel){animal}<=25)/length(late_peakDelay.(tgLabel){animal})-sum(early_peakDelay.(tgLabel){animal}>=1 & early_peakDelay.(tgLabel){animal}<=25)/length(early_peakDelay.(tgLabel){animal})];
    end         
    change_pos_all_animals.(tgLabel) = sum([late_peakDelay.(tgLabel){:}]>=27 & [late_peakDelay.(tgLabel){:}]<=51)/length(late_peakDelay.(tgLabel){animal})-sum([early_peakDelay.(tgLabel){:}]>=27 & [early_peakDelay.(tgLabel){:}]<=51)/length([early_peakDelay.(tgLabel){:}]);
    change_neg_all_animals.(tgLabel) = sum([late_peakDelay.(tgLabel){:}]>=1 & [late_peakDelay.(tgLabel){:}]<=25)/length(late_peakDelay.(tgLabel){animal})-sum([early_peakDelay.(tgLabel){:}]>=1 & [early_peakDelay.(tgLabel){:}]<=25)/length([early_peakDelay.(tgLabel){:}]);

    subplot(1,2,trialGroupIdx)
    hold on;
    for animal = 1:length(early_pairs_per_animal.(tgLabel))-1
        plot([1 2],[change_neg.(tgLabel)(animal) change_pos.(tgLabel)(animal)],'color',animal_colors(animal,:),'LineWidth',1)
    end
    errorbar(1,nanmean(change_neg.(tgLabel)),nanstd(change_neg.(tgLabel))/sqrt(sum(~isnan(change_neg.(tgLabel)))),'LineWidth',2,'color','k');            
    errorbar(2,nanmean(change_pos.(tgLabel)),nanstd(change_pos.(tgLabel))/sqrt(sum(~isnan(change_pos.(tgLabel)))),'LineWidth',2,'color','k');
    plot([.5 2.5],[0 0],'color','k','LineStyle','--')
    xlim([.5 2.5])
    [h p1, ~, stats] = ttest(change_pos.(tgLabel),change_neg.(tgLabel));
    disp(['all animal channel pair change: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                                

    title([tgLabel ' | ' num2str(p1)])

    subplot(2,2,trialGroupIdx+2)
    hold on;
    for animal = 1:length(early_pairs_per_animal.(tgLabel))-1
        plot([1 2],[change_neg.(tgLabel)(animal) change_pos.(tgLabel)(animal)],'color',animal_colors(animal,:),'LineWidth',1)
    end
    errorbar(1,nanmean(change_neg.(tgLabel)),nanstd(change_neg.(tgLabel))/sqrt(sum(~isnan(change_neg.(tgLabel)))),'LineWidth',2,'color','k');            
    errorbar(2,nanmean(change_pos.(tgLabel)),nanstd(change_pos.(tgLabel))/sqrt(sum(~isnan(change_pos.(tgLabel)))),'LineWidth',2,'color','k');
    plot([.5 2.5],[0 0],'color','k','LineStyle','--')
    xlim([.5 2.5])
    [h p1, ~, stats] = ttest(change_pos.(tgLabel),change_neg.(tgLabel));
    disp(['all animal channel pair change: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p1)]);                                                

    title([tgLabel ' | ' num2str(p1)])
end
sgtitle(['pooled'])
if plot_params.save == 1
    print([plot_params.save_path '\PID_sharedInfo_delays_by_animal_pooled.eps'],'-painters','-depsc');
end    

%%% PLOT Time-Delay max stat (MC July '23)
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};
    info_M1_DLS_early = squeeze(mean(early_info_all.(tgLabel)(:,:,27:51),3));
    info_DLS_M1_early = squeeze(mean(early_info_all.(tgLabel)(:,:,1:25),3));
    info_M1_DLS_late = squeeze(mean(late_info_all.(tgLabel)(:,:,27:51),3));
    info_DLS_M1_late = squeeze(mean(late_info_all.(tgLabel)(:,:,1:25),3));

    % take maximum over time per pair
    peaks_M1_DLS_early.(tgLabel) = mean(info_M1_DLS_early,2);
    peaks_DLS_M1_early.(tgLabel) = mean(info_DLS_M1_early,2);
    peaks_M1_DLS_late.(tgLabel) = mean(info_M1_DLS_late,2);
    peaks_DLS_M1_late.(tgLabel) = mean(info_DLS_M1_late,2);
    %sgtitle([params.reachFeatures{fidx},' '])
end

%% HISTOGRAM

figure();
clear max_early_DLStoM1 max_early_M1toDLS max_late_DLStoM1 max_late_M1toDLS
for trialGroupIdx = 1:length(trial_group_names)
    tgLabel = trial_group_names{trialGroupIdx};

    max_early_DLStoM1.(tgLabel) = [];
    for pair = 1:size(early_info_all.(tgLabel),1)
        max_early_DLStoM1.(tgLabel) = [max_early_DLStoM1.(tgLabel) max(squeeze(early_info_all.(tgLabel)(pair,reach_bins,1:25)),[],'all')];
    end
    max_early_M1toDLS.(tgLabel) = [];
    for pair = 1:size(early_info_all.(tgLabel),1)
        max_early_M1toDLS.(tgLabel) = [max_early_M1toDLS.(tgLabel) max(squeeze(early_info_all.(tgLabel)(pair,reach_bins,27:51)),[],'all')];
    end
    max_late_DLStoM1.(tgLabel) = [];
    for pair = 1:size(late_info_all.(tgLabel),1)
        max_late_DLStoM1.(tgLabel) = [max_late_DLStoM1.(tgLabel) max(squeeze(late_info_all.(tgLabel)(pair,reach_bins,1:25)),[],'all')];
    end
    max_late_M1toDLS.(tgLabel) = [];
    for pair = 1:size(late_info_all.(tgLabel),1)
        max_late_M1toDLS.(tgLabel) = [max_late_M1toDLS.(tgLabel) max(squeeze(late_info_all.(tgLabel)(pair,reach_bins,27:51)),[],'all')];
    end

    subplot(2,4,1+4*(trialGroupIdx-1)); hold on;
    histogram(max_early_M1toDLS.(tgLabel),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
    histogram(max_early_DLStoM1.(tgLabel),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
    [p,~,stats] = signrank(max_early_M1toDLS.(tgLabel),max_early_DLStoM1.(tgLabel));
    title([num2str(p) '' tgLabel]);
    xlim([0 0.1]);

    subplot(2,4,3+4*(trialGroupIdx-1));
    hold on;
    cdfplot(max_early_M1toDLS.(tgLabel));
    cdfplot(max_early_DLStoM1.(tgLabel));
    title(['Naive ' tgLabel])
    xlim([0 0.1]);
    grid off;

    subplot(2,4,2+4*(trialGroupIdx-1)); hold on;
    histogram(max_late_M1toDLS.(tgLabel),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');    
    histogram(max_late_DLStoM1.(tgLabel),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
    [p,~,stats] = signrank(max_late_M1toDLS.(tgLabel),max_late_DLStoM1.(tgLabel));
    title([num2str(p) '' tgLabel]);
    xlim([0 0.1]);

    subplot(2,4,4+4*(trialGroupIdx-1));
    hold on;
    cdfplot(max_late_M1toDLS.(tgLabel));
    cdfplot(max_late_DLStoM1.(tgLabel));
    xlim([0 0.1]);
    title(['Skilled ' tgLabel])
    grid off;
end
sgtitle('pooled')

max_early_diff.success = (max_early_DLStoM1.success-max_early_M1toDLS.success);
max_late_diff.success = (max_late_DLStoM1.success-max_late_M1toDLS.success);
max_early_diff.failure = (max_early_DLStoM1.failure-max_early_M1toDLS.failure);
max_late_diff.failure = (max_late_DLStoM1.failure-max_late_M1toDLS.failure);

% Delays reversal
figure;
lowPrc = 25;
highPrc = 75;
dpt = 0.001;
hold on
% success
meanPlot1 = mean(max_early_diff.success,'omitnan');
h(1)=bar([1],meanPlot1,'b');
meanPlot2 = mean(max_late_diff.success,'omitnan');
h(2)=bar([2],meanPlot2,'r');
error1H = std(max_early_diff.success,'omitnan')/sqrt(sum(~isnan(max_early_diff.success)));
error1L = std(max_early_diff.success,'omitnan')/sqrt(sum(~isnan(max_early_diff.success)));
error2H = std(max_late_diff.success,'omitnan')/sqrt(sum(~isnan(max_late_diff.success)));
error2L = std(max_late_diff.success,'omitnan')/sqrt(sum(~isnan(max_late_diff.success)));

errorbar([1,2], [meanPlot1,meanPlot2],[error1L,error2L],[error1H,error2H],'ko');
[pval(1)]=ranksum(max_early_diff.success,max_late_diff.success);
maxY = max([nanmean(max_early_diff.success),nanmean(max_late_diff.success)]);
pvalues_plot(pval(1),1.5,1.3*maxY,dpt,1,14,0.5,'k',0)
[p1]=signrank(max_early_diff.success);
[p2]=signrank(max_late_diff.success);
pvalues_plot(p1,1,1.1*maxY,dpt,1,11,0,'k',0)
pvalues_plot(p2,2,1.1*maxY,dpt,1,11,0,'k',0)
legend([h(1) h(2)],'naive','skilled','AutoUpdate','off')

% failure
meanPlot1 = mean(max_early_diff.failure,'omitnan');
h(1)=bar([4],meanPlot1,'b');
meanPlot2 = mean(max_late_diff.failure,'omitnan');
h(2)=bar([5],meanPlot2,'r');
error1H = std(max_early_diff.failure,'omitnan')/sqrt(sum(~isnan(max_early_diff.failure)));
error1L = std(max_early_diff.failure,'omitnan')/sqrt(sum(~isnan(max_early_diff.failure)));
error2H = std(max_late_diff.failure,'omitnan')/sqrt(sum(~isnan(max_late_diff.failure)));
error2L = std(max_late_diff.failure,'omitnan')/sqrt(sum(~isnan(max_late_diff.failure)));
errorbar([4,5], [meanPlot1,meanPlot2],[error1L,error2L],[error1H,error2H],'ko');

[pval(2)]=ranksum(max_early_diff.failure,max_late_diff.failure);
maxY = max([nanmean(max_early_diff.failure),nanmean(max_late_diff.failure)]);
pvalues_plot(pval(2),4.5,2*maxY,dpt,1,14,0.5,'k',0)
[p1]=signrank(max_early_diff.failure);
[p2]=signrank(max_late_diff.failure);
pvalues_plot(p1,4,1.3*maxY,dpt,1,11,0,'k',0)
pvalues_plot(p2,5,1.3*maxY,dpt,1,11,0,'k',0)

%succ-fail 
max_early_diff.diff = max_early_diff.success-max_early_diff.failure;
max_late_diff.diff = max_late_diff.success-max_late_diff.failure;

meanPlot1 = mean(max_early_diff.diff,'omitnan');
h(1)=bar([7],meanPlot1,'b');
meanPlot2 = mean(max_late_diff.diff,'omitnan');
h(2)=bar([8],meanPlot2,'r');
error1H = std(max_early_diff.diff,'omitnan')/sqrt(sum(~isnan(max_early_diff.diff)));
error1L = std(max_early_diff.diff,'omitnan')/sqrt(sum(~isnan(max_early_diff.diff)));
error2H = std(max_late_diff.diff,'omitnan')/sqrt(sum(~isnan(max_late_diff.diff)));
error2L = std(max_late_diff.diff,'omitnan')/sqrt(sum(~isnan(max_late_diff.diff)));

errorbar([7,8], [meanPlot1,meanPlot2],[error1L,error2L],[error1H,error2H],'ko');

% [pval(3)]=ranksum(max_early_diff.diff,max_late_diff.diff);
maxY = max([nanmean(max_early_diff.diff),nanmean(max_late_diff.diff)]);
% pvalues_plot(pval(3),7.5,1.6*maxY,dpt,1,14,0.5,'k',0)
% [p1]=signrank(max_early_diff.diff);
% [p2]=signrank(max_late_diff.diff);
% pvalues_plot(p1,7,1.2*maxY,dpt,1,11,0,'k',0)
% pvalues_plot(p2,8,1.2*maxY,dpt,1,11,0,'k',0)

ylabel('M1->DLS - DLS->M1 [bits]')
set(gca, 'XTick',[1.5,4.5,7.5], 'XTickLabel',{'Success','Failure','Succ-Fail'})
title('Info flow directionality within success and failure trials')
% ylim([-70,70])

if plot_params.save == 1
    print([plot_params.save_path '\PID_sharedInfo_pooled_early_late_Histogram.svg'],'-painters','-dsvg');
    print([plot_params.save_path '\PID_sharedInfo_pooled_early_late_Histogram.png'],'-dpng');
end  

%% Reversal magnitude success vs failure (Supplemental figure 3 R)

% Info flow 
nBoot = 100000;
reversal.success = mean(max_late_diff.success)-mean(max_early_diff.success);
reversal.failure = mean(max_late_diff.failure)-mean(max_early_diff.failure);

clear bootSuccE bootSuccL bootFailE bootFailL bootReversal_succ bootReversal_fail
for booti=1:nBoot
    bootSuccE(booti) = mean(datasample(max_early_diff.success,numel(max_early_diff.success)));
    bootSuccL(booti) = mean(datasample(max_late_diff.success,numel(max_late_diff.success)));
    bootReversal_succ(booti) = bootSuccL(booti)-bootSuccE(booti);

    bootFailE(booti) = mean(datasample(max_early_diff.failure,numel(max_early_diff.failure)));
    bootFailL(booti) = mean(datasample(max_late_diff.failure,numel(max_late_diff.failure)));
    bootReversal_fail(booti) = bootFailL(booti)-bootFailE(booti);
end 

figure()
hold on
bar([1],reversal.success,'g')
bar([2],reversal.failure,'FaceColor',[0.7,0.7,0.7])
legend('Success','Failure','AutoUpdate','off')
error1H = prctile(bootReversal_succ,97.5)-reversal.success;
error1L = reversal.success-prctile(bootReversal_succ,2.5);
error2H = prctile(bootReversal_fail,97.5)-reversal.failure;
error2L = reversal.failure-prctile(bootReversal_fail,2.5);
errorbar([1,2],[reversal.success,reversal.failure],[error1L,error2L],[error1H,error2H],'ko')

maxY = max([reversal.success+error1H, reversal.failure+error2H]);
p1 = 2*min([mean(bootReversal_succ<0), mean(bootReversal_succ>0)]); % evaluate the magnitude of both tails and multiply by 2 to do a two tailed test
p2 = 2*min([mean(bootReversal_fail<0), mean(bootReversal_fail>0)]); % evaluate the magnitude of both tails and multiply by 2 to do a two tailed test
p_diff = 2*min([mean(bootReversal_succ-bootReversal_fail<0), mean(bootReversal_succ-bootReversal_fail>0)]); % evaluate the magnitude of both tails and multiply by 2 to do a two tailed test
pvalues_plot(p1,1,1.1*maxY,maxY/10,1,12,0,'k',0)
pvalues_plot(p2,2,1.1*maxY,maxY/10,1,12,0,'k',0)
pvalues_plot(p_diff,1.5,1.2*maxY,maxY/10,1,14,0.5,'k',0)
ylim([0,1.4*maxY])

ylabel('Info flow reversal [bits]')
title('Within succes vs failure reversal (info flow)')

% Delays 
reversal.success = mean(early_delays_succ,'omitnan')-mean(late_delays_succ,'omitnan');
reversal.failure = mean(early_delays_fail,'omitnan')-mean(late_delays_fail,'omitnan');

clear bootSuccE bootSuccL bootFailE bootFailL bootReversal_succ bootReversal_fail
for booti=1:nBoot
    bootSuccE(booti) = mean(datasample(early_delays_succ,numel(early_delays_succ)),'omitnan');
    bootSuccL(booti) = mean(datasample(late_delays_succ,numel(late_delays_succ)),'omitnan');
    bootReversal_succ(booti) = bootSuccE(booti)-bootSuccL(booti);

    bootFailE(booti) = mean(datasample(early_delays_fail,numel(early_delays_fail)),'omitnan');
    bootFailL(booti) = mean(datasample(late_delays_fail,numel(late_delays_fail)),'omitnan');
    bootReversal_fail(booti) = bootFailE(booti)-bootFailL(booti);
end 

figure()
hold on
bar([1],reversal.success,'g')
bar([2],reversal.failure,'FaceColor',[0.7,0.7,0.7])
legend('Success','Failure','AutoUpdate','off')
error1H = prctile(bootReversal_succ,97.5)-reversal.success;
error1L = reversal.success-prctile(bootReversal_succ,2.5);
error2H = prctile(bootReversal_fail,97.5)-reversal.failure;
error2L = reversal.failure-prctile(bootReversal_fail,2.5);
errorbar([1,2],[reversal.success,reversal.failure],[error1L,error2L],[error1H,error2H],'ko')

maxY = max([reversal.success+error1H, reversal.failure+error2H]);
p1 = 2*min([mean(bootReversal_succ<0), mean(bootReversal_succ>0)]); % evaluate the magnitude of both tails and multiply by 2 to do a two tailed test
p2 = 2*min([mean(bootReversal_fail<0), mean(bootReversal_fail>0)]); % evaluate the magnitude of both tails and multiply by 2 to do a two tailed test
p_diff = 2*min([mean(bootReversal_succ-bootReversal_fail<0), mean(bootReversal_succ-bootReversal_fail>0)]); % evaluate the magnitude of both tails and multiply by 2 to do a two tailed test
pvalues_plot(p1,1,1.1*maxY,maxY/10,1,12,0,'k',0)
pvalues_plot(p2,2,1.1*maxY,maxY/10,1,12,0,'k',0)
pvalues_plot(p_diff,1.5,1.2*maxY,maxY/10,1,14,0.5,'k',0)
ylim([0,1.4*maxY])
ylabel('Delays reversal [ms]')
title('Within succes vs failure reversal (delays)')

%%

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end

end