function [] = plot_PID_Spikes_newChans_examples(PIDout,params,reach_bins,features_to_plot,plot_params,clusterParams)

    animal_colors = distinguishable_colors(numel(params.animals));
    
    for fidx = features_to_plot%1:length(params.reachFeatures)
        
        early_info_all = [];
        late_info_all = [];
        early_info_all_noZero = [];
        late_info_all_noZero = [];        
        for animal = 1:length(PIDout)
        
            if ~isempty(PIDout{animal})
        
                % early
                for day = 1:params.num_earlylate_days{animal}
                    if day<=length(PIDout{animal}.spikes)
                        if isfield(PIDout{animal}.spikes{day},(params.reachFeatures{fidx}))
                            for unit = 1:size(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared,1)
                                if sum(sum(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:)))>0
                                    tmp_info = squeeze(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:));
                                    tmp_infoSh = squeeze(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).sharedSh(unit,:,:,:));
                                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                                    early_info_all_noZero = cat(3,early_info_all_noZero,tmp_info);                                    
                                    tmp_info(~tmp_sig) = 0;
                                    early_info_all = cat(3,early_info_all,tmp_info);
                                end
                            end
                        end
                    end
                end
        
                % late
                for day = length(params.days{animal})-params.num_earlylate_days{animal}+1:length(params.days{animal})
                    if day<=length(PIDout{animal}.spikes)
                        if isfield(PIDout{animal}.spikes{day},(params.reachFeatures{fidx}))
                            for unit = 1:size(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared,1)
                                if sum(sum(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:)))>0
                                    tmp_info = squeeze(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).shared(unit,:,:));
                                    tmp_infoSh = squeeze(PIDout{animal}.spikes{day}.(params.reachFeatures{fidx}).sharedSh(unit,:,:,:));
                                    tmp_sig = clusterStat_v2(tmp_info,tmp_infoSh,clusterParams(1),clusterParams(2));
                                    tmp_info = tmp_info-squeeze(mean(tmp_infoSh,3));
                                    late_info_all_noZero = cat(3,late_info_all_noZero,tmp_info);                                    
                                    tmp_info(~tmp_sig) = 0;
                                    late_info_all = cat(3,late_info_all,tmp_info);
                                end
                            end
                        end
                    end
                end
            end
        end
              
%         for pair = 1:size(early_info_all,3)
%             [v, i] = max(median(squeeze(early_info_all(reach_bins,:,pair)),1));
%             if v>0 & i>1& i<51
%                 figure;
%                 hold on;
%                 imagesc(squeeze(early_info_all_noZero(:,:,pair))')
%                 plot([0.5 39.5],[25.5 25.5],'color','r')
%                 plot([0.5 39.5],[26.5 26.5],'color','r')
%                 plot([26 26],[.5 51.5],'color','k')
%                 colorbar;
%                 xlim([6 36])
%                 xticks([6 16 26 36])
%                 xticklabels([-1 -.5 0 .5])
%                 xlabel('time from pellet touch (ms)')
%                 ylim([.5 51.5])
%                 yticks(1:2:51)
%                 yticklabels(-250:20:250)
%                 ylabel('DLS time delay (ms)')
%                 title(['naive example ' num2str(pair)])
%             end
%         end
%         
%         for pair = 1:size(late_info_all,3)
%             [v, i] = max(median(squeeze(late_info_all(reach_bins,:,pair)),1));
%             if v>0 & i>1 & i<51
%                 figure;
%                 hold on;
%                 imagesc(squeeze(late_info_all_noZero(:,:,pair))')
%                 plot([0.5 39.5],[25.5 25.5],'color','r')
%                 plot([0.5 39.5],[26.5 26.5],'color','r')
%                 plot([26 26],[.5 51.5],'color','k')
%                 colorbar;
%                 xlim([6 36])
%                 xticks([6 16 26 36])
%                 xticklabels([-1 -.5 0 .5])
%                 xlabel('time from pellet touch (ms)')
%                 ylim([.5 51.5])
%                 yticks(1:2:51)
%                 yticklabels(-250:20:250)
%                 ylabel('DLS time delay (ms)')
%                 title(['skilled example ' num2str(pair)])
%             end
%         end

        if plot_params.save == 1

            %%% EARLY 1
                figure;
                hold on;
                imagesc(squeeze(early_info_all_noZero(:,:,221))',[0 0.013])
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
                title(['naive example ' num2str(221)])
%                 print([plot_params.save_path '\PID_naive_example_' num2str(221) '.eps'],'-painters','-depsc');
  
            %%% EARLY 2
                figure;
                hold on;
                imagesc(squeeze(early_info_all_noZero(:,:,224))',[0 0.015])
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
                title(['naive example ' num2str(224)])
                print([plot_params.save_path '\PID_naive_example_' num2str(224) '.eps'],'-painters','-depsc');
  
            %%% EARLY 3
                figure;
                hold on;
                imagesc(squeeze(early_info_all_noZero(:,:,228))',[0 0.0075])
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
                title(['naive example ' num2str(228)])
                print([plot_params.save_path '\PID_naive_example_' num2str(228) '.eps'],'-painters','-depsc');
  
            %%% SKILLED 1
                figure;
                hold on;
                imagesc(squeeze(late_info_all_noZero(:,:,197))',[0 0.012])
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
                title(['skilled example ' num2str(197)])
                print([plot_params.save_path '\PID_skilled_example_' num2str(197) '.eps'],'-painters','-depsc');
  
            %%% SKILLED 2
                figure;
                hold on;
                imagesc(squeeze(late_info_all_noZero(:,:,196))',[0 0.01])
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
                title(['skilled example ' num2str(196)])
                print([plot_params.save_path '\PID_skilled_example_' num2str(196) '.eps'],'-painters','-depsc');
  
            %%% SKILLED 3
                figure;
                hold on;
                imagesc(squeeze(late_info_all_noZero(:,:,258))',[0 0.016])
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
                title(['skilled example ' num2str(258)])
                print([plot_params.save_path '\PID_skilled_example_' num2str(258) '.eps'],'-painters','-depsc');
  
        end




    end
end
