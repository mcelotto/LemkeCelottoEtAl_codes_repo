function [] = plot_beh_match_earlyTOlate(reachFeatures_early_late_noBinNoTrialMatch,reachFeatures_early_late_noBin,params,animals_to_run,features_to_run,plot_params)

for fidx = features_to_run

    all_early_match = [];
    all_late_match = [];
    all_early_all = [];
    all_late_all = [];

    for animal = animals_to_run

        if fidx==9
            dist_thresh = .54;
        elseif fidx==13
            dist_thresh = 100;
        elseif fidx==12
            dist_thresh = 100;            
        elseif fidx==14
            dist_thresh = .2;
        end

        all_early_all = [all_early_all reachFeatures_early_late_noBin{animal}{1}.(params.reachFeatures{fidx})];
        all_late_all = [all_late_all reachFeatures_early_late_noBin{animal}{2}.(params.reachFeatures{fidx})]; 
                    
        %%% match early TO late
        all_early_vals = reachFeatures_early_late_noBinNoTrialMatch{animal}{1}.(params.reachFeatures{fidx});
        matched_trial_vals = cell(1,2);
        for i = 1:length(reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx}))
            tmp_late_val = reachFeatures_early_late_noBinNoTrialMatch{animal}{2}.(params.reachFeatures{fidx})(i);
            [val, idx] = min(abs(all_early_vals-tmp_late_val));
            if val<dist_thresh
                matched_trial_vals{2} = [matched_trial_vals{2} tmp_late_val];
                matched_trial_vals{1} = [matched_trial_vals{1} all_early_vals(idx)];
            end
            all_early_vals(idx) = NaN;
        end
        all_early_match = [all_early_match matched_trial_vals{1}];
        all_late_match = [all_late_match matched_trial_vals{2}];

    end

    if fidx==9
        figure;
            subplot(1,2,1); hold on;
                histogram(all_early_all,[0:0.2:4],'DisplayStyle','Stairs');
                histogram(all_late_all,[0:0.2:4],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_all,all_late_all);
                title(['kstest: p = ' num2str(p) '|  trials: ' num2str(length(all_early_all))])
            subplot(1,2,2); hold on;
                histogram(all_early_match,[0:0.2:4],'DisplayStyle','Stairs');
                histogram(all_late_match,[0:0.2:4],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_match,all_late_match);
                title(['kstest: p = ' num2str(p) '|  trials: ' num2str(length(all_early_match))])               
            sgtitle(['matched early TO late | ' params.reachFeatures{fidx}])
            if plot_params.save == 1
                print([plot_params.save_path '\behMatch_histogram_maxVel.eps'],'-painters','-depsc');
            end
    elseif fidx==13
        figure;
            subplot(1,2,1); hold on;
                histogram(all_early_all,[0:20:500],'DisplayStyle','Stairs');
                histogram(all_late_all,[0:20:500],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_all,all_late_all);
                title(['kstest: p = ' num2str(p) '|  trials: ' num2str(length(all_early_all))])               
            subplot(1,2,2); hold on;
                histogram(all_early_match,[0:20:500],'DisplayStyle','Stairs');
                histogram(all_late_match,[0:20:500],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_match,all_late_match);
                title(['kstest: p = ' num2str(p) '|  trials: ' num2str(length(all_early_match))])               
            sgtitle(['matched early TO late | ' params.reachFeatures{fidx}])  
            if plot_params.save == 1
                print([plot_params.save_path '\behMatch_histogram_distTrav.eps'],'-painters','-depsc');
            end

    elseif fidx==12
        figure;
            subplot(1,2,1); hold on;
                histogram(all_early_all,[0:50:1000],'DisplayStyle','Stairs');
                histogram(all_late_all,[0:50:1000],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_all,all_late_all);
                title(['kstest: p = ' num2str(p) '|  trials: ' num2str(length(all_early_all))])              
            subplot(1,2,2); hold on;
                histogram(all_early_match,[0:50:1000],'DisplayStyle','Stairs');
                histogram(all_late_match,[0:50:1000],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_match,all_late_match);
                title(['kstest: p = ' num2str(p) '|  trials: ' num2str(length(all_early_match))])               
            sgtitle(['matched early TO late | ' params.reachFeatures{fidx}])                
            if plot_params.save == 1
                print([plot_params.save_path '\behMatch_histogram_movDur.eps'],'-painters','-depsc');
            end

    elseif fidx==14
        figure;
            subplot(1,2,1); hold on;
                histogram(all_early_all,[-5:0.2:5],'DisplayStyle','Stairs');
                histogram(all_late_all,[-5:0.2:5],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_all,all_late_all);
                title(['kstest: p = ' num2str(p)])                
            subplot(1,2,2); hold on;
                histogram(all_early_match,[-5:0.2:5],'DisplayStyle','Stairs');
                histogram(all_late_match,[-5:0.2:5],'DisplayStyle','Stairs');
                [h p] = kstest2(all_early_match,all_late_match);
                title(['kstest: p = ' num2str(p)])                
            sgtitle(['matched early TO late | ' params.reachFeatures{fidx}])                
    end    
end

       