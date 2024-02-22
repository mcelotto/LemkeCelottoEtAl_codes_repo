function plot_linearRegression_LFP(params, LFP_early_late, reachFeatures_early_late, save_params, features_to_plot)

for fidx = features_to_plot

    all_animal_early_pVal_M1 = [];
    all_animal_late_pVal_M1 = [];
    all_animal_early_pVal_DLS = [];
    all_animal_late_pVal_DLS = [];

    all_animal_early_rVal_M1 = [];
    all_animal_late_rVal_M1 = [];
    all_animal_early_rVal_DLS = [];
    all_animal_late_rVal_DLS = [];

    for animal = 1:8
        
        disp(['Animal ' num2str(animal)]);

        pVal_M1 = nan(2,size(LFP_early_late{animal}{1}.reach_related_lfp_ref.M1,1),250);
        rVal_M1 = nan(2,size(LFP_early_late{animal}{1}.reach_related_lfp_ref.M1,1),250);
        for timeBin = 51:200
            for early_late = 1:2
                for channel = 1:size(LFP_early_late{animal}{early_late}.reach_related_lfp_ref.M1,1)
                    LFPsignal = squeeze(LFP_early_late{animal}{early_late}.reach_related_lfp_ref.M1(channel,timeBin-2:timeBin,:))';
                    movFeature = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx});
                    mdl = fitlm(LFPsignal,movFeature);
                    
                    rVal_M1(early_late,channel,timeBin) = mdl.Rsquared.Ordinary;    
                    pVal_M1(early_late,channel,timeBin) = coefTest(mdl);    
                end
            end
        end
        all_animal_early_rVal_M1 = [all_animal_early_rVal_M1; squeeze(rVal_M1(1,:,:))];
        all_animal_late_rVal_M1 = [all_animal_late_rVal_M1; squeeze(rVal_M1(2,:,:))];
        all_animal_early_pVal_M1 = [all_animal_early_pVal_M1; squeeze(pVal_M1(1,:,:))];
        all_animal_late_pVal_M1 = [all_animal_late_pVal_M1; squeeze(pVal_M1(2,:,:))];
        
        pVal_DLS = nan(2,size(LFP_early_late{animal}{1}.reach_related_lfp_ref.DLS,1),250);
        rVal_DLS = nan(2,size(LFP_early_late{animal}{1}.reach_related_lfp_ref.DLS,1),250);
        for timeBin = 51:200
            for early_late = 1:2
                for channel = 1:size(LFP_early_late{animal}{early_late}.reach_related_lfp_ref.DLS,1)
                    LFPsignal = squeeze(LFP_early_late{animal}{early_late}.reach_related_lfp_ref.DLS(channel,timeBin-2:timeBin,:))';
                    movFeature = reachFeatures_early_late{animal}{early_late}.(params.reachFeatures{fidx});
                    mdl = fitlm(LFPsignal,movFeature);
                    
                    rVal_DLS(early_late,channel,timeBin) = mdl.Rsquared.Ordinary;    
                    pVal_DLS(early_late,channel,timeBin) = coefTest(mdl);      
                end
            end
        end
        all_animal_early_rVal_DLS = [all_animal_early_rVal_DLS; squeeze(rVal_DLS(1,:,:))];
        all_animal_late_rVal_DLS = [all_animal_late_rVal_DLS; squeeze(rVal_DLS(2,:,:))];
        all_animal_early_pVal_DLS = [all_animal_early_pVal_DLS; squeeze(pVal_DLS(1,:,:))];
        all_animal_late_pVal_DLS = [all_animal_late_pVal_DLS; squeeze(pVal_DLS(2,:,:))];

    end
    
	%%% PEAK SHIFTS

        [edge_M1_early_val, edge_M1_early] = min(all_animal_early_pVal_M1(:,51:200)');
        edge_M1_early = edge_M1_early(edge_M1_early_val<0.01);
        
        [edge_M1_late_val, edge_M1_late] = min(all_animal_late_pVal_M1(:,51:200)');
        edge_M1_late = edge_M1_late(edge_M1_late_val<0.01);

        [edge_DLS_early_val, edge_DLS_early] = min(all_animal_early_pVal_DLS(:,51:200)');
        edge_DLS_early = edge_DLS_early(edge_DLS_early_val<0.01);

        [edge_DLS_late_val, edge_DLS_late] = min(all_animal_late_pVal_DLS(:,51:200)');
        edge_DLS_late = edge_DLS_late(edge_DLS_late_val<0.01);

    %%% PEAK SHIFT M1

        figure('units','normalized','outerposition',[0 0 1 1]);
            clearvars g;

            peak1_i = edge_M1_early;
            peak2_i = edge_M1_late;

            grouping = [ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))];
            ID = cell(1,length([ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))]));
            ID(1:length(peak1_i)) = {'Naive'};
            ID(1+length(peak1_i):length(grouping)) = {'Learned'};

            g(1,1)=gramm('x',ID,'y',[peak1_i peak2_i],'color',grouping);
            g(1,1).stat_boxplot();
            g(1,1).coord_flip();
            g(1,1).axe_property('YLim',[1 150]);
            g(2,1)=gramm('x',[peak1_i peak2_i],'color',grouping);
            g(2,1).stat_density('bandwidth',4);
            g(2,1).axe_property('XLim',[1 150]);

            [h, p, ci, stats] = ttest2(peak1_i,peak2_i);
            [p2, ~, ~] = ranksum(peak1_i,peak2_i);
            disp(['M1 shift: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

            g.set_title(['M1 shift | ' params.reachFeatures{fidx} ' | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
            g.draw();
            if save_params.save == 1
                print([save_params.save_path '\LFP_M1_naiveSkilled_linearRegression_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
            end

    %%% PEAK SHIFT DLS

        figure('units','normalized','outerposition',[0 0 1 1]);
            clearvars g;

            peak1_i = edge_DLS_early;
            peak2_i = edge_DLS_late;

            grouping = [ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))];
            ID = cell(1,length([ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))]));
            ID(1:length(peak1_i)) = {'Naive'};
            ID(1+length(peak1_i):length(grouping)) = {'Learned'};

            g(1,1)=gramm('x',ID,'y',[peak1_i peak2_i],'color',grouping);
            g(1,1).stat_boxplot();
            g(1,1).coord_flip();
            g(1,1).axe_property('YLim',[1 150]); 
            g(2,1)=gramm('x',[peak1_i peak2_i],'color',grouping);
            g(2,1).stat_density('bandwidth',4);
            g(2,1).axe_property('XLim',[1 150]); 

            [h, p, ci, stats] = ttest2(peak1_i,peak2_i);
            [p2, ~, ~] = ranksum(peak1_i,peak2_i);
            disp(['DLS shift: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p) ' - ranksum p: ' num2str(p2)]);

            g.set_title(['DLS shift | ' params.reachFeatures{fidx} ' | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
            g.draw();
            if save_params.save == 1
                print([save_params.save_path '\LFP_DLS_naiveSkilled_linearRegression_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
            end  

    %%% PEAK SHIFT NAIVE

%         figure('units','normalized','outerposition',[0 0 1 1]);
%             clearvars g;
% 
%             peak1_i = edge_M1_early;
%             peak2_i = edge_DLS_early;
% 
%             grouping = [ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))];
%             ID = cell(1,length([ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))]));
%             ID(1:length(peak1_i)) = {'M1'};
%             ID(1+length(peak1_i):length(grouping)) = {'DLS'};
% 
%             g(1,1)=gramm('x',ID,'y',[peak1_i peak2_i],'color',grouping);
%             g(1,1).stat_boxplot();
%             g(1,1).coord_flip();
%             g(1,1).axe_property('YLim',[1 150]); 
%             g(2,1)=gramm('x',[peak1_i peak2_i],'color',grouping);
%             g(2,1).stat_density('bandwidth',4);
%             g(2,1).axe_property('XLim',[1 150]); 
% 
%             [h, p, ci, stats] = ttest2(peak1_i,peak2_i);
%             [p2, ~, ~] = ranksum(peak1_i,peak2_i);
%             disp(['NAIVE shift: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
% 
%             g.set_title(['NAIVE shift | ' params.reachFeatures{fidx} ' | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
%             g.draw();
%             if save_params.save == 1
%                 print([save_params.save_path '\LFP_naive_M1DLS_linearRegression_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
%             end
    
    %%% PEAK SHIFT SKILLED
    
%         figure('units','normalized','outerposition',[0 0 1 1]);
%             clearvars g;
% 
%             peak1_i = edge_M1_late;
%             peak2_i = edge_DLS_late;
% 
%             grouping = [ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))];
%             ID = cell(1,length([ones(1,length(peak1_i)) 2*ones(1,length(peak2_i))]));
%             ID(1:length(peak1_i)) = {'M1'};
%             ID(1+length(peak1_i):length(grouping)) = {'DLS'};
% 
%             g(1,1)=gramm('x',ID,'y',[peak1_i peak2_i],'color',grouping);
%             g(1,1).stat_boxplot();
%             g(1,1).coord_flip();
%             g(1,1).axe_property('YLim',[1 150]); 
%             g(2,1)=gramm('x',[peak1_i peak2_i],'color',grouping);
%             g(2,1).stat_density('bandwidth',4);
%             g(2,1).axe_property('XLim',[1 150]); 
% 
%             [h, p, ci, stats] = ttest2(peak1_i,peak2_i);
%             [p2, ~, ~] = ranksum(peak1_i,peak2_i);
%             disp(['SKILLED shift: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
% 
%             g.set_title(['SKILLED shift | ' params.reachFeatures{fidx} ' | p=' num2str(p) ' - ranksum p: ' num2str(p2)]);
%             g.draw();
%             if save_params.save == 1
%                 print([save_params.save_path '\LFP_skilled_M1DLS_linearRegression_' params.reachFeatures{fidx} '.eps'],'-painters','-depsc');
%             end

    pause(0.1);


end   