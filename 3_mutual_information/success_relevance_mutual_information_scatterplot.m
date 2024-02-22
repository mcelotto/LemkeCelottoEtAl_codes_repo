function [] = success_relevance_mutual_information_scatterplot(reach_features,MIout,LFP_early_late,params,clusterParams,save_params,newBinTimes)

%%% GLME ALL ANIMAL MODEL

    feature_pVal = [];
    feature_tStat = [];
    feature_deviance = [];
    feature_rSquared = [];
    feature_name = {};

    for rfIdx = 1:14

        feature_name = [feature_name params.reachFeatures{rfIdx}];
        
        early_success = [];
        early_features = [];
        early_ratID = [];

        late_success = [];
        late_features = [];
        late_ratID = [];

        for animalIdx = 1:length(params.animals)

            tmp_early_success = [];
            tmp_early_features = [];
            tmp_early_ratID = [];
    
            tmp_late_success = [];
            tmp_late_features = [];
            tmp_late_ratID = [];

            for dayIdx = 1:params.num_earlylate_days{animalIdx}
                tmp_early_success = [tmp_early_success reach_features{animalIdx}(dayIdx).success_rate];
                tmp_early_features = [tmp_early_features reach_features{animalIdx}(dayIdx).(params.reachFeatures{rfIdx})];
                tmp_early_ratID = [tmp_early_ratID animalIdx*ones(1,length(reach_features{animalIdx}(dayIdx).success_rate))];
            end

            for dayIdx = numel(reach_features{animalIdx})-params.num_earlylate_days{animalIdx}+1:numel(reach_features{animalIdx})
                tmp_late_success = [tmp_late_success reach_features{animalIdx}(dayIdx).success_rate];
                tmp_late_features = [tmp_late_features reach_features{animalIdx}(dayIdx).(params.reachFeatures{rfIdx})];
                tmp_late_ratID = [tmp_late_ratID animalIdx*ones(1,length(reach_features{animalIdx}(dayIdx).success_rate))];
            end

            early_idx = (~isnan(tmp_early_success) & ~isnan(tmp_early_features) & ~isnan(tmp_early_ratID));
            tmp_early_success = tmp_early_success(early_idx);
            tmp_early_features = tmp_early_features(early_idx);
            tmp_early_ratID = tmp_early_ratID(early_idx);

            late_idx = (~isnan(tmp_late_success) & ~isnan(tmp_late_features) & ~isnan(tmp_late_ratID));
            tmp_late_success = tmp_late_success(late_idx);
            tmp_late_features = tmp_late_features(late_idx);
            tmp_late_ratID = tmp_late_ratID(late_idx);

            min_trials = min([sum(~isnan(tmp_early_features)) sum(~isnan(tmp_late_features))]);

            early_success = [early_success tmp_early_success(1:min_trials)];
            early_features = [early_features tmp_early_features(1:min_trials)];
            early_ratID = [early_ratID tmp_early_ratID(1:min_trials)];
    
            late_success = [late_success tmp_late_success(end-min_trials+1:end)];
            late_features = [late_features tmp_late_features(end-min_trials+1:end)];
            late_ratID = [late_ratID tmp_late_ratID(end-min_trials+1:end)];

        end

        T = table(early_success',early_features',early_ratID','VariableNames',{'success','reachFeature','ratID'});
        Eglme = fitglme(T,'success ~ reachFeature + (1|ratID)','Distribution','binomial');

        T = table(late_success',late_features',late_ratID','VariableNames',{'success','reachFeature','ratID'});
        Lglme = fitglme(T,'success ~ reachFeature + (1|ratID)','Distribution','binomial');

        disp([params.reachFeatures{rfIdx} ' | ' ...
              'naive: t(' num2str(Eglme.Coefficients(2,5).DF) ')=' num2str(Eglme.Coefficients(2,4).tStat) ', P=' num2str(Eglme.Coefficients(2,6).pValue) ...
              ' | skilled: t(' num2str(Lglme.Coefficients(2,5).DF) ')=' num2str(Lglme.Coefficients(2,4).tStat) ', P=' num2str(Lglme.Coefficients(2,6).pValue) ...
            ]);
        
        feature_pVal = [feature_pVal; Eglme.Coefficients(2,6).pValue Lglme.Coefficients(2,6).pValue];
        feature_tStat = [feature_tStat; Eglme.Coefficients(2,4).tStat Lglme.Coefficients(2,4).tStat];
        feature_deviance = [feature_deviance; Eglme.ModelCriterion(1,4).Deviance Lglme.ModelCriterion(1,4).Deviance];
        feature_rSquared = [feature_rSquared; Eglme.Rsquared.Ordinary Lglme.Rsquared.Ordinary];

    end

%%% REACH FEATURES

    allFeatureInfo_M1_naive = cell(2,length(params.reachFeatures)-1);
    allFeatureInfo_M1_skilled = cell(2,length(params.reachFeatures)-1);
    allFeatureInfo_DLS_naive = cell(2,length(params.reachFeatures)-1);
    allFeatureInfo_DLS_skilled = cell(2,length(params.reachFeatures)-1);
    
    for n_lfp = 1:length(params.lfpFeatures)
        for fidx = 1:length(params.reachFeatures)-1 % not computed for success rate
          
        info_mean_all.M1 = cell(1,2);
        info_mean_all.DLS = cell(1,2);        
        info_mean_sig.M1 = cell(1,2);
        info_mean_sig.DLS = cell(1,2);
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_all.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                                           
                        info_mean_sig.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
                    end
                end
            end
        end
        
        m1_early_info_all = [];
        m1_early_info = [];
        for unit = 1:size(info_mean_sig.M1{1},1)
            tmp_info = info_mean_sig.M1{1}(unit,:);
            tmp_info2 = info_mean_all.M1{1}(unit,:);
            [v,i] = max(tmp_info);        
            m1_early_info_all = [m1_early_info_all; tmp_info];
            if v>0.001 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info];
            end
        end
    
        m1_late_info_all = [];
        m1_late_info = [];
        for unit = 1:size(info_mean_sig.M1{2},1)
            tmp_info = info_mean_sig.M1{2}(unit,:);
            tmp_info2 = info_mean_all.M1{2}(unit,:);
            [v,i] = max(tmp_info);
            m1_late_info_all = [m1_late_info_all; tmp_info];        
            if v>0.001 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info];
            end
        end
    
        dls_early_info_all = [];
        dls_early_info = [];
        for unit = 1:size(info_mean_sig.DLS{1},1)
            tmp_info = info_mean_sig.DLS{1}(unit,:);
            tmp_info2 = info_mean_all.DLS{1}(unit,:);
            [v,i] = max(tmp_info);        
            dls_early_info_all = [dls_early_info_all; tmp_info];
            if v>0.001 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info];
            end
        end
    
        dls_late_info_all = [];
        dls_late_info = [];
        for unit = 1:size(info_mean_sig.DLS{2},1)
            tmp_info = info_mean_sig.DLS{2}(unit,:);
            tmp_info2 = info_mean_all.DLS{2}(unit,:);
            [v,i] = max(tmp_info);        
            dls_late_info_all = [dls_late_info_all; tmp_info];
            if v>0.001 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info];
            end
        end         
        
        allFeatureInfo_M1_naive{1,fidx} = m1_early_info_all;
        allFeatureInfo_M1_naive{2,fidx} = m1_early_info;
        allFeatureInfo_M1_skilled{1,fidx} = m1_late_info_all;
        allFeatureInfo_M1_skilled{2,fidx} = m1_late_info;
           
        allFeatureInfo_DLS_naive{1,fidx} = dls_early_info_all;
        allFeatureInfo_DLS_naive{2,fidx} = dls_early_info;
        allFeatureInfo_DLS_skilled{1,fidx} = dls_late_info_all;
        allFeatureInfo_DLS_skilled{2,fidx} = dls_late_info;
          
        disp(['LFP | ' params.reachFeatures{fidx} ' | M1 | early: ' num2str(size(m1_early_info,1)) ' out of ' num2str(size(m1_early_info_all,1)) ' - ' num2str(100*size(m1_early_info,1)/size(m1_early_info_all,1))  ...
                                                      '% | late: ' num2str(size(m1_late_info,1)) ' out of ' num2str(size(m1_late_info_all,1)) ' - ' num2str(100*size(m1_late_info,1)/size(m1_late_info_all,1))  '%'])
        disp(['LFP | ' params.reachFeatures{fidx} ' | DLS | early: ' num2str(size(dls_early_info,1)) ' out of ' num2str(size(dls_early_info_all,1)) ' - ' num2str(100*size(dls_early_info,1)/size(dls_early_info_all,1))  ...
                                                      '% | late: ' num2str(size(dls_late_info,1)) ' out of ' num2str(size(dls_late_info_all,1)) ' - ' num2str(100*size(dls_late_info,1)/size(dls_late_info_all,1))  '%'])
        end
    end

%%% PLOT

    all_m1_early_mean = [];
    for n = 1:14
        all_m1_early_mean = [all_m1_early_mean mean(max(allFeatureInfo_M1_naive{1,n}(:,12:86)'))];
    end
    all_m1_late_mean = [];
    for n = 1:14
        all_m1_late_mean = [all_m1_late_mean mean(max(allFeatureInfo_M1_skilled{1,n}(:,12:86)'))];
    end
    all_dls_early_mean = [];
    for n = 1:14
        all_dls_early_mean = [all_dls_early_mean mean(max(allFeatureInfo_DLS_naive{1,n}(:,12:86)'))];
    end
    all_dls_late_mean = [];
    for n = 1:14
        all_dls_late_mean = [all_dls_late_mean mean(max(allFeatureInfo_DLS_skilled{1,n}(:,12:86)'))];
    end
    
    figure;
    hold on;
    count = 1;
    all_x = [];
    all_y = [];
    for n = [13 9 3 4 10 6 5 12 2 1]
        scatter(all_m1_early_mean(n),feature_tStat(n,1),30,[count/11 0 1-(count/11)])
        scatter(all_dls_early_mean(n),feature_tStat(n,1),30,[count/11 0 1-(count/11)],'diamond')
        scatter(all_m1_late_mean(n),feature_tStat(n,2),30,[count/11 0 1-(count/11)],'filled')
        scatter(all_dls_late_mean(n),feature_tStat(n,2),30,[count/11 0 1-(count/11)],'filled','diamond')
        all_x = [all_x all_m1_early_mean(n) all_dls_early_mean(n) all_m1_late_mean(n) all_dls_late_mean(n)]
        all_y = [all_y feature_tStat(n,1) feature_tStat(n,1) feature_tStat(n,2) feature_tStat(n,2)]
        count = 1 + count; 
%         pause;
    end
    [R,P] = corrcoef(all_x,all_y)
    [p,S] = polyfit(all_x,all_y,1)
    [y_fit,delta] = polyval(p,all_x,S);
    plot(all_x,y_fit,'k')
    title(num2str(P(1,2)))

    if save_params.save == 1
        print([save_params.save_path '\info_success_relevance_scatter_LFP.eps'],'-painters','-depsc');
    end
    
%     scatter(all_m1_early_mean([9]),feature_tStat([9],1),30,[1 0 0],'filled')
%     scatter(all_dls_early_mean([9]),feature_tStat([9],1),30,[.75 0 0],'filled')
%     scatter(all_m1_late_mean([9]),feature_tStat([9],2),30,[1 0 0],'filled')
%     scatter(all_dls_late_mean([9]),feature_tStat([9],2),30,[.75 0 0],'filled')
% 
%     scatter(all_m1_early_mean([13]),feature_tStat([13],1),30,[0 1 0],'filled')
%     scatter(all_dls_early_mean([13]),feature_tStat([13],1),30,[0 .75 0],'filled')
%     scatter(all_m1_late_mean([13]),feature_tStat([13],2),30,[0 1 0],'filled')
%     scatter(all_dls_late_mean([13]),feature_tStat([13],2),30,[0 .75 0],'filled')
% 
%     scatter(all_m1_early_mean(10),feature_tStat(10,1),30,[0 0 1],'filled')
%     scatter(all_dls_early_mean(10),feature_tStat(10,1),30,[0 0 .5],'filled')
%     scatter(all_m1_late_mean(10),feature_tStat(10,2),30,[0 0 1],'filled')
%     scatter(all_dls_late_mean(10),feature_tStat(10,2),30,[0 0 .5],'filled')
% 
%     scatter(all_m1_early_mean(12),feature_tStat(12,1),30,[0 0 0],'filled')
%     scatter(all_dls_early_mean(12),feature_tStat(12,1),30,[.75 .75 .75],'filled')
%     scatter(all_m1_late_mean(12),feature_tStat(12,2),30,[0 0 0],'filled')
%     scatter(all_dls_late_mean(12),feature_tStat(12,2),30,[.75 .75 .75],'filled')
% 
%     figure;
%     hold on;
%     scatter(all_m1_early_mean([1:13]),feature_tStat([1:13],1),30,[.5 .5 .5],'filled')
%     scatter(all_dls_early_mean([1:13]),feature_tStat([1:13],1),30,[1 .5 .5],'filled')
%     scatter(all_m1_late_mean([1:13]),feature_tStat([1:13],2),30,[0 0 0],'filled')
%     scatter(all_dls_late_mean([1:13]),feature_tStat([1:13],2),30,[1 0 0],'filled')

end

