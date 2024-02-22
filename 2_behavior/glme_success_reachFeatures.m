function [] = glme_success_reachFeatures(reach_features,params,save_params)

% all animal model

    feature_pVal = [];
    feature_tStat = [];
    feature_deviance = [];
    feature_rSquared = [];
    feature_name = {};

    for rfIdx = [1 2 3 4 9 5 6 10 13 12]%1:numel(params.reachFeatures)

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

    figure;
        hold on;
        [~,idx] = sort(mean([feature_tStat(:,1) feature_tStat(:,2)],2)); 
        plot(feature_tStat(flip(idx),1),'k'); 
        scatter(1:10,feature_tStat(flip(idx),1),30,[0 0 0],'filled');
        plot(feature_tStat(flip(idx),2),'r');
        scatter(1:10,feature_tStat(flip(idx),2),30,[1 0 0],'filled');
        ylim([-5 15]);
        xticklabels(feature_name(flip(idx)))
        if save_params.save == 1
            print([save_params.save_path '\all_animal_glme.eps'],'-painters','-depsc');
        end

    figure;
        hold on;
        plot(log10(feature_pVal(flip(idx),1)),'k'); 
        scatter(1:10,log10(feature_pVal(flip(idx),1)),30,[0 0 0],'filled');
        plot(log10(feature_pVal(flip(idx),2)),'r');
        scatter(1:10,log10(feature_pVal(flip(idx),2)),30,[1 0 0],'filled');
        xticklabels(feature_name(flip(idx)));     
        if save_params.save == 1
            print([save_params.save_path '\all_animal_glme_logpVal.eps'],'-painters','-depsc');
        end

% individual animal models

    features_used = [1 2 3 4 9 5 6 10 13 12];
    features_used = features_used(flip(idx));
    
    figure
    hold on

    all_feature_p = [];
    countPlot = 1;
    for rfIdx = features_used
        
        feature_pVal = [];
        feature_tStat = [];
        feature_deviance = [];
        feature_rSquared = [];

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
    
            T = table(tmp_early_success(1:min_trials)',tmp_early_features(1:min_trials)',tmp_early_ratID(1:min_trials)','VariableNames',{'success','reachFeature','ratID'});
            Eglme = fitglme(T,'success ~ reachFeature + (1|ratID)','Distribution','binomial');
    
            T = table(tmp_late_success(end-min_trials+1:end)',tmp_late_features(end-min_trials+1:end)',tmp_late_ratID(end-min_trials+1:end)','VariableNames',{'success','reachFeature','ratID'});
            Lglme = fitglme(T,'success ~ reachFeature + (1|ratID)','Distribution','binomial');

            feature_pVal = [feature_pVal; Eglme.Coefficients(2,6).pValue Lglme.Coefficients(2,6).pValue];
            feature_tStat = [feature_tStat; Eglme.Coefficients(2,4).tStat Lglme.Coefficients(2,4).tStat];
            feature_deviance = [feature_deviance; Eglme.ModelCriterion(1,4).Deviance Lglme.ModelCriterion(1,4).Deviance];
            feature_rSquared = [feature_rSquared; Eglme.Rsquared.Ordinary Lglme.Rsquared.Ordinary];

        end

        scatter(countPlot*ones(size(feature_tStat,1),1),feature_tStat(:,1),30,[.5 .5 .5],'filled');
        errorbar(countPlot,mean(feature_tStat(:,1)),std(feature_tStat(:,1))/sqrt(size(feature_tStat,1)),'color','k','LineWidth',3,'CapSize',20);
        scatter(countPlot*ones(size(feature_tStat,1),1),feature_tStat(:,2),30,[0.75 0 0],'filled');
        errorbar(countPlot,mean(feature_tStat(:,2)),std(feature_tStat(:,2))/sqrt(size(feature_tStat,1)),'color',[.8 0 0],'LineWidth',3,'CapSize',20);
        [p h] = signrank(feature_tStat(:,1),feature_tStat(:,2));

        params.reachFeatures{rfIdx};
        [h p ci stats] = ttest(feature_tStat(:,1),feature_tStat(:,2));

        disp([params.reachFeatures{rfIdx} ' | ' ...
              'naive: t(' num2str(stats.df) ')=' num2str(stats.tstat) ', P=' num2str(p)]);
        
        all_feature_p = [all_feature_p p];
        if p<0.05
            scatter(countPlot,10,50,[0 0 0],'filled');
        end
        countPlot = countPlot+1;

    end

    xticklabels(feature_name(flip(idx)));

    if save_params.save == 1
        print([save_params.save_path '\by_animal_glme.eps'],'-painters','-depsc');
    end

end