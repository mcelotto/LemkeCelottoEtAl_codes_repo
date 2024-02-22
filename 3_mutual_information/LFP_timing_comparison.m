function [] = LFP_timing_comparison(MIout,LFP_early_late,params,clusterParams,f2plot,newBinTimes)
% added all time features in a single plot
% Results for individual comparison differ from the originally submitted
% paper since now we do paired statistics

sigTime_start = 12; % -1s minimum time point in the plot
sigTime_end = 86; % 0.5s minimum time point in the plot

% Time window used to select timing values
% minTime = 36; % submitted paper = 33
% maxTime = 71; % submitted paper = 85
minTime = 36; % submitted paper = 33
maxTime = 71; % submitted paper = 85
baseline_time = 1:21;

for n_lfp = 1:length(params.lfpFeatures)
    for fidx = f2plot

    %%% LOAD DATA
    
        info_mean_sig.M1 = cell(1,2);
        info_mean_sig.DLS = cell(1,2);
        raw_LFP.M1 = cell(1,2);
        raw_LFP.M1{1}=cell(1,1);
        raw_LFP.M1{2}=cell(1,1);
        raw_LFP.DLS = cell(1,2);
        raw_LFP.DLS{1}=cell(1,1);
        raw_LFP.DLS{2}=cell(1,1);
    
        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
                        infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant(~tmp_sig) = 0;      % COMMENT OUT FOR NON SHIFT PLOTS                   
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        info_mean_sig.(params.areas{aidx}){early_late} = [info_mean_sig.(params.areas{aidx}){early_late}; infQuant];
                        raw_LFP.(params.areas{aidx}){early_late}{end+1} =  zscore(squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:)))';
                    end
                end
            end
        end

        timeBins4LFP = round(newBinTimes/10)+150;

        m1_early_info = [];
        m1_early_info_all = [];
        m1_early_info_LFP = [];
        m1_early_info_LFP_z = [];
        m1_early_info_LFP_STD = [];
        m1_early_noInfo_LFP = [];
        m1_early_noInfo_LFP_z = [];
        m1_early_noInfo_LFP_STD = [];
        for unit = 1:size(info_mean_sig.M1{1},1)
            tmp_info = info_mean_sig.M1{1}(unit,:);
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.M1{1}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            tmp_LFP_STD = [];
            for bins = 1:size(timeBins4LFP,1)
                trial_sum = sum(raw_LFP.M1{1}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2))');
                tmp_LFP_STD = [tmp_LFP_STD std(trial_sum)];
            end
            m1_early_info_all = [m1_early_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_early_info = [m1_early_info; tmp_info];
                m1_early_info_LFP = [m1_early_info_LFP; tmp_LFP];
                m1_early_info_LFP_z = [m1_early_info_LFP_z; zscore(tmp_LFP)];
                m1_early_info_LFP_STD = [m1_early_info_LFP_STD; tmp_LFP_STD];
            else
                m1_early_noInfo_LFP = [m1_early_noInfo_LFP; tmp_LFP];
                m1_early_noInfo_LFP_z = [m1_early_noInfo_LFP_z; zscore(tmp_LFP)];
                m1_early_noInfo_LFP_STD = [m1_early_noInfo_LFP_STD; tmp_LFP_STD];
            end
        end

        m1_late_info = [];
        m1_late_info_all = [];
        m1_late_info_LFP = [];
        m1_late_info_LFP_z = [];
        m1_late_info_LFP_STD = [];
        m1_late_noInfo_LFP = [];
        m1_late_noInfo_LFP_z = [];
        m1_late_noInfo_LFP_STD = [];
        for unit = 1:size(info_mean_sig.M1{2},1)
            tmp_info = info_mean_sig.M1{2}(unit,:);
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.M1{2}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            tmp_LFP_STD = [];
            for bins = 1:size(timeBins4LFP,1)
                trial_sum = sum(raw_LFP.M1{2}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2))');
                tmp_LFP_STD = [tmp_LFP_STD std(trial_sum)];
            end
            m1_late_info_all = [m1_late_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                m1_late_info = [m1_late_info; tmp_info];
                m1_late_info_LFP = [m1_late_info_LFP; tmp_LFP];
                m1_late_info_LFP_z = [m1_late_info_LFP_z; zscore(tmp_LFP)];
                m1_late_info_LFP_STD = [m1_late_info_LFP_STD; tmp_LFP_STD];
            else
                m1_late_noInfo_LFP = [m1_late_noInfo_LFP; tmp_LFP];
                m1_late_noInfo_LFP_z = [m1_late_noInfo_LFP_z; zscore(tmp_LFP)];
                m1_late_noInfo_LFP_STD = [m1_late_noInfo_LFP_STD; tmp_LFP_STD];
            end
        end

        dls_early_info = [];
        dls_early_info_all = [];
        dls_early_info_LFP = [];
        dls_early_info_LFP_z = [];
        dls_early_info_LFP_STD = [];
        dls_early_noInfo_LFP = [];
        dls_early_noInfo_LFP_z = [];
        dls_early_noInfo_LFP_STD = [];
        for unit = 1:size(info_mean_sig.DLS{1},1)
            tmp_info = info_mean_sig.DLS{1}(unit,:);
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.DLS{1}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            tmp_LFP_STD = [];
            for bins = 1:size(timeBins4LFP,1)
                trial_sum = sum(raw_LFP.DLS{1}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2))');
                tmp_LFP_STD = [tmp_LFP_STD std(trial_sum)];
            end
            dls_early_info_all = [dls_early_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_early_info = [dls_early_info; tmp_info];
                dls_early_info_LFP = [dls_early_info_LFP; tmp_LFP];
                dls_early_info_LFP_z = [dls_early_info_LFP_z; zscore(tmp_LFP)];
                dls_early_info_LFP_STD = [dls_early_info_LFP_STD; tmp_LFP_STD];
            else
                dls_early_noInfo_LFP = [dls_early_noInfo_LFP; tmp_LFP];
                dls_early_noInfo_LFP_z = [dls_early_noInfo_LFP_z; zscore(tmp_LFP)];
                dls_early_noInfo_LFP_STD = [dls_early_noInfo_LFP_STD; tmp_LFP_STD];
            end
        end

        dls_late_info = [];
        dls_late_info_all = [];
        dls_late_info_LFP = [];
        dls_late_info_LFP_z = [];
        dls_late_info_LFP_STD = [];
        dls_late_noInfo_LFP = [];
        dls_late_noInfo_LFP_z = [];
        dls_late_noInfo_LFP_STD = [];
        for unit = 1:size(info_mean_sig.DLS{2},1)
            tmp_info = info_mean_sig.DLS{2}(unit,:);
            tmp_LFP = [];
            for bins = 1:size(timeBins4LFP,1)
                tmp_LFP = [tmp_LFP mean(raw_LFP.DLS{2}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
            end
            tmp_LFP_STD = [];
            for bins = 1:size(timeBins4LFP,1)
                trial_sum = sum(raw_LFP.DLS{2}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2))');
                tmp_LFP_STD = [tmp_LFP_STD std(trial_sum)];
            end
            dls_late_info_all = [dls_late_info_all; tmp_info];
            [v,i] = max(tmp_info);
            if v>0 && i>=12 && i<=86
                dls_late_info = [dls_late_info; tmp_info];
                dls_late_info_LFP = [dls_late_info_LFP; tmp_LFP];
                dls_late_info_LFP_z = [dls_late_info_LFP_z; zscore(tmp_LFP)];
                dls_late_info_LFP_STD = [dls_late_info_LFP_STD; tmp_LFP_STD];
            else
                dls_late_noInfo_LFP = [dls_late_noInfo_LFP; tmp_LFP];
                dls_late_noInfo_LFP_z = [dls_late_noInfo_LFP_z; zscore(tmp_LFP)];
                dls_late_noInfo_LFP_STD = [dls_late_noInfo_LFP_STD; tmp_LFP_STD];
            end
        end




    
    % Initialize structure used to track the identity of significant
    % channels
    sig_channels_id = cell(f2plot(end),2,2);
    for fidx = f2plot
        for aidx = 1:length(params.areas)
            for early_late = 1:2 
                sig_channels_id{fidx,aidx,early_late}=[];
            end
        end
    end
    
    for fidx = f2plot
    
    %%% LOAD DATA

        info_mean_sig.m1 = cell(1,2);
        info_mean_sig.dls = cell(1,2);
        raw_LFP_traces.m1 = cell(1,2);
        raw_LFP_traces.m1{1}=cell(1,1); 
        raw_LFP_traces.m1{2}=cell(1,1);
        raw_LFP_traces.dls = cell(1,2);        
        raw_LFP_traces.dls{1}=cell(1,1); 
        raw_LFP_traces.dls{2}=cell(1,1);

        for animal = 1:numel(params.animals)
            for aidx = 1:length(params.areas)
                for early_late = 1:2
                    sig_chan = []; all_chan = [];
                    for chan = 1:size(MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info,1)
%                         if sig_chans{animal,n_lfp,fidx}{aidx,early_late}(chan)==1      
                        infQuant = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).info(chan,:);
                        infQuantSh = MIout{animal}.LFP_early_late{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx}).(params.reachFeatures{fidx}).infoSh(chan,:,:);
                        tmp_sig = clusterStat_v2(infQuant,infQuantSh,clusterParams(1),clusterParams(2));
                        infQuant = infQuant-mean(squeeze(infQuantSh),2)';
                        infQuant(~tmp_sig) = 0;
                        [~,i] = max(infQuant);
                        info_mean_sig.(areas_lab{aidx}){early_late} = [info_mean_sig.(areas_lab{aidx}){early_late}; infQuant];
                        raw_LFP_traces.(areas_lab{aidx}){early_late}{end+1} =  zscore(squeeze(LFP_early_late{animal}{early_late}.(params.lfpFeatures{n_lfp}).(params.areas{aidx})(chan,:,:)))';
                        
                        % check significace of information for the channel
                        
                        sig_chan = [sig_chan i>=sigTime_start & i<=sigTime_end];
                        all_chan = [all_chan 1];
                    end
                    sig_channels(animal,n_lfp,fidx,aidx,early_late) = sum(sig_chan);
                    all_channels(animal,n_lfp,fidx,aidx,early_late) = sum(all_chan);
                    sig_channels_id{fidx,aidx,early_late} = cat(2,sig_channels_id{fidx,aidx,early_late},sig_chan);
%                         end
                end
            end
        end

        timeBins4LFP = round(newBinTimes/10)+150;
        
        % Fill activity and info structs
        for aidx = 1:2
            aLab = areas_lab{aidx};
            for elidx = 1:2
                elLab = early_late_lab{elidx};
                info.(aLab).(elLab) = [];
                info_all.(aLab).(elLab) = [];
                info_LFP.(aLab).(elLab) = [];
                info_LFP_z.(aLab).(elLab) = [];
                info_LFP_STD.(aLab).(elLab) = [];
                raw_LFP.(aLab).(elLab) = [];
                raw_LFP_z.(aLab).(elLab) = [];
                raw_LFP_STD.(aLab).(elLab) = [];
                for unit = 1:size(info_mean_sig.(aLab){elidx},1)
                    tmp_info = info_mean_sig.(aLab){elidx}(unit,:);
                    tmp_LFP = [];
                    for bins = 1:size(timeBins4LFP,1)
                        tmp_LFP = [tmp_LFP mean(raw_LFP_traces.(aLab){elidx}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2)),'all')];
                    end
                    tmp_LFP_STD = [];
                    for bins = 1:size(timeBins4LFP,1)
                        trial_sum = sum(raw_LFP_traces.(aLab){elidx}{unit+1}(:,timeBins4LFP(bins,1):timeBins4LFP(bins,2))');
                        tmp_LFP_STD = [tmp_LFP_STD std(trial_sum)];
                    end
                    info_all.(aLab).(elLab) = [info_all.(aLab).(elLab); tmp_info];
                    [v,i] = max(tmp_info);
                    if v>0 && i>=sigTime_start && i<=sigTime_end
                        info.(aLab).(elLab) = [info.(aLab).(elLab); tmp_info];
                        info_LFP.(aLab).(elLab) = [info_LFP.(aLab).(elLab); tmp_LFP];
                        info_LFP_z.(aLab).(elLab) = [info_LFP_z.(aLab).(elLab); zscore(tmp_LFP)];
                        info_LFP_STD.(aLab).(elLab) = [info_LFP_STD.(aLab).(elLab); tmp_LFP_STD];

                        % We append raw profiles for both MI sig and non-sig
                        % channels
                        raw_LFP.(aLab).(elLab) = [raw_LFP.(aLab).(elLab); tmp_LFP];
                        raw_LFP_z.(aLab).(elLab) = [raw_LFP_z.(aLab).(elLab); zscore(tmp_LFP)];
                    else
                        raw_LFP.(aLab).(elLab) = [raw_LFP.(aLab).(elLab); tmp_LFP];
                        raw_LFP_z.(aLab).(elLab) = [raw_LFP_z.(aLab).(elLab); zscore(tmp_LFP)];
                        raw_LFP_STD.(aLab).(elLab) = [raw_LFP_STD.(aLab).(elLab); tmp_LFP_STD];
                    end
                end
            end
        end
        

    disp(['LFP | ' params.reachFeatures{fidx} ' | M1 | early: ' num2str(size(info.m1.early,1)) ' out of ' num2str(size(info_all.m1.early,1)) ' - ' num2str(100*size(info.m1.early,1)/size(info_all.m1.early,1))  ...
                                                  '% | late: ' num2str(size(info.m1.late,1)) ' out of ' num2str(size(info_all.m1.late,1)) ' - ' num2str(100*size(info.m1.late,1)/size(info_all.m1.late,1))  '%'])
    disp(['LFP | ' params.reachFeatures{fidx} ' | DLS | early: ' num2str(size(info.dls.early,1)) ' out of ' num2str(size(info_all.dls.early,1)) ' - ' num2str(100*size(info.dls.early,1)/size(info_all.dls.early,1))  ...
                                                  '% | late: ' num2str(size(info.dls.late,1)) ' out of ' num2str(size(info_all.dls.late,1)) ' - ' num2str(100*size(info.dls.late,1)/size(info_all.dls.late,1))  '%'])
    %%% PLOT RAW MEAN / SEM
    
        figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(raw_LFP_z.m1.early); abs(raw_LFP_z.m1.late)];
                c = [ones(1,size(raw_LFP_z.m1.early,1)) 2*ones(1,size(raw_LFP_z.m1.late,1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 2]);           
                g(1,1).set_title('m1 lfp'); 

                y = [abs(raw_LFP_z.dls.early); abs(raw_LFP_z.dls.late)];
                c = [ones(1,size(raw_LFP_z.dls.early,1)) 2*ones(1,size(raw_LFP_z.dls.late,1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 2]);     
                g(1,2).set_title('dls lfp');   
                
                y = [abs(raw_LFP_z.m1.early); abs(raw_LFP_z.dls.early)];
                c = [ones(1,size(raw_LFP_z.m1.early,1)) 2*ones(1,size(raw_LFP_z.dls.early,1))];
                g(2,1)=gramm('x',x,'y',y,'color',c);
                g(2,1).stat_summary('type','sem');
                g(2,1).axe_property('XLim',[12 86]);
                g(2,1).axe_property('YLim',[0 2]);           
                g(2,1).set_title('early lfp');    
                
                y = [abs(raw_LFP_z.m1.late); abs(raw_LFP_z.dls.late)];
                c = [ones(1,size(raw_LFP_z.m1.late,1)) 2*ones(1,size(raw_LFP_z.dls.late,1))];
                g(2,2)=gramm('x',x,'y',y,'color',c);
                g(2,2).stat_summary('type','sem');
                g(2,2).axe_property('XLim',[12 86]);
                g(2,2).axe_property('YLim',[0 2]);     
                g(2,2).set_title('late lfp');  
                
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\LFP_raw_' params.reachFeatures{fidx} ' ' date '_mean_sem_2.svg'],'-painters','-dsvg');
                end   

    %%% PEAK RAW SHIFT

        for aidx = 1:2
            aLab = areas_lab{aidx};
            for elidx = 1:2
                elLab = early_late_lab{elidx};
                timing.(aLab).(elLab).std2 = []; timing.(aLab).(elLab).std15 = []; timing.(aLab).(elLab).std3 = []; timing.(aLab).(elLab).com = []; timing.(aLab).(elLab).max = [];  timing.(aLab).(elLab).median = [];
                sel_chans.(aLab).(elLab) = timing.(aLab).(elLab);
                for unit = 1:size((raw_LFP_z.(aLab).(elLab)),1)

                    timing.(aLab).(elLab).std15 = [timing.(aLab).(elLab).std15 min([find(raw_LFP_z.(aLab).(elLab)(unit,:)>1.5) find(raw_LFP_z.(aLab).(elLab)(unit,:)<-1.5)])];   
                    sel_chans.(aLab).(elLab).std15 = [sel_chans.(aLab).(elLab).std15 ~isempty(min([find(raw_LFP_z.(aLab).(elLab)(unit,:)>1.5) find(raw_LFP_z.(aLab).(elLab)(unit,:)<-1.5)]))];
                    timing.(aLab).(elLab).std2 = [timing.(aLab).(elLab).std2 min([find(raw_LFP_z.(aLab).(elLab)(unit,:)>2) find(raw_LFP_z.(aLab).(elLab)(unit,:)<-2)])];            
                    sel_chans.(aLab).(elLab).std2 = [sel_chans.(aLab).(elLab).std2 ~isempty(min([find(raw_LFP_z.(aLab).(elLab)(unit,:)>2) find(raw_LFP_z.(aLab).(elLab)(unit,:)<-2)]))];
                    timing.(aLab).(elLab).std3 = [timing.(aLab).(elLab).std3 min([find(raw_LFP_z.(aLab).(elLab)(unit,:)>3) find(raw_LFP_z.(aLab).(elLab)(unit,:)<-3)])];
                    sel_chans.(aLab).(elLab).std3 = [sel_chans.(aLab).(elLab).std3 ~isempty(min([find(raw_LFP_z.(aLab).(elLab)(unit,:)>3) find(raw_LFP_z.(aLab).(elLab)(unit,:)<-3)]))];

                    timing.(aLab).(elLab).com = [timing.(aLab).(elLab).com center_of_mass_time(abs(raw_LFP_z.(aLab).(elLab)(unit,:)))];            
                    sel_chans.(aLab).(elLab).com = [sel_chans.(aLab).(elLab).com 1];

                    [~,peak_time] = max(abs(raw_LFP_z.(aLab).(elLab)(unit,:)));
                    timing.(aLab).(elLab).max = [timing.(aLab).(elLab).max peak_time];
                    sel_chans.(aLab).(elLab).max = [sel_chans.(aLab).(elLab).max 1];

                    timing.(aLab).(elLab).median = [timing.(aLab).(elLab).median median_time(abs(raw_LFP_z.(aLab).(elLab)(unit,:)))];
                    sel_chans.(aLab).(elLab).median = [sel_chans.(aLab).(elLab).median 1];
                end
            end
        end
        
        % Pick only values in the selected time window
        count = 0;
        for timeFeat = time_features_to_plot
            count = count + 1;
            for aidx = 1:2
                aLab = areas_lab{aidx};
                % Pick only values in the selected time window
                sel_chans.(aLab).early.(timeFeat{1}) = find_sel_chans_idxs(sel_chans.(aLab).early.(timeFeat{1}), timing.(aLab).early.(timeFeat{1}), maxTime, minTime);
                timing.(aLab).early.(timeFeat{1}) = timing.(aLab).early.(timeFeat{1})(timing.(aLab).early.(timeFeat{1})>minTime & timing.(aLab).early.(timeFeat{1})<maxTime);      
                sel_chans.(aLab).late.(timeFeat{1}) = find_sel_chans_idxs(sel_chans.(aLab).late.(timeFeat{1}), timing.(aLab).late.(timeFeat{1}), maxTime, minTime);
                timing.(aLab).late.(timeFeat{1}) = timing.(aLab).late.(timeFeat{1})(timing.(aLab).late.(timeFeat{1})>minTime & timing.(aLab).late.(timeFeat{1})<maxTime);     
                
                % find paired idxs of channels selected both early and late
                tmp_early_idxs = 1:size(raw_LFP_z.(aLab).early,1); tmp_late_idxs = 1:size(raw_LFP_z.(aLab).late,1);
                tmp_early_idxs = tmp_early_idxs(sel_chans.(aLab).early.(timeFeat{1})); tmp_late_idxs = tmp_late_idxs(sel_chans.(aLab).late.(timeFeat{1}));
                tmp_overlap = intersect(tmp_early_idxs,tmp_late_idxs); 
                paired_EL_chans_idxs{aidx,1}.(timeFeat{1}) = ismember(tmp_early_idxs,tmp_overlap);
                paired_EL_chans_idxs{aidx,2}.(timeFeat{1}) = ismember(tmp_late_idxs,tmp_overlap);
               
                % Remove non paired significant M1 and DLS channels
                tmp_idxs = find(sel_chans.(aLab).early.(timeFeat{1}));
                tmp_idxs = tmp_idxs(~paired_EL_chans_idxs{aidx,1}.(timeFeat{1}));
                sel_chans.(aLab).early.(timeFeat{1})(tmp_idxs) = 0;
                tmp_idxs = find(sel_chans.(aLab).late.(timeFeat{1}));
                tmp_idxs = tmp_idxs(~paired_EL_chans_idxs{aidx,2}.(timeFeat{1}));
                sel_chans.(aLab).late.(timeFeat{1})(tmp_idxs) = 0;
                
            end
        end
        
        % plot raw timing comparisons (M1 vs DLS)
        count = 0;
        figure('units','normalized','position',[0.15,-0.02,0.79,0.91]);
        clearvars g;
        nTimeFields = numel(time_features_to_plot); p_naive = []; p_skilled = [];
        grouping_M1_DLS_naive = []; ID_M1_DLS_naive = {}; y_M1_DLS_naive = []; title_M1_DLS_naive = 'M1-DLS pairs naive';
        grouping_M1_DLS_skilled = []; ID_M1_DLS_skilled = {}; y_M1_DLS_skilled = []; title_M1_DLS_skilled = 'M1-DLS pairs skilled';
        for timeFeat = time_features_to_plot
            
            count = count + 1;
            
            %%% Plot pooled animals
            % Channels per animal (used in the M1 vs DLS comparison, to
            % olny pick pairs of M1 and DLS channels within animals - with
            % both M1 and DLS channels are carrying significant early and late info)

            % M1 vs DLS (pairs)
            peak_M1_DLS_naive{2*count-1} = []; peak_M1_DLS_naive{2*count} = [];
            peak_M1_DLS_naive{2*count-1} = timing.m1.early.(timeFeat{1}); peak_M1_DLS_naive{2*count} = timing.dls.early.(timeFeat{1});
            
            grouping_M1_DLS_naive = [grouping_M1_DLS_naive ones(1,length(peak_M1_DLS_naive{2*count-1})) 2*ones(1,length(peak_M1_DLS_naive{2*count}))];
            y_M1_DLS_naive = [y_M1_DLS_naive peak_M1_DLS_naive{2*count-1} peak_M1_DLS_naive{2*count}];
            
            startId = numel(ID_M1_DLS_naive);
            ID_M1_DLS_naive(startId+1:startId+length(peak_M1_DLS_naive{2*count-1})) = {['M1 ' timeFeat{1}]};
            ID_M1_DLS_naive(startId+1+length(peak_M1_DLS_naive{2*count-1}):length(grouping_M1_DLS_naive)) = {['DLS ' timeFeat{1}]};
            
            peak_M1_DLS_skilled{2*count-1} = []; peak_M1_DLS_skilled{2*count} = [];
            peak_M1_DLS_skilled{2*count-1} = timing.m1.late.(timeFeat{1}); peak_M1_DLS_skilled{2*count} = timing.dls.late.(timeFeat{1});
            
            grouping_M1_DLS_skilled = [grouping_M1_DLS_skilled ones(1,length(peak_M1_DLS_skilled{2*count-1})) 2*ones(1,length(peak_M1_DLS_skilled{2*count}))];
            y_M1_DLS_skilled = [y_M1_DLS_skilled peak_M1_DLS_skilled{2*count-1} peak_M1_DLS_skilled{2*count}];
            
            startId = numel(ID_M1_DLS_skilled);
            ID_M1_DLS_skilled(startId+1:startId+length(peak_M1_DLS_skilled{2*count-1})) = {['M1 ' timeFeat{1}]};
            ID_M1_DLS_skilled(startId+1+length(peak_M1_DLS_skilled{2*count-1}):length(grouping_M1_DLS_skilled)) = {['DLS ' timeFeat{1}]};
            
            % paired test
            [p_naive(count), ~, ~] = ranksum(peak_M1_DLS_naive{2*count-1},peak_M1_DLS_naive{2*count});
            title_M1_DLS_naive = [title_M1_DLS_naive ', rpp_{' timeFeat{1} '}=',num2str(p_naive(count),2)];
            [p_skilled(count), ~, ~] = ranksum(peak_M1_DLS_skilled{2*count-1},peak_M1_DLS_skilled{2*count});
            title_M1_DLS_skilled = [title_M1_DLS_skilled ', rpp_{' timeFeat{1} '}=',num2str(p_skilled(count),2)];
        end

        g(1,1)=gramm('x',ID_M1_DLS_naive,'y',y_M1_DLS_naive,'color',grouping_M1_DLS_naive);
        g(1,1).stat_boxplot();
        g(1,1).set_order_options('x',0);
        g(1,1).coord_flip();
        g(1,1).axe_property('YLim',[11 86]);
        g(1,1).axe_property('YTick',[11:25:86]);
        g(1,1).axe_property('YTickLabels',[-1:0.5:0.5]);
        g(1,1).set_title(title_M1_DLS_naive)
        
        g(2,1)=gramm('x',ID_M1_DLS_skilled,'y',y_M1_DLS_skilled,'color',grouping_M1_DLS_skilled);
        g(2,1).stat_boxplot();
        g(2,1).set_order_options('x',0);
        g(2,1).coord_flip();
        g(2,1).axe_property('YLim',[11 86]);
        g(2,1).axe_property('YTick',[11:25:86]);
        g(2,1).axe_property('YTickLabels',[-1:0.5:0.5]);
        g(2,1).set_title(title_M1_DLS_skilled)
        
        g.set_title(['LFP chan activity timing, '  params.reachFeatures{fidx}])
        g.draw();

        if plot_params.save == 1
            print([plot_params.save_path '\LFP_raw_' params.reachFeatures{fidx} ' ' date '_timingShift.svg'],'-painters','-dsvg');
            print([plot_params.save_path '\LFP_raw_' params.reachFeatures{fidx} ' ' date '_timingShift.png'],'-dpng');
        end
        
        % plot raw timing comparisons (differences)
        count = 0;
        figure('units','normalized','position',[0.15,-0.02,0.79,0.91]);
        clearvars g;
        nTimeFields = numel(time_features_to_plot); p_EL = []; p_areas = []; stat_M1_DLS = [];
        grouping_EL = []; ID_EL = {}; y_EL = []; title_EL = 'Naive-Skilled channels';
        grouping_M1_DLS = []; ID_M1_DLS = {}; y_M1_DLS = []; title_M1_DLS = 'M1-DLS pairs';
        for timeFeat = time_features_to_plot
            
            count = count + 1;
            
            %%% Plot pooled animals
            % Channels per animal (used in the M1 vs DLS comparison, to
            % olny pick pairs of M1 and DLS channels within animals - with
            % both M1 and DLS channels are carrying significant early and late info)

            chans_per_animal.m1.early = sum_per_sess(sel_chans.m1.early.(timeFeat{1}),squeeze(all_channels(:,1,fidx,1,1))');
            chans_per_animal.m1.late = sum_per_sess(sel_chans.m1.late.(timeFeat{1}),squeeze(all_channels(:,1,fidx,1,2))');
            chans_per_animal.dls.early = sum_per_sess(sel_chans.dls.early.(timeFeat{1}),squeeze(all_channels(:,1,fidx,2,1))');
            chans_per_animal.dls.late = sum_per_sess(sel_chans.dls.late.(timeFeat{1}),squeeze(all_channels(:,1,fidx,2,2))');
            start_chan.m1.early = [1,cumsum(chans_per_animal.m1.early(1:end-1))+1]; end_chan.m1.early = cumsum(chans_per_animal.m1.early);
            start_chan.m1.late = [1,cumsum(chans_per_animal.m1.late(1:end-1))+1]; end_chan.m1.late = cumsum(chans_per_animal.m1.late);
            start_chan.dls.early = [1,cumsum(chans_per_animal.dls.early(1:end-1))+1]; end_chan.dls.early = cumsum(chans_per_animal.dls.early);
            start_chan.dls.late = [1,cumsum(chans_per_animal.dls.late(1:end-1))+1]; end_chan.dls.late = cumsum(chans_per_animal.dls.late);
            
            % Naive vs skilled
            peak_EL{2*count-1} = timing.m1.early.(timeFeat{1})(paired_EL_chans_idxs{1,1}.(timeFeat{1}))-timing.m1.late.(timeFeat{1})(paired_EL_chans_idxs{1,2}.(timeFeat{1}));
            peak_EL{2*count} = timing.dls.early.(timeFeat{1})(paired_EL_chans_idxs{2,1}.(timeFeat{1}))-timing.dls.late.(timeFeat{1})(paired_EL_chans_idxs{2,2}.(timeFeat{1}));
            grouping_EL = [grouping_EL ones(1,length(peak_EL{2*count-1})) 2*ones(1,length(peak_EL{2*count}))];
            y_EL = [y_EL peak_EL{2*count-1} peak_EL{2*count}];
            
            startId = numel(ID_EL);
            ID_EL(startId+1:startId+length(peak_EL{2*count-1})) = {['M1 ' timeFeat{1}]};
            ID_EL(startId+1+length(peak_EL{2*count-1}):length(grouping_EL)) = {['DLS ' timeFeat{1}]};
            
            % sig test
            [p, ~, ~] = ranksum(peak_EL{2*count-1},peak_EL{2*count});
            title_EL = [title_EL ', rp_{' timeFeat{1} '}=',num2str(p,2)];
            [p_EL(2*count-1), ~, stat_EL(2*count-1)] = signrank(peak_EL{2*count-1}, zeros(1,numel(peak_EL{2*count-1})), 'method' , 'approximate');
            [p_EL(2*count), ~, stat_EL(2*count)] = signrank(peak_EL{2*count}, zeros(1,numel(peak_EL{2*count})), 'method' , 'approximate');
            
            % M1 vs DLS (pairs)
            peak_M1_DLS{2*count-1} = []; peak_M1_DLS{2*count} = [];
            for animal = 1:numel(chans_per_animal.m1.early)
                tmp_peak1_animal = [timing.m1.early.(timeFeat{1})(start_chan.m1.early(animal):end_chan.m1.early(animal))...
                    - timing.dls.early.(timeFeat{1})(start_chan.dls.early(animal):end_chan.dls.early(animal))'];
                peak_M1_DLS{2*count-1} = [peak_M1_DLS{2*count-1}, tmp_peak1_animal(:)'];
                tmp_peak2_animal = [timing.m1.late.(timeFeat{1})(start_chan.m1.late(animal):end_chan.m1.late(animal))...
                    - timing.dls.late.(timeFeat{1})(start_chan.dls.late(animal):end_chan.dls.late(animal))'];
                peak_M1_DLS{2*count} = [peak_M1_DLS{2*count}, tmp_peak2_animal(:)'];
            end
            grouping_M1_DLS = [grouping_M1_DLS ones(1,length(peak_M1_DLS{2*count-1})) 2*ones(1,length(peak_M1_DLS{2*count}))];
            y_M1_DLS = [y_M1_DLS peak_M1_DLS{2*count-1} peak_M1_DLS{2*count}];
            
            startId = numel(ID_M1_DLS);
            ID_M1_DLS(startId+1:startId+length(peak_M1_DLS{2*count-1})) = {['Naive ' timeFeat{1}]};
            ID_M1_DLS(startId+1+length(peak_M1_DLS{2*count-1}):length(grouping_M1_DLS)) = {['Skilled ' timeFeat{1}]};
            
            % paired test
            [p2, ~, ~] = signrank(peak_M1_DLS{2*count-1},peak_M1_DLS{2*count}, 'method', 'approximate');
            title_M1_DLS = [title_M1_DLS ', rpp_{' timeFeat{1} '}=',num2str(p2,2), ' (' num2str(numel(peak_M1_DLS{2*count-1})) ') pairs'];
            [p_areas(2*count-1), ~, tmp_stat] = signrank(peak_M1_DLS{2*count-1}, zeros(1,numel(peak_M1_DLS{2*count-1})), 'method', 'approximate');
            stat_M1_DLS(2*count-1) = tmp_stat.zval;
            [p_areas(2*count), ~, tmp_stat] = signrank(peak_M1_DLS{2*count}, zeros(1,numel(peak_M1_DLS{2*count})), 'method', 'approximate');
            stat_M1_DLS(2*count) = tmp_stat.zval;
        end

        g(1,1)=gramm('x',ID_EL,'y',y_EL,'color',grouping_EL);
        g(1,1).stat_boxplot();
        g(1,1).set_order_options('x',0);
%         g(count,1).coord_flip();
        g(1,1).axe_property('YLim',[-50 50]);
        g(1,1).set_title(title_EL)
%         g(1,1).geom_label(1:count,45,num2str(p_EL));

        
        g(2,1)=gramm('x',ID_M1_DLS,'y',y_M1_DLS,'color',grouping_M1_DLS);
        g(2,1).stat_boxplot();
        g(2,1).set_order_options('x',0);
%         g(2,1).coord_flip();
        g(2,1).axe_property('YLim',[-50 50]);
        g(2,1).set_title(title_M1_DLS)
        
        g.set_title(['LFP chan pairs activity timing differences, '  params.reachFeatures{fidx}])
        g.draw();
        p_EL = arrayfun(@(z) num2str(z,2), p_EL, 'UniformOutput', 0);
        p_areas = arrayfun(@(z) num2str(z,2), p_areas, 'UniformOutput', 0);
        stat_M1_DLS = arrayfun(@(z) num2str(z,2), stat_M1_DLS, 'UniformOutput', 0);
    
        text(1:2*count,45*ones(1,2*count),p_EL,'Parent',g(1,1).facet_axes_handles(1),'FontName','Courier');
        yline(0,'k--','Parent',g(1,1).facet_axes_handles(1));
        text(1:2*count,45*ones(1,2*count),p_areas,'Parent',g(2,1).facet_axes_handles(1),'FontName','Courier');
        yline(0,'k--','Parent',g(2,1).facet_axes_handles(1));
        text(1:2*count,40*ones(1,2*count),stat_M1_DLS,'Parent',g(2,1).facet_axes_handles(1),'FontName','Courier');
        yline(0,'k--','Parent',g(2,1).facet_axes_handles(1));

        if plot_params.save == 1
            print([plot_params.save_path '\LFPpairs_raw_diffs_' params.reachFeatures{fidx} ' ' date '_timingShift.svg'],'-painters','-dsvg');
            print([plot_params.save_path '\LFPpairs_raw_diffs_' params.reachFeatures{fidx} ' ' date '_timingShift.png'],'-dpng');
        end
        
        %%% PLOT RAW MEAN / SEM paired cells
    
        figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(raw_LFP_z.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:)); abs(raw_LFP_z.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:))];
                c = [ones(1,size(raw_LFP_z.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:),1)) 2*ones(1,size(raw_LFP_z.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:),1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 2]);           
                g(1,1).set_title('m1 lfp'); 

                y = [abs(raw_LFP_z.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:)); abs(raw_LFP_z.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:))];
                c = [ones(1,size(raw_LFP_z.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:),1)) 2*ones(1,size(raw_LFP_z.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:),1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 2]);     
                g(1,2).set_title('dls lfp');   
                
                y = [abs(raw_LFP_z.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:)); abs(raw_LFP_z.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:))];
                c = [ones(1,size(raw_LFP_z.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:),1)) 2*ones(1,size(raw_LFP_z.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:),1))];
                g(2,1)=gramm('x',x,'y',y,'color',c);
                g(2,1).stat_summary('type','sem');
                g(2,1).axe_property('XLim',[12 86]);
                g(2,1).axe_property('YLim',[0 2]);           
                g(2,1).set_title('early lfp');    
                
                y = [abs(raw_LFP_z.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:)); abs(raw_LFP_z.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:))];
                c = [ones(1,size(raw_LFP_z.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:),1)) 2*ones(1,size(raw_LFP_z.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:),1))];
                g(2,2)=gramm('x',x,'y',y,'color',c);
                g(2,2).stat_summary('type','sem');
                g(2,2).axe_property('XLim',[12 86]);
                g(2,2).axe_property('YLim',[0 2]);     
                g(2,2).set_title('late lfp');  
                
                g.set_title('raw activity paired naive and skilled')
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\LFP_raw_paired_' params.reachFeatures{fidx} ' ' date '_mean_sem_2.svg'],'-painters','-dsvg');
                end   

        
    %%% PLOT INFO CASCADES AND MEAN / SEM
        figure('units','normalized','outerposition',[0 0 1 1])
        clearvars g
            x = [1:98];
            y = [abs(info.m1.early); abs(info.m1.late)];
            c = [ones(1,size(info.m1.early,1)) 2*ones(1,size(info.m1.late,1))];
            g(1,1)=gramm('x',x,'y',y,'color',c);
            g(1,1).stat_summary('type','sem');
            g(1,1).axe_property('XLim',[12 86]);
            g(1,1).axe_property('YLim',[0 0.05]);           
            g(1,1).set_title('m1 lfp info');        
            y = [abs(info.dls.early); abs(info.dls.late)];
            c = [ones(1,size(info.dls.early,1)) 2*ones(1,size(info.dls.late,1))];
            g(1,2)=gramm('x',x,'y',y,'color',c);
            g(1,2).stat_summary('type','sem');
            g(1,2).axe_property('XLim',[12 86]);
            g(1,2).axe_property('YLim',[0 0.05]);     
            g(1,2).set_title('dls lfp info');       

            y = [abs(info.m1.early); abs(info.dls.early)];
            c = [ones(1,size(info.m1.early,1)) 2*ones(1,size(info.dls.early,1))];
            g(2,1)=gramm('x',x,'y',y,'color',c);
            g(2,1).stat_summary('type','sem');
            g(2,1).axe_property('XLim',[12 86]);
            g(2,1).axe_property('YLim',[0 0.05]);           
            g(2,1).set_title('early lfp info');        
            y = [abs(info.m1.late); abs(info.dls.late)];
            c = [ones(1,size(info.m1.late,1)) 2*ones(1,size(info.dls.late,1))];
            g(2,2)=gramm('x',x,'y',y,'color',c);
            g(2,2).stat_summary('type','sem');
            g(2,2).axe_property('XLim',[12 86]);
            g(2,2).axe_property('YLim',[0 0.05]);     
            g(2,2).set_title('late lfp info');       

            g.draw();
            if plot_params.save == 1
                print([plot_params.save_path '\LFP_information_' params.reachFeatures{fidx} ' ' date '_mean_sem.svg'],'-painters','-dsvg');
                print([plot_params.save_path '\LFP_information_' params.reachFeatures{fidx} ' ' date '_mean_sem.png'],'-dpng');
            end 
                
    %%% PEAK SHIFT INFO

        for aidx = 1:2
            aLab = areas_lab{aidx};
            for elidx = 1:2
                elLab = early_late_lab{elidx};
                timing.(aLab).(elLab).std2 = []; timing.(aLab).(elLab).std15 = []; timing.(aLab).(elLab).std3 = []; timing.(aLab).(elLab).com = []; timing.(aLab).(elLab).max = [];  timing.(aLab).(elLab).median = [];
                sel_chans.(aLab).(elLab) = timing.(aLab).(elLab);
                for unit = 1:size((info.(aLab).(elLab)),1)

                    base2SD = 2*std(info.(aLab).(elLab)(unit,:));
                    base15SD = 2*std(info.(aLab).(elLab)(unit,:));
                    base3SD = 3*std(info.(aLab).(elLab)(unit,:));
                    baseMean = mean(info.(aLab).(elLab)(unit,:));
                    timing.(aLab).(elLab).std15 = [timing.(aLab).(elLab).std15 min(find(info.(aLab).(elLab)(unit,:)>baseMean+base15SD))];
                    timing.(aLab).(elLab).std2 = [timing.(aLab).(elLab).std2 min(find(info.(aLab).(elLab)(unit,:)>baseMean+base2SD))];
                    timing.(aLab).(elLab).std3 = [timing.(aLab).(elLab).std3 min(find(info.(aLab).(elLab)(unit,:)>baseMean+base3SD))];
                    sel_chans.(aLab).(elLab).std15 = [sel_chans.(aLab).(elLab).std15 ~isempty(min(find(info.(aLab).(elLab)(unit,:)>baseMean+base15SD)))];
                    sel_chans.(aLab).(elLab).std2 = [sel_chans.(aLab).(elLab).std2 ~isempty(min(find(info.(aLab).(elLab)(unit,:)>baseMean+base2SD)))];
                    sel_chans.(aLab).(elLab).std3 = [sel_chans.(aLab).(elLab).std3 ~isempty(min(find(info.(aLab).(elLab)(unit,:)>baseMean+base3SD)))];

                    timing.(aLab).(elLab).com = [timing.(aLab).(elLab).com center_of_mass_time(abs(info.(aLab).(elLab)(unit,:)))];            
                    sel_chans.(aLab).(elLab).com = [sel_chans.(aLab).(elLab).com 1];

                    [~,peak_time] = max(abs(info.(aLab).(elLab)(unit,:)));
                    timing.(aLab).(elLab).max = [timing.(aLab).(elLab).max peak_time];
                    sel_chans.(aLab).(elLab).max = [sel_chans.(aLab).(elLab).max 1];

                    timing.(aLab).(elLab).median = [timing.(aLab).(elLab).median median_time(abs(info.(aLab).(elLab)(unit,:)))];
                    sel_chans.(aLab).(elLab).median = [sel_chans.(aLab).(elLab).median 1];
                end
            end
        end
        
        % Pick only values in the selected time window
        count = 0;
        for timeFeat = time_features_to_plot
            count = count + 1
            for aidx = 1:2
                aLab = areas_lab{aidx};
                % Pick only values in the selected time window
                sel_chans.(aLab).early.(timeFeat{1}) = find_sel_chans_idxs(sel_chans.(aLab).early.(timeFeat{1}), timing.(aLab).early.(timeFeat{1}), maxTime, minTime);
                timing.(aLab).early.(timeFeat{1}) = timing.(aLab).early.(timeFeat{1})(timing.(aLab).early.(timeFeat{1})>minTime & timing.(aLab).early.(timeFeat{1})<maxTime);      
                sel_chans.(aLab).late.(timeFeat{1}) = find_sel_chans_idxs(sel_chans.(aLab).late.(timeFeat{1}), timing.(aLab).late.(timeFeat{1}), maxTime, minTime);
                timing.(aLab).late.(timeFeat{1}) = timing.(aLab).late.(timeFeat{1})(timing.(aLab).late.(timeFeat{1})>minTime & timing.(aLab).late.(timeFeat{1})<maxTime);     
                
                % find paired idxs of channels selected both early and late
                tmp_early_idxs = find(sig_channels_id{fidx,aidx,1}); tmp_late_idxs = find(sig_channels_id{fidx,aidx,2});
                tmp_early_idxs = tmp_early_idxs(sel_chans.(aLab).early.(timeFeat{1})); tmp_late_idxs = tmp_late_idxs(sel_chans.(aLab).late.(timeFeat{1}));
                tmp_overlap = intersect(tmp_early_idxs,tmp_late_idxs); 
                paired_EL_chans_idxs{aidx,1}.(timeFeat{1}) = ismember(tmp_early_idxs,tmp_overlap);
                paired_EL_chans_idxs{aidx,2}.(timeFeat{1}) = ismember(tmp_late_idxs,tmp_overlap);
                
                % Remove non paired significant M1 and DLS channels
                tmp_idxs = find(sel_chans.(aLab).early.(timeFeat{1}));
                tmp_idxs = tmp_idxs(~paired_EL_chans_idxs{aidx,1}.(timeFeat{1}));
                sel_chans.(aLab).early.(timeFeat{1})(tmp_idxs) = 0;
                tmp_idxs = find(sel_chans.(aLab).late.(timeFeat{1}));
                tmp_idxs = tmp_idxs(~paired_EL_chans_idxs{aidx,2}.(timeFeat{1}));
                sel_chans.(aLab).late.(timeFeat{1})(tmp_idxs) = 0;
                
            end
        end
        
        % plot info timing comparisons (M1 vs DLS)
        count = 0;
        figure('units','normalized','position',[0.15,-0.02,0.79,0.91]);
        clearvars g;
        nTimeFields = numel(time_features_to_plot); p_naive = []; p_skilled = [];
        grouping_M1_DLS_naive = []; ID_M1_DLS_naive = {}; y_M1_DLS_naive = []; title_M1_DLS_naive = 'M1-DLS pairs naive';
        grouping_M1_DLS_skilled = []; ID_M1_DLS_skilled = {}; y_M1_DLS_skilled = []; title_M1_DLS_skilled = 'M1-DLS pairs skilled';
        for timeFeat = time_features_to_plot
            
            count = count + 1;
            
            %%% Plot pooled animals
            % Channels per animal (used in the M1 vs DLS comparison, to
            % olny pick pairs of M1 and DLS channels within animals - with
            % both M1 and DLS channels are carrying significant early and late info)

            % M1 vs DLS (pairs)
            peak_M1_DLS_naive{2*count-1} = []; peak_M1_DLS_naive{2*count} = [];
            peak_M1_DLS_naive{2*count-1} = timing.m1.early.(timeFeat{1}); peak_M1_DLS_naive{2*count} = timing.dls.early.(timeFeat{1});
            
            grouping_M1_DLS_naive = [grouping_M1_DLS_naive ones(1,length(peak_M1_DLS_naive{2*count-1})) 2*ones(1,length(peak_M1_DLS_naive{2*count}))];
            y_M1_DLS_naive = [y_M1_DLS_naive peak_M1_DLS_naive{2*count-1} peak_M1_DLS_naive{2*count}];
            
            startId = numel(ID_M1_DLS_naive);
            ID_M1_DLS_naive(startId+1:startId+length(peak_M1_DLS_naive{2*count-1})) = {['M1 ' timeFeat{1}]};
            ID_M1_DLS_naive(startId+1+length(peak_M1_DLS_naive{2*count-1}):length(grouping_M1_DLS_naive)) = {['DLS ' timeFeat{1}]};
            
            peak_M1_DLS_skilled{2*count-1} = []; peak_M1_DLS_skilled{2*count} = [];
            peak_M1_DLS_skilled{2*count-1} = timing.m1.late.(timeFeat{1}); peak_M1_DLS_skilled{2*count} = timing.dls.late.(timeFeat{1});
            
            grouping_M1_DLS_skilled = [grouping_M1_DLS_skilled ones(1,length(peak_M1_DLS_skilled{2*count-1})) 2*ones(1,length(peak_M1_DLS_skilled{2*count}))];
            y_M1_DLS_skilled = [y_M1_DLS_skilled peak_M1_DLS_skilled{2*count-1} peak_M1_DLS_skilled{2*count}];
            
            startId = numel(ID_M1_DLS_skilled);
            ID_M1_DLS_skilled(startId+1:startId+length(peak_M1_DLS_skilled{2*count-1})) = {['M1 ' timeFeat{1}]};
            ID_M1_DLS_skilled(startId+1+length(peak_M1_DLS_skilled{2*count-1}):length(grouping_M1_DLS_skilled)) = {['DLS ' timeFeat{1}]};
            
            % paired test
            [p_naive(count), ~, ~] = ranksum(peak_M1_DLS_naive{2*count-1},peak_M1_DLS_naive{2*count});
            title_M1_DLS_naive = [title_M1_DLS_naive ', rpp_{' timeFeat{1} '}=',num2str(p_naive(count),2)];
            [p_skilled(count), ~, ~] = ranksum(peak_M1_DLS_skilled{2*count-1},peak_M1_DLS_skilled{2*count});
            title_M1_DLS_skilled = [title_M1_DLS_skilled ', rpp_{' timeFeat{1} '}=',num2str(p_skilled(count),2)];
        end

        g(1,1)=gramm('x',ID_M1_DLS_naive,'y',y_M1_DLS_naive,'color',grouping_M1_DLS_naive);
        g(1,1).stat_boxplot();
        g(1,1).set_order_options('x',0);
        g(1,1).coord_flip();
        g(1,1).axe_property('YLim',[11 86]);
        g(1,1).axe_property('YTick',[11:25:86]);
        g(1,1).axe_property('YTickLabels',[-1:0.5:0.5]);
        g(1,1).set_title(title_M1_DLS_naive)
        
        g(2,1)=gramm('x',ID_M1_DLS_skilled,'y',y_M1_DLS_skilled,'color',grouping_M1_DLS_skilled);
        g(2,1).stat_boxplot();
        g(2,1).set_order_options('x',0);
        g(2,1).coord_flip();
        g(2,1).axe_property('YLim',[11 86]);
        g(2,1).axe_property('YTick',[11:25:86]);
        g(2,1).axe_property('YTickLabels',[-1:0.5:0.5]);
        g(2,1).set_title(title_M1_DLS_skilled)
        
        g.set_title(['LFP chan info timing, '  params.reachFeatures{fidx}])
        g.draw();

        if plot_params.save == 1
            print([plot_params.save_path '\LFP_info_' params.reachFeatures{fidx} ' ' date '_timingShift.svg'],'-painters','-dsvg');
            print([plot_params.save_path '\LFP_info_' params.reachFeatures{fidx} ' ' date '_timingShift.png'],'-dpng');
        end
        
        
        % plot info timing comparisons (differences)
        count = 0;
        figure('units','normalized','position',[0.15,-0.02,0.79,0.91]);
        clearvars g;
        nTimeFields = numel(time_features_to_plot); p_EL = []; p_areas = []; stat_M1_DLS = [];
        grouping_EL = []; ID_EL = {}; y_EL = []; title_EL = 'Naive-Skilled channels';
        grouping_M1_DLS = []; ID_M1_DLS = {}; y_M1_DLS = []; title_M1_DLS = 'M1-DLS pairs';
        for timeFeat = time_features_to_plot
            
            count = count + 1;
            
            %%% Plot pooled animals
            % Channels per animal (used in the M1 vs DLS comparison, to
            % olny pick pairs of M1 and DLS channels within animals - with
            % both M1 and DLS channels are carrying significant early and late info)

            chans_per_animal.m1.early = sum_per_sess(sel_chans.m1.early.(timeFeat{1}),squeeze(sig_channels(:,1,fidx,1,1))');
            chans_per_animal.m1.late = sum_per_sess(sel_chans.m1.late.(timeFeat{1}),squeeze(sig_channels(:,1,fidx,1,2))');
            chans_per_animal.dls.early = sum_per_sess(sel_chans.dls.early.(timeFeat{1}),squeeze(sig_channels(:,1,fidx,2,1))');
            chans_per_animal.dls.late = sum_per_sess(sel_chans.dls.late.(timeFeat{1}),squeeze(sig_channels(:,1,fidx,2,2))');
            start_chan.m1.early = [1,cumsum(chans_per_animal.m1.early(1:end-1))+1]; end_chan.m1.early = cumsum(chans_per_animal.m1.early);
            start_chan.m1.late = [1,cumsum(chans_per_animal.m1.late(1:end-1))+1]; end_chan.m1.late = cumsum(chans_per_animal.m1.late);
            start_chan.dls.early = [1,cumsum(chans_per_animal.dls.early(1:end-1))+1]; end_chan.dls.early = cumsum(chans_per_animal.dls.early);
            start_chan.dls.late = [1,cumsum(chans_per_animal.dls.late(1:end-1))+1]; end_chan.dls.late = cumsum(chans_per_animal.dls.late);
            
            % Naive vs skilled
            
            peak_EL{2*count-1} = timing.m1.early.(timeFeat{1})(paired_EL_chans_idxs{1,1}.(timeFeat{1}))-timing.m1.late.(timeFeat{1})(paired_EL_chans_idxs{1,2}.(timeFeat{1}));
            peak_EL{2*count} = timing.dls.early.(timeFeat{1})(paired_EL_chans_idxs{2,1}.(timeFeat{1}))-timing.dls.late.(timeFeat{1})(paired_EL_chans_idxs{2,2}.(timeFeat{1}));
            grouping_EL = [grouping_EL ones(1,length(peak_EL{2*count-1})) 2*ones(1,length(peak_EL{2*count}))];
            y_EL = [y_EL peak_EL{2*count-1} peak_EL{2*count}];
            
            startId = numel(ID_EL);
            ID_EL(startId+1:startId+length(peak_EL{2*count-1})) = {['M1 ' timeFeat{1}]};
            ID_EL(startId+1+length(peak_EL{2*count-1}):length(grouping_EL)) = {['DLS ' timeFeat{1}]};
            
            % sig test
            [p, ~, ~] = ranksum(peak_EL{2*count-1},peak_EL{2*count});
            title_EL = [title_EL ', rp_{' timeFeat{1} '}=',num2str(p,2)];
            [p_EL(2*count-1), ~, ~] = signrank(peak_EL{2*count-1}, zeros(1,numel(peak_EL{2*count-1})), 'method', 'approximate');
            [p_EL(2*count), ~, ~] = signrank(peak_EL{2*count}, zeros(1,numel(peak_EL{2*count})), 'method', 'approximate');
            
            % M1 vs DLS (pairs)
            peak_M1_DLS{2*count-1} = []; peak_M1_DLS{2*count} = [];
            for animal = 1:numel(chans_per_animal.m1.early)
                tmp_peak1_animal = [timing.m1.early.(timeFeat{1})(start_chan.m1.early(animal):end_chan.m1.early(animal))...
                    - timing.dls.early.(timeFeat{1})(start_chan.dls.early(animal):end_chan.dls.early(animal))'];
                peak_M1_DLS{2*count-1} = [peak_M1_DLS{2*count-1}, tmp_peak1_animal(:)'];
                tmp_peak2_animal = [timing.m1.late.(timeFeat{1})(start_chan.m1.late(animal):end_chan.m1.late(animal))...
                    - timing.dls.late.(timeFeat{1})(start_chan.dls.late(animal):end_chan.dls.late(animal))'];
                peak_M1_DLS{2*count} = [peak_M1_DLS{2*count}, tmp_peak2_animal(:)'];
            end
            grouping_M1_DLS = [grouping_M1_DLS ones(1,length(peak_M1_DLS{2*count-1})) 2*ones(1,length(peak_M1_DLS{2*count}))];
            y_M1_DLS = [y_M1_DLS peak_M1_DLS{2*count-1} peak_M1_DLS{2*count}];
            
            startId = numel(ID_M1_DLS);
            ID_M1_DLS(startId+1:startId+length(peak_M1_DLS{2*count-1})) = {['Naive ' timeFeat{1}]};
            ID_M1_DLS(startId+1+length(peak_M1_DLS{2*count-1}):length(grouping_M1_DLS)) = {['Skilled ' timeFeat{1}]};
            
            % paired test
            [p2, ~, ~] = signrank(peak_M1_DLS{2*count-1},peak_M1_DLS{2*count}, 'method', 'approximate');
            title_M1_DLS = [title_M1_DLS ', rpp_{' timeFeat{1} '}=',num2str(p2,2) ' (' num2str(numel(peak_M1_DLS{2*count-1})) ') pairs'];
            [p_areas(2*count-1), ~, tmp_stat] = signrank(peak_M1_DLS{2*count-1}, zeros(1,numel(peak_M1_DLS{2*count-1})), 'method' , 'approximate');
            stat_M1_DLS(2*count-1) = tmp_stat.zval;
            [p_areas(2*count), ~, tmp_stat] = signrank(peak_M1_DLS{2*count}, zeros(1,numel(peak_M1_DLS{2*count-1})), 'method' , 'approximate');
            stat_M1_DLS(2*count) = tmp_stat.zval;
        end

        g(1,1)=gramm('x',ID_EL,'y',y_EL,'color',grouping_EL);
        g(1,1).stat_boxplot();
        g(1,1).set_order_options('x',0);
%         g(count,1).coord_flip();
        g(1,1).axe_property('YLim',[-50 50]);
        g(1,1).set_title(title_EL)
        
        g(2,1)=gramm('x',ID_M1_DLS,'y',y_M1_DLS,'color',grouping_M1_DLS);
        g(2,1).stat_boxplot();
        g(2,1).set_order_options('x',0);
%         g(2,1).coord_flip();
        g(2,1).axe_property('YLim',[-50 50]);
        g(2,1).set_title(title_M1_DLS)
        
        g.set_title(['LFP chan pairs info timing differences, '  params.reachFeatures{fidx}])
        g.draw();
        p_EL = arrayfun(@(z) num2str(z,2), p_EL, 'UniformOutput', 0);
        p_areas = arrayfun(@(z) num2str(z,2), p_areas, 'UniformOutput', 0);
        stat_M1_DLS = arrayfun(@(z) num2str(z,2), stat_M1_DLS, 'UniformOutput', 0);

        text(1:2*count,45*ones(1,2*count),p_EL,'Parent',g(1,1).facet_axes_handles(1),'FontName','Courier');
        yline(0,'k--','Parent',g(1,1).facet_axes_handles(1));
        text(1:2*count,45*ones(1,2*count),p_areas,'Parent',g(2,1).facet_axes_handles(1),'FontName','Courier');
        yline(0,'k--','Parent',g(2,1).facet_axes_handles(1));
        text(1:2*count,40*ones(1,2*count),stat_M1_DLS,'Parent',g(2,1).facet_axes_handles(1),'FontName','Courier');
        yline(0,'k--','Parent',g(2,1).facet_axes_handles(1));

        if plot_params.save == 1
            print([plot_params.save_path '\LFPpairs_info_' params.reachFeatures{fidx} ' ' date '_timingShift.svg'],'-painters','-dsvg');
            print([plot_params.save_path '\LFPpairs_info_' params.reachFeatures{fidx} ' ' date '_timingShift.png'],'-dpng');
        end
        
        %%% PLOT INFO MEAN / SEM paired cells
        if params.plot_paired_profiles
            figure('units','normalized','outerposition',[0 0 1 1])
            clearvars g
                x = [1:98];
                y = [abs(info.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:)); abs(info.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:))];
                c = [ones(1,size(info.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:),1)) 2*ones(1,size(info.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:),1))];
                g(1,1)=gramm('x',x,'y',y,'color',c);
                g(1,1).stat_summary('type','sem');
                g(1,1).axe_property('XLim',[12 86]);
                g(1,1).axe_property('YLim',[0 0.05]);           
                g(1,1).set_title('m1 lfp'); 

                y = [abs(info.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:)); abs(info.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:))];
                c = [ones(1,size(info.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:),1)) 2*ones(1,size(info.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:),1))];
                g(1,2)=gramm('x',x,'y',y,'color',c);
                g(1,2).stat_summary('type','sem');
                g(1,2).axe_property('XLim',[12 86]);
                g(1,2).axe_property('YLim',[0 0.05]);     
                g(1,2).set_title('dls lfp');   
                
                y = [abs(info.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:)); abs(info.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:))];
                c = [ones(1,size(info.m1.early(paired_EL_chans_idxs{1,1}.(timeFeat{1}),:),1)) 2*ones(1,size(info.dls.early(paired_EL_chans_idxs{2,1}.(timeFeat{1}),:),1))];
                g(2,1)=gramm('x',x,'y',y,'color',c);
                g(2,1).stat_summary('type','sem');
                g(2,1).axe_property('XLim',[12 86]);
                g(2,1).axe_property('YLim',[0 0.05]);           
                g(2,1).set_title('early lfp');    
                
                y = [abs(info.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:)); abs(info.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:))];
                c = [ones(1,size(info.m1.late(paired_EL_chans_idxs{1,2}.(timeFeat{1}),:),1)) 2*ones(1,size(info.dls.late(paired_EL_chans_idxs{2,2}.(timeFeat{1}),:),1))];
                g(2,2)=gramm('x',x,'y',y,'color',c);
                g(2,2).stat_summary('type','sem');
                g(2,2).axe_property('XLim',[12 86]);
                g(2,2).axe_property('YLim',[0 0.05]);     
                g(2,2).set_title('late lfp');  
                
                g.set_title('info profiles (paired naive and skilled channels)')
                g.draw();
                if plot_params.save == 1
                    print([plot_params.save_path '\LFP_info_paired_' params.reachFeatures{fidx} ' ' date '_mean_sem_2.svg'],'-painters','-dsvg');
                end   
        end
    %%% NEW NAIVE SKILLED COMPARISON HISTOGRAM/CDF

        figure;
            subplot(2,2,1); hold on;
                histogram(max(info_all.m1.early'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(info_all.m1.late'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(info_all.m1.early'),max(info_all.m1.late'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,2);
                hold on;
                cdfplot(max(info_all.m1.early'));
                cdfplot(max(info_all.m1.late'));
                xlim([0 0.25]);
                grid off;
            subplot(2,2,2); hold on;
                histogram(max(info_all.dls.early'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                histogram(max(info_all.dls.late'),[0:0.01:.1],'Normalization','Probability','DisplayStyle','Stairs');
                [p,~,stats] = ranksum(max(info_all.dls.early'),max(info_all.dls.late'));
                title(p);
                xlim([0 0.25]);
            subplot(2,2,4);
                hold on;
                cdfplot(max(info_all.dls.early'));
                cdfplot(max(info_all.dls.late'));
                xlim([0 0.25]);
                grid off;
            if plot_params.save == 1
                print([plot_params.save_path '\LFP_information_new_skilled_' params.reachFeatures{fidx} ' ' date '.svg'],'-painters','-dsvg');
                print([plot_params.save_path '\LFP_information_new_skilled_' params.reachFeatures{fidx} ' ' date '.png'],'-dpng');
            end

    end
end
end

