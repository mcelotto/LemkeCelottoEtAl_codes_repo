function [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = extract_reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot)
% spiking activity (seconds): blocks x channels x units
% lfp activity (lfp_samp_rate): blocks x channels
% beh markers (seconds)" blocks x 1

%% TRIAL INFORMATION

    count = 1;
    for block = 1:length(beh_markers)
        for trial = 1:length(beh_markers{block})
            if isnan(beh_markers{block}(trial))
                continue
            end
            trial_information(count).trial_video = trial_video{block}{trial};
            trial_information(count).timelock_frame = timelock_frame{block}(trial,:);
            trial_information(count).first_last_frame_time = first_last_frame_time{block}(trial,:);
            trial_information(count).lfp_same_rate = lfp_samp_rate;
            count = count + 1;
        end
    end

%% SPIKING ACTIVITY

    %%% M1

        unit_count = 1;
        for chan = 1:size(m1_spiking,2)
            for unit = 1:size(m1_spiking,3)
                if isempty(m1_spiking{1,chan,unit})
                    continue
                end

                raster = [];
                baseline = [];
                for block = 1:length(beh_markers)
                    for trial = 1:length(beh_markers{block})
                        if isnan(beh_markers{block}(trial))
                            continue
                        end
                        raster = [raster; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.001:5])];
                        if beh_markers{block}(trial)-15<0
                            baseline = [baseline; nan(1,10000)];
                        else
                            baseline = [baseline; histcounts(m1_spiking{block,chan,unit}(m1_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & m1_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.001:5])];
                        end
                    end
                end
                m1_reach_spiking(unit_count).m1_chan_unit = [chan unit];
                m1_reach_spiking(unit_count).reach_related_spiking = raster;
                m1_reach_spiking(unit_count).baseline_spiking = baseline;
                unit_count = unit_count + 1;

            end
        end
        
    %%% DLS

        unit_count = 1;
        for chan = 1:size(dls_spiking,2)
            for unit = 1:size(dls_spiking,3)
                if isempty(dls_spiking{1,chan,unit})
                    continue
                end

                raster = [];
                baseline = [];
                for block = 1:length(beh_markers)
                    for trial = 1:length(beh_markers{block})
                        if isnan(beh_markers{block}(trial))
                            continue
                        end
                        raster = [raster; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-5 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)+5)-beh_markers{block}(trial),[-5:0.001:5])];
                        if beh_markers{block}(trial)-15<0
                            baseline = [baseline; nan(1,10000)];
                        else
                            baseline = [baseline; histcounts(dls_spiking{block,chan,unit}(dls_spiking{block,chan,unit}>beh_markers{block}(trial)-15 & dls_spiking{block,chan,unit}<beh_markers{block}(trial)-5)-beh_markers{block}(trial)-10,[-5:0.001:5])];
                        end

                    end
                end
                dls_reach_spiking(unit_count).dls_chan_unit = [chan unit];
                dls_reach_spiking(unit_count).reach_related_spiking = raster;
                dls_reach_spiking(unit_count).baseline_spiking = baseline;
                unit_count = unit_count + 1;

            end
        end
    
%% LFP ACTIVITY 

    %%% GENERATE REFERENCED LFP 
    
        m1_lfp_ref = cell(size(m1_lfp));
        dls_lfp_ref = cell(size(dls_lfp));

        for block = 1:size(m1_lfp,1)

            %%% M1

                tmp_m1 = [];
                for m1_chan = 1:size(m1_lfp,2)
                    tmp_m1 = [tmp_m1; m1_lfp{block,m1_chan}];
                end
                tmp_m1 = median(tmp_m1);

                for m1_chan = 1:size(m1_lfp,2)
                    m1_lfp_ref{block,m1_chan} = m1_lfp{block,m1_chan}-tmp_m1;
                end      

            %%% DLS

                tmp_dls = [];
                for dls_chan = 1:size(dls_lfp,2)
                    tmp_dls = [tmp_dls; dls_lfp{block,dls_chan}];
                end
                tmp_dls = median(tmp_dls);

                for dls_chan = 1:size(dls_lfp,2)
                    dls_lfp_ref{block,dls_chan} = dls_lfp{block,dls_chan}-tmp_dls;
                end    

        end
    
    %%% M1
    
        chan_count = 1;
        for chan = 1:size(m1_lfp,2)

            raster = [];
            baseline = [];
            for block = 1:length(beh_markers)
                for trial = 1:length(beh_markers{block})
                    if isnan(beh_markers{block}(trial))
                        continue
                    end
                    raster = [raster; m1_lfp{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)+(5*lfp_samp_rate))];
                    if ((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate))<0
                        baseline = [baseline; nan(1,10173)];
                    else
                        baseline = [baseline; m1_lfp{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate))];
                    end
                end
            end

            m1_reach_lfp(chan_count).reach_related_lfp = raster;
            m1_reach_lfp(chan_count).baseline_lfp = baseline;
            chan_count = chan_count + 1;
            
        end    
          
        chan_count = 1;
        for chan = 1:size(m1_lfp_ref,2)

            raster = [];
            baseline = [];
            for block = 1:length(beh_markers)
                for trial = 1:length(beh_markers{block})
                    if isnan(beh_markers{block}(trial))
                        continue
                    end
                    raster = [raster; m1_lfp_ref{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)+(5*lfp_samp_rate))];
                    if ((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate))<0
                        baseline = [baseline; nan(1,10173)];
                    else
                        baseline = [baseline; m1_lfp_ref{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate))];
                    end
                end
            end

            m1_reach_lfp(chan_count).reach_related_lfp_ref = raster;
            m1_reach_lfp(chan_count).baseline_lfp_ref = baseline;
            chan_count = chan_count + 1;
            
        end

    %%% DLS
    
        chan_count = 1;
        for chan = 1:size(dls_lfp,2)

            raster = [];
            baseline = [];
            for block = 1:length(beh_markers)
                for trial = 1:length(beh_markers{block})
                    if isnan(beh_markers{block}(trial))
                        continue
                    end
                    raster = [raster; dls_lfp{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)+(5*lfp_samp_rate))];
                    if ((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate))<0
                        baseline = [baseline; nan(1,10173)];
                    else
                        baseline = [baseline; dls_lfp{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate))];
                    end
                end
            end

            dls_reach_lfp(chan_count).reach_related_lfp = raster;
            dls_reach_lfp(chan_count).baseline_lfp = baseline;
            chan_count = chan_count + 1;
            
        end
        
        chan_count = 1;
        for chan = 1:size(dls_lfp_ref,2)

            raster = [];
            baseline = [];
            for block = 1:length(beh_markers)
                for trial = 1:length(beh_markers{block})
                    if isnan(beh_markers{block}(trial))
                        continue
                    end
                    raster = [raster; dls_lfp_ref{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)+(5*lfp_samp_rate))];
                    if ((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate))<0
                        baseline = [baseline; nan(1,10173)];
                    else
                        baseline = [baseline; dls_lfp_ref{block,chan}((beh_markers{block}(trial)*lfp_samp_rate)-(15*lfp_samp_rate):(beh_markers{block}(trial)*lfp_samp_rate)-(5*lfp_samp_rate))];
                    end
                end
            end

            dls_reach_lfp(chan_count).reach_related_lfp_ref = raster;
            dls_reach_lfp(chan_count).baseline_lfp_ref = baseline;
            chan_count = chan_count + 1;
            
        end
               
%% PLOT

    if make_plot==1

        m1_mean_lfp =[];
        for chan = 1:length(m1_reach_lfp)
            m1_mean_lfp =[m1_mean_lfp; mean(m1_reach_lfp(chan).reach_related_lfp)];
        end

        dls_mean_lfp =[];
        for chan = 1:length(dls_reach_lfp)
            dls_mean_lfp =[dls_mean_lfp; mean(dls_reach_lfp(chan).reach_related_lfp)];
        end

        for unit = 1:length(m1_reach_spiking)

            figure;
                subplot(3,1,1); hold on;
                    imagesc(m1_reach_spiking(unit).reach_related_spiking)    
                    xlim([3000 7000])
                    ylim([1 size(m1_reach_spiking(unit).reach_related_spiking,1)])
                    colormap(flipud(gray))
                    title(['Day ' num2str(day) ' M1 channel ' num2str(m1_reach_spiking(unit).m1_chan_unit(1)) ' unit ' num2str(m1_reach_spiking(unit).m1_chan_unit(2))])
                subplot(3,1,2); hold on;
                    plot(smooth(mean(m1_reach_spiking(unit).reach_related_spiking),20),'color','k')
                    xlim([3000 7000])
                subplot(3,1,3); hold on;
                    plot(mean(m1_mean_lfp),'color',[0.5 0 0])
                    plot(mean(dls_mean_lfp),'color',[0 0 0.5])
                    xlim([3*lfp_samp_rate 7*lfp_samp_rate])
                    legend('M1 LFP','DLS LFP')

        end

        for unit = 1:length(dls_reach_spiking)

            figure;
                subplot(3,1,1); hold on;
                    imagesc(dls_reach_spiking(unit).reach_related_spiking)    
                    xlim([3000 7000])
                    ylim([1 size(dls_reach_spiking(unit).reach_related_spiking,1)])
                    colormap(flipud(gray))
                    title(['Day ' num2str(day) ' DLS channel ' num2str(dls_reach_spiking(unit).dls_chan_unit(1)) ' unit ' num2str(dls_reach_spiking(unit).dls_chan_unit(2))])
                subplot(3,1,2); hold on;
                    plot(smooth(mean(dls_reach_spiking(unit).reach_related_spiking),20),'color','k')
                    xlim([3000 7000])
                subplot(3,1,3); hold on;
                    plot(mean(m1_mean_lfp),'color',[0.5 0 0])
                    plot(mean(dls_mean_lfp),'color',[0 0 0.5])
                    xlim([3*lfp_samp_rate 7*lfp_samp_rate])
                    legend('M1 LFP','DLS LFP')

        end

    end

end
