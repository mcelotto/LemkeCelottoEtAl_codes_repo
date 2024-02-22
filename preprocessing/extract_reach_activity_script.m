%% EXTRACT DATA FROM RAW FILES
% generates a .mat file with reach-related neural data for each animal for each day of training

    clear all; clc; close all;

%% T102

day_blocks = {[2 3 7 8 9 11 12 13], ...
              [14 16 18 19 20 21], ...
              [22 24 25 26 27 29], ...
              [30 31 32 33 34 35 36 37 38], ...
              [39 43 44 45 47], ...
              [48 49 51 53 54 56], ...
              [57 58 59 60 61], ...
              [62 63 64 65 66]};
          
reach_blocks = {[3 4 6 7], [2 3 5],[2 3 6],[2 3 5 6 7 9],[2 3 5],[3 4 6],[2 3],[2 3 5]};

plx_files = {'T102_blocks_2_3_7_8_9_11_12_13-01', ...
             'T102_blocks_14_16_18_19_20_21-01', ...
             'T102_blocks_22_24_25_26_27_29-01', ...
             'T102_blocks_30_31_32_33_34_35_36_37_38-01', ...
             'T102_blocks_39_43_44_45_47-01', ...
             'T102_blocks_48_49_51_53_54_56-01', ...
             'T102_blocks_57_58_59_60_61-01', ...
             'T102_blocks_62_63_64_65_66-01'};
         
day_dirs = {['L:\videos\DLC\T102_kin3\day_1\'], ...
            ['L:\videos\DLC\T102_kin3\day_2\'], ...
            ['L:\videos\DLC\T102_kin3\day_3\'], ...
            ['L:\videos\DLC\T102_kin3\day_4\'], ...
            ['L:\videos\DLC\T102_kin3\day_5\'], ...
            ['L:\videos\DLC\T102_kin3\day_6\'], ...
            ['L:\videos\DLC\T102_kin3\day_7\'], ...
            ['L:\videos\DLC\T102_kin3\day_8\']};

time_each_block = {{'13h26m','13h42m','16h19m','16h32m'}, ...
                   {'13h44m','14h21m','16h38m'}, ...
                   {'12h50m','13h4m','16h4m'}, ...
                   {'12h2m','12h15m','15h21m','15h34m','15h45m','18h6m'}, ...
                   {'12h43m','13h4m','16h25m'}, ...
                   {'13h56m','14h16m','16h58m'}, ...
                   {'12h14m','12h29m','12h42m'}, ...
                   {'11h52m','12h6m','12h19m'}};
        
save_path = '\\mycloudpr4100\Public\SL_Backups\MEGA\data\T102\neural_data_reach';
make_plot = 0;

for day = 1:8
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
            for chan = 1:32
                tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
            end
            block_count = block_count + 1;
        end
        
        tmp_length = zeros(32,length(tmp_day_blocks));
        for chan = 1:32
            for block = 1:length(tmp_day_blocks)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(32,1) tmp_length];
        
        tmp_sort = cell(32,length(tmp_day_blocks));
        for chan = 1:32
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_day_blocks)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end

        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:16
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,count,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end
        
    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),16);
        dls_lfp = cell(length(reach_blocks{day}),16);
        
        for blocks = 1:length(reach_blocks{day})  
            tmp_LFP = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(tmp_day_blocks(reach_blocks{day}(blocks)))],'Type',4);
            
            count = 1;
            for m1_chan = 1:16
                m1_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(m1_chan,:);
                count = count + 1;
            end

            count = 1;
            for dls_chan = 17:32
                dls_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(dls_chan,:);
                count = count + 1;
            end

        end        
        lfp_samp_rate = tmp_LFP.streams.LFPs.fs;
 
	%%% LOAD BEHAVIORAL MARKERS
    
    beh_markers = cell(length(reach_blocks{day}),1);
    trial_video = cell(length(reach_blocks{day}),1);
    timelock_frame = cell(length(reach_blocks{day}),1);
    first_last_frame_time = cell(length(reach_blocks{day}),1);
    
    for blocks = 1:length(reach_blocks{day})
                
        %%% GET TDT TRIAL START TIMES
            clearvars final_pulses pellet_drops start_times
            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T102\neural_data\Block-' num2str(day_blocks{day}(reach_blocks{day}(blocks)))],'TYPE',4);
            fs = wave.streams.Wave.fs;
            wave = wave.streams.Wave.data(1,:);
            thr = max(wave)*0.4;
            temp=wave>thr;
            pulses = find(temp==1);
            pulse_diff = diff(pulses);
            start_pulses = find(pulse_diff>2);
            start_pulses = start_pulses+1;
            start_pulses = [1 start_pulses];
            for n=1:length(start_pulses)
                final_pulses(n) = pulses(start_pulses(n));
            end
            diff_final_pulses = diff(final_pulses);
            num_pulses = length(diff_final_pulses);
            diff_final_pulses = [301 diff_final_pulses 301 301];
            p_count= 1;
            s_count = 1;
            for n=1:num_pulses+1
                % pellet drops
                if (diff_final_pulses(n) > 200 && diff_final_pulses(n+1) >300)
                    pellet_drops(p_count) = final_pulses(n);
                    p_count = p_count + 1;
                end
                % start times
                if (diff_final_pulses(n) > 300 && diff_final_pulses(n+1) < 20 && diff_final_pulses(n+2)>200)
                    start_times(s_count) = final_pulses(n+1);
                    s_count = s_count + 1;
                end
            end
        
        %%% GET KINEMATICS
            dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\T102\reach_data\trajectories\day_' num2str(day) '\trajectories\*' time_each_block{day}{blocks} '*corrected*']);
            [~, i] = sort_nat({dlc_files(:).name});
            dlc_files = dlc_files(i);
            pellet_position = [];
            for file = 1:length(dlc_files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T102\reach_data\trajectories\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);
    
            mat_files = dir(['\\MyCloudPR4100\Public\SL_Backups\T102\reach_data\mat_videos\Day' num2str(day) '\*' time_each_block{day}{blocks} '*.mat']);
            [~, i] = sort_nat({mat_files(:).name});
            mat_files = mat_files(i);

            tmp_dlc_frame = [];
            tmp_wave_time = []; 
            
            tmp_beh_markers = [];            
            tmp_lock_frame = []; 
            tmp_trial_video = cell(1,length(mat_files));
            tmp_first_last_frame_time = [];
        
            for file = 1:length(mat_files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T102\reach_data\trajectories\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
                frame_times = load(['\\MyCloudPR4100\Public\SL_Backups\T102\reach_data\mat_videos\Day' num2str(day) '\' mat_files(file).name]);
                frame_times = frame_times.t_lat;
                frame_times_sec = frame_times;
                frame_times = frame_times*fs + start_times(file);
                %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;
                %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                    if (i<31 || i>size(trial_kin,1)-31)
                        tmp_beh_markers = [tmp_beh_markers nan];
                        tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                    else
                        tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                        tmp_lock_frame = [tmp_lock_frame; 1 i length(frame_times)];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; frame_times_sec(1)-frame_times_sec(i) frame_times_sec(end)-frame_times_sec(i)];
                    end
                tmp_dlc_frame = [tmp_dlc_frame size(trial_kin,1)];
                tmp_wave_time = [tmp_wave_time pellet_drops(file)-start_times(file)];
                tmp_trial_video{file} = dlc_files(file).name;
            end

            timelock_frame{blocks} = tmp_lock_frame;
            beh_markers{blocks} = tmp_beh_markers/fs;
            trial_video{blocks} = tmp_trial_video;
            first_last_frame_time{blocks} = tmp_first_last_frame_time;

%             figure 
%                 subplot(1,3,[1 2])
%                     hold on
%                     plot(wave);
%                     scatter(start_times,ones(1,length(start_times)))
%                     scatter(pellet_drops,ones(1,length(pellet_drops)))
%                     scatter(tmp_beh_markers,2*ones(1,length(tmp_beh_markers)))
%                     title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - vids: ' num2str(length(mat_files)) ' - dlc files: ' num2str(length(dlc_files)) ' - wave starts: ' num2str(length(start_times))])
%                 subplot(1,3,3)
%                     hold on
%                     plot(zscore(tmp_dlc_frame));
%                     plot(zscore(tmp_wave_time));
%                     title('DLC frames vs. wave length');
    end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = extract_reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['T102_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
    
end

%% T107

clc; clear all; close all;

day_blocks = {[1 2], ...
              [3 4 5 6 8], ...
              [9 10 11 12 13 14 15], ...
              [16 17 18 19 20 21 22 24 25], ...
              [26 27 30 31 36 39], ...
              [41 42 43 47 49 51 52]};
          
reach_blocks = {[], [2 5],[2 4 5 7],[4 5 6 7 9],[3 4 6],[2 5 6 7]};

plx_files = {'T107_blocks_1_2-01', ...
             'T107_blocks_3_4_5_6_8-01', ...
             'T107_blocks_9_10_11_12_13_14_15-01', ...
             'T107_blocks_16_17_18_19_20_21_22_24_25-01', ...
             'T107_blocks_26_27_30_31_36_39-01', ...
             'T107_blocks_41_42_43_47_49_51_52-01'};

day_dirs = {[], ...
    ['L:\videos\DLC\T107_kin\day_2\'], ...
    ['L:\videos\DLC\T107_kin\day_3\'], ...
    ['L:\videos\DLC\T107_kin\day_4\'], ...
    ['L:\videos\DLC\T107_kin\day_5\'], ...
    ['L:\videos\DLC\T107_kin\day_6\']};

time_each_block = {{}, ...
                   {'10h43m','15h21m'}, ...
                   {'11h3m','14h55m','15h15m','15h30m'}, ...
                   {'12h10m','12h25m','12h42m','12h59m','15h11m'}, ...
                   {'14h12m','14h27m','16h56m'}, ...
                   {'13h33m','16h11m','16h23m','16h33m'}};
               
save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\T107\neural_data_reach';
make_plot = 0;

for day = 2:6
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
            for chan = 1:32
                tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
            end
            block_count = block_count + 1;
        end
        
        tmp_length = zeros(32,length(tmp_day_blocks));
        for chan = 1:32
            for block = 1:length(tmp_day_blocks)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(32,1) tmp_length];
        
        tmp_sort = cell(32,length(tmp_day_blocks));
        for chan = 1:32
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_day_blocks)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end

        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:16
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,count,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end
        
    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),16);
        dls_lfp = cell(length(reach_blocks{day}),16);
        
        for blocks = 1:length(reach_blocks{day})  
            tmp_LFP = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(tmp_day_blocks(reach_blocks{day}(blocks)))],'Type',4);
            
            count = 1;
            for m1_chan = 1:16
                m1_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(m1_chan,:);
                count = count + 1;
            end

            count = 1;
            for dls_chan = 17:32
                dls_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(dls_chan,:);
                count = count + 1;
            end

        end        
        lfp_samp_rate = tmp_LFP.streams.LFPs.fs;
 
	%%% LOAD BEHAVIORAL MARKERS
    
    beh_markers = cell(length(reach_blocks{day}),1);
    trial_video = cell(length(reach_blocks{day}),1);
    timelock_frame = cell(length(reach_blocks{day}),1);
    first_last_frame_time = cell(length(reach_blocks{day}),1);

    for blocks = 1:length(reach_blocks{day})
                
        %%% GET TDT TRIAL START TIMES
        clearvars final_pulses pellet_drops start_times
        wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T107\neural_data\Block-' num2str(day_blocks{day}(reach_blocks{day}(blocks)))],'TYPE',4);
        fs = wave.streams.Wave.fs;
        wave = wave.streams.Wave.data(1,:);
        thr = max(wave)*0.4;
        temp=wave>thr;
        pulses = find(temp==1);
        pulse_diff = diff(pulses);
        start_pulses = find(pulse_diff>2);
        start_pulses = start_pulses+1;
        start_pulses = [1 start_pulses];
        for n=1:length(start_pulses)
            final_pulses(n) = pulses(start_pulses(n));
        end
        diff_final_pulses = diff(final_pulses);
        num_pulses = length(diff_final_pulses);
        diff_final_pulses = [301 diff_final_pulses 301 301];
        p_count= 1;
        s_count = 1;
        for n=1:num_pulses+1
            % pellet drops
            if (diff_final_pulses(n) > 200 && diff_final_pulses(n+1) >300)
                pellet_drops(p_count) = final_pulses(n);
                p_count = p_count + 1;
            end
            % start times
            if (diff_final_pulses(n) > 300 && diff_final_pulses(n+1) < 20 && diff_final_pulses(n+2)>200)
                start_times(s_count) = final_pulses(n+1);
                s_count = s_count + 1;
            end
        end
        
        %%% GET KINEMATICS
        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\T107\reach_data\trajectories\day_' num2str(day) '\trajectories\*' time_each_block{day}{blocks} '*corrected*']);
        [~, i] = sort_nat({dlc_files(:).name});
        dlc_files = dlc_files(i);
        pellet_position = [];
        for file = 1:length(dlc_files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T107\reach_data\trajectories\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
            tmp_pellet_x = trial_kin(1:25,5);
            tmp_pellet_y = trial_kin(1:25,6);
            tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
            tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
            pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
            pellet_position = [pellet_position; pellet_location];
        end
        pellet_position = nanmean(pellet_position);
    
        mat_files = dir(['\\MyCloudPR4100\Public\SL_Backups\T107\reach_data\mat_videos\Day' num2str(day) '\*' time_each_block{day}{blocks} '*.mat']);
        [~, i] = sort_nat({mat_files(:).name});
        mat_files = mat_files(i);
        
        tmp_dlc_frame = [];
        tmp_wave_time = []; 
        tmp_beh_markers = [];            
        tmp_lock_frame = []; 
        tmp_trial_video = cell(1,length(mat_files));
        tmp_first_last_frame_time = [];
        for file = 1:length(mat_files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T107\reach_data\trajectories\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
            frame_times = load(['\\MyCloudPR4100\Public\SL_Backups\T107\reach_data\mat_videos\Day' num2str(day) '\' mat_files(file).name]);
            frame_times = frame_times.t_lat;
            frame_times_sec = frame_times;
            frame_times = frame_times*fs + start_times(file);
            %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;
            %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);
            %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
            %%% skip bad trials
                if (i<31 || i>size(trial_kin,1)-31)
                    tmp_beh_markers = [tmp_beh_markers nan];
                    tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                    tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                else
                    tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                    tmp_lock_frame = [tmp_lock_frame; 1 i length(frame_times)];
                    tmp_first_last_frame_time = [tmp_first_last_frame_time; frame_times_sec(1)-frame_times_sec(i) frame_times_sec(end)-frame_times_sec(i)];
                end
            tmp_dlc_frame = [tmp_dlc_frame size(trial_kin,1)];
            tmp_wave_time = [tmp_wave_time pellet_drops(file)-start_times(file)];
            tmp_trial_video{file} = dlc_files(file).name;
        end
        timelock_frame{blocks} = tmp_lock_frame;
        beh_markers{blocks} = tmp_beh_markers/fs;
        trial_video{blocks} = tmp_trial_video;
        first_last_frame_time{blocks} = tmp_first_last_frame_time;
        
%         figure 
%             subplot(1,3,[1 2])
%                 hold on
%                 plot(wave);
%                 scatter(start_times,ones(1,length(start_times)))
%                 scatter(pellet_drops,ones(1,length(pellet_drops)))
%                 scatter(tmp_beh_markers,2*ones(1,length(tmp_beh_markers)))
%                 title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - vids: ' num2str(length(mat_files)) ' - dlc files: ' num2str(length(dlc_files)) ' - wave starts: ' num2str(length(start_times))])
%             subplot(1,3,3)
%                 hold on
%                 plot(zscore(tmp_dlc_frame));
%                 plot(zscore(tmp_wave_time));
%                 title('DLC frames vs. wave length');

    end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['T107_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
    
end

%% T110

clc; clear all; close all;

day_blocks = {[1],[2 5 7 9],[10 12 13 14],[15 16 17 18],[19 20 21 26],[27 28 31 32 33 34 35],[41 42 43 44 45]};

reach_blocks = {[],[2 3],[2 4],[2],[2 4],[3 4 5 7],[2 3 4]};

plx_files = {'T110_blocks_1-01','T110_blocks_2_5_7_9-01','T110_blocks_10_12_13_14-01','T110_blocks_15_16_17_18-01','T110_blocks_19_20_21_26-01','T110_blocks_27_28_31_32_33_34_35-01','T110_blocks_41_42_43_44_45-01'};

time_each_block = {{''}, ...
                   {'12h33m','13h36m'}, ...
                   {'13h54m','14h8m'}, ...
                   {'10h23m'}, ...
                   {'12h19m','15h14m'}, ...
                   {'11h42m','11h49m','12h0m','13h32m'}, ...
                   {'12h6m','14h55m','15h11m'}};
        
save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\T110\neural_data_reach';
make_plot = 0;

for day = 2:7
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T110\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T110\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
            for chan = 1:32
                tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
            end
            block_count = block_count + 1;
        end
        
        tmp_length = zeros(32,length(tmp_day_blocks));
        for chan = 1:32
            for block = 1:length(tmp_day_blocks)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(32,1) tmp_length];
        
        tmp_sort = cell(32,length(tmp_day_blocks));
        for chan = 1:32
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_day_blocks)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end

        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:16
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 17:32
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,count,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end
        
    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),16);
        dls_lfp = cell(length(reach_blocks{day}),16);
        
        for blocks = 1:length(reach_blocks{day})  
            tmp_LFP = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T110\neural_data\Block-' num2str(tmp_day_blocks(reach_blocks{day}(blocks)))],'Type',4);
            
            count = 1;
            for m1_chan = 1:16
                m1_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(m1_chan,:);
                count = count + 1;
            end

            count = 1;
            for dls_chan = 17:32
                dls_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(dls_chan,:);
                count = count + 1;
            end

        end        
        lfp_samp_rate = tmp_LFP.streams.LFPs.fs;
 
	%%% LOAD BEHAVIORAL MARKERS

        beh_markers = cell(length(reach_blocks{day}),1);
        trial_video = cell(length(reach_blocks{day}),1);
        timelock_frame = cell(length(reach_blocks{day}),1);
        first_last_frame_time = cell(length(reach_blocks{day}),1);

        for blocks = 1:length(reach_blocks{day})

            %%% GET TDT TRIAL START TIMES
            clearvars final_pulses pellet_drops start_times
            wave = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T110\neural_data\Block-' num2str(day_blocks{day}(reach_blocks{day}(blocks)))],'TYPE',4);
            fs = wave.streams.Wave.fs;
            wave = wave.streams.Wave.data(1,:);
            thr = max(wave)*0.4;
            temp=wave>thr;
            pulses = find(temp==1);
            pulse_diff = diff(pulses);
            start_pulses = find(pulse_diff>2);
            start_pulses = start_pulses+1;
            start_pulses = [1 start_pulses];
            for n=1:length(start_pulses)
                final_pulses(n) = pulses(start_pulses(n));
            end
            diff_final_pulses = diff(final_pulses);
            num_pulses = length(diff_final_pulses);
            diff_final_pulses = [301 diff_final_pulses 301 301];
            p_count= 1;
            s_count = 1;
            for n=1:num_pulses+1
                % pellet drops
                if (diff_final_pulses(n) > 200 && diff_final_pulses(n+1) >300)
                    pellet_drops(p_count) = final_pulses(n);
                    p_count = p_count + 1;
                end
                % start times
                if (diff_final_pulses(n) > 300 && diff_final_pulses(n+1) < 20 && diff_final_pulses(n+2)>200)
                    start_times(s_count) = final_pulses(n+1);
                    s_count = s_count + 1;
                end
            end

            %%% GET KINEMATICS
            dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\T110\reach_data\trajectories\day_' num2str(day) '\trajectories\*' time_each_block{day}{blocks} '*.csv']);
            [~, i] = sort_nat({dlc_files(:).name});
            dlc_files = dlc_files(i);
            pellet_position = [];
            for file = 1:length(dlc_files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T110\reach_data\trajectories\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            mat_files = dir(['\\MyCloudPR4100\Public\SL_Backups\T110\reach_data\mat_videos\day_' num2str(day) '\*' time_each_block{day}{blocks} '*.mat']);
            [~, i] = sort_nat({mat_files(:).name});
            mat_files = mat_files(i);

            if day==2 & blocks==1
                mat_files = mat_files(1:15);
            end

            if day==2 & blocks==2
                dlc_files = dlc_files(2:end);
                mat_files = mat_files(2:end);
            end     

            if day==3 & blocks==1
                pellet_drops = pellet_drops(2:end);
            end

            if day==4 & blocks==1
                dlc_files = dlc_files(1:28);
                mat_files = mat_files(1:28);
            end        

            tmp_dlc_frame = [];
            tmp_wave_time = []; 
            tmp_beh_markers = [];            
            tmp_lock_frame = []; 
            tmp_trial_video = cell(1,length(mat_files));
            tmp_first_last_frame_time = [];
            for file = 1:length(mat_files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T110\reach_data\trajectories\day_' num2str(day) '\trajectories\' dlc_files(file).name]);
                frame_times = load(['\\MyCloudPR4100\Public\SL_Backups\T110\reach_data\mat_videos\day_' num2str(day) '\' mat_files(file).name]);
                frame_times = frame_times.t_lat;
                frame_times_sec = frame_times;
                frame_times = frame_times*fs + start_times(file);
                %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;
                %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                    if (i<31 || i>size(trial_kin,1)-31)
                        tmp_beh_markers = [tmp_beh_markers nan];
                        tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                    else
                        tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                        tmp_lock_frame = [tmp_lock_frame; 1 i length(frame_times)];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; frame_times_sec(1)-frame_times_sec(i) frame_times_sec(end)-frame_times_sec(i)];
                    end
                tmp_dlc_frame = [tmp_dlc_frame size(trial_kin,1)];
                tmp_wave_time = [tmp_wave_time pellet_drops(file)-start_times(file)];
                tmp_trial_video{file} = dlc_files(file).name;
            end
            timelock_frame{blocks} = tmp_lock_frame;
            beh_markers{blocks} = tmp_beh_markers/fs;
            trial_video{blocks} = tmp_trial_video;
            first_last_frame_time{blocks} = tmp_first_last_frame_time;

%             figure 
%                 subplot(1,3,[1 2])
%                     hold on
%                     plot(wave);
%                     scatter(start_times,ones(1,length(start_times)))
%                     scatter(pellet_drops,ones(1,length(pellet_drops)))
%                     scatter(tmp_beh_markers,2*ones(1,length(tmp_beh_markers)))
%                     title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - vids: ' num2str(length(mat_files)) ' - dlc files: ' num2str(length(dlc_files)) ' - wave starts: ' num2str(length(start_times))])
%                 subplot(1,3,3)
%                     hold on
%                     plot(zscore(tmp_dlc_frame));
%                     plot(zscore(tmp_wave_time));
%                     title('DLC frames vs. wave length');
        end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['T110_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
    
end

%% T200

clc; clear all; close all;

day_blocks = {[1], ...
              [2 7 8], ...
              [9 11 12 13 14], ...
              [15 16 17 18 19], ...
              [20 21 22 23 24], ...
              [25 26 27], ...
              [28 29 30], ...
              [32 33 34 35 36], ...
              [37 38 39 40 41], ...
              [42 43 44 45 46 47]};
          
reach_blocks = {[],[2],[3 4],[3 4],[2 3],[2],[2],[2 3],[2 3],[2]};

plx_files = {'T200_blocks_1-01', ...
             'T200_blocks_2_7_8-01', ...
             'T200_blocks_9_11_12_13_14-01', ...
             'T200_blocks_15_16_17_18_19-01', ...
             'T200_blocks_20_21_22_23_24-01', ...
             'T200_blocks_25_26_27-01', ...
             'T200_blocks_28_29_30-01', ...
             'T200_blocks_32_33_34_35_36-01', ...
             'T200_blocks_37_38_39_40_41-01', ...
             'T200_blocks_42_43_44_45_46_47-01'};

save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\T200\neural_data_reach';
make_plot = 0;

for day = 2:10
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(32,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
            for chan = 1:32
                tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
            end
            block_count = block_count + 1;
        end
        
        tmp_length = zeros(32,length(tmp_day_blocks));
        for chan = 1:32
            for block = 1:length(tmp_day_blocks)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(32,1) tmp_length];
        
        tmp_sort = cell(32,length(tmp_day_blocks));
        for chan = 1:32
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_day_blocks)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end
        
        m1_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 17:32
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 1:16
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end
 
    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),16);
        dls_lfp = cell(length(reach_blocks{day}),16);
        
        for blocks = 1:length(reach_blocks{day})  
            tmp_LFP = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_day_blocks(reach_blocks{day}(blocks)))],'Type',4);
            
            count = 1;
            for m1_chan = 17:32
                m1_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(m1_chan,:);
                count = count + 1;
            end

            count = 1;
            for dls_chan = 1:16
                dls_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(dls_chan,:);
                count = count + 1;
            end

        end        
        lfp_samp_rate = tmp_LFP.streams.LFPs.fs;
        
	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);
        trial_video = cell(length(reach_blocks{day}),1);
        timelock_frame = cell(length(reach_blocks{day}),1);
        first_last_frame_time = cell(length(reach_blocks{day}),1);

        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\Lemke_etal_sleep\data\T200\behavioral_data\day_' num2str(day) '\trajectories\*corrected*']);
        [~, i] = sort_nat({dlc_files.name});
        dlc_files = dlc_files(i);     

        for blocks = 1:length(reach_blocks{day})

            tmp_block = tmp_day_blocks(reach_blocks{day}(blocks));

            %%% LOAD FRAME TIMES
            Vid = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T200\neural_data\Block-' num2str(tmp_block)],'STORE','Vid0');
            [~,i] = find(Vid.scalars.Vid0.data>1e8);
            for n_i = 1:length(i)
                Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
            end  
            Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
            first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
            trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
            trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];      

            %%% GET KINEMATICS
            kin_all = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\Lemke_etal_sleep\data\T200\behavioral_data\day_' num2str(day) '\trajectories\' dlc_files(blocks).name]);

            pellet_position = [];
            for n_trial = 1:length(trial_starts_idx)
                %%% get timesstamps (ts) for frames in trial
                trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                %%% interpolate to get missing ts
                interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));
                %%% get trial kinematics
                trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);
                %%% get pellet location from first 25 frames
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            tmp_beh_markers = [];
            tmp_dlc_frame = [];
            tmp_wave_time = []; 
            tmp_lock_frame = []; 
            tmp_trial_video = cell(1,length(trial_starts_idx));
            tmp_first_last_frame_time = [];
            for n_trial = 1:length(trial_starts_idx)
                %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));
                %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);
                %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;
                %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                    if (i<31 || i>size(trial_kin,1)-31)
                        tmp_beh_markers = [tmp_beh_markers nan];
                        tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                    else
                        tmp_beh_markers = [tmp_beh_markers interp_trial_ts(i)];
                        tmp_lock_frame = [tmp_lock_frame; trial_kin(1,1) trial_kin(i,1) trial_kin(end,1)];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; interp_trial_ts(1)-interp_trial_ts(i) interp_trial_ts(end)-interp_trial_ts(i)];
                    end
                tmp_trial_video{n_trial} = dlc_files(blocks).name;
            end
            timelock_frame{blocks} = tmp_lock_frame;
            beh_markers{blocks} = tmp_beh_markers;
            trial_video{blocks} = tmp_trial_video;
            first_last_frame_time{blocks} = tmp_first_last_frame_time;            
        end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['T200_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)

end

%% T201

clc; clear all; close all;

day_blocks = {[1], ...
              [2 3], ...
              [4 5 6 7], ...
              [8 9 10], ...
              [11 12 13], ...
              [14 15 16 17], ...
              [18 19 20], ...
              [21 22 23 24 25 26], ...
              [27 28 29], ...
              [30 31 32], ...
              [33 34 35 36 37], ...
              [38 39 40], ...
              [41 42 43 44], ...
              [45 46 47 48], ...
              [49 50 51 52], ...
              [53 54 55 56], ...
              [57 58 59 60 61]};

reach_blocks = {[], ...
                [], ...
                [3], ...
                [2], ...
                [2], ...
                [3], ...
                [2], ...
                [2 3 4], ...
                [2], ...
                [2], ...
                [3], ...
                [2], ...
                [2], ...
                [2], ...
                [2], ...
                [2], ...
                [3]};

plx_files = {'T201_blocks_1-01', ...
             'T201_blocks_2_3-01', ...
             'T201_blocks_4_5_6_7-01', ...
             'T201_blocks_8_9_10-01', ...
             'T201_blocks_11_12_13-01', ...
             'T201_blocks_14_15_16_17-01', ...
             'T201_blocks_18_19_20-01', ...
             'T201_blocks_21_22_23_24_25_26-01', ...
             'T201_blocks_27_28_29-01', ... 
             'T201_blocks_30_31_32-01', ...
             'T201_blocks_33_34_35_36_37-01', ...
             'T201_blocks_38_39_40-01', ...
             'T201_blocks_41_42_43_44-01', ...
             'T201_blocks_45_46_47_48-01', ...
             'T201_blocks_49_50_51_52-01', ...
             'T201_blocks_53_54_55_56-01', ...
             'T201_blocks_57_58_59_60_61-01'};

save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\T201\neural_data_reach';
make_plot = 0;

for day = 4:17
    
    tmp_day_blocks = day_blocks{day};
    
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_ts = cell(64,length(tmp_day_blocks));
        block_count = 1;
        for block = tmp_day_blocks
            if block == 41
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(block)],'Type',3,'T2',7200);
            else
                tmp_num = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(block)],'Type',3,'NODATA',1);
            end
            for chan = 1:64
                tmp_ts{chan,block_count} = tmp_num.snips.eNeu.ts(tmp_num.snips.eNeu.chan==chan);
            end
            block_count = block_count + 1;
        end
        
        tmp_length = zeros(64,length(tmp_day_blocks));
        for chan = 1:64
            for block = 1:length(tmp_day_blocks)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(64,1) tmp_length];
        
        tmp_sort = cell(64,length(tmp_day_blocks));
        for chan = 1:64
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_day_blocks)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end
        
        m1_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = 1:2:63
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),16,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = 34:2:64
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end

    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),32);
        dls_lfp = cell(length(reach_blocks{day}),16);
        
        for blocks = 1:length(reach_blocks{day})  
            tmp_LFP = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_day_blocks(reach_blocks{day}(blocks)))],'Type',4);
            
            count = 1;
            for m1_chan = 1:2:64
                m1_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(m1_chan,:);
                count = count + 1;
            end

            count = 1;
            for dls_chan = 34:2:64
                dls_lfp{blocks,count} = tmp_LFP.streams.LFPs.data(dls_chan,:);
                count = count + 1;
            end

        end        
        lfp_samp_rate = tmp_LFP.streams.LFPs.fs;
        
	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);
        trial_video = cell(length(reach_blocks{day}),1);
        timelock_frame = cell(length(reach_blocks{day}),1);
        first_last_frame_time = cell(length(reach_blocks{day}),1);
        
        dlc_files = dir(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\Lemke_etal_sleep\data\T201\behavioral_data\day_' num2str(day) '\trajectories\*corrected*']);
        [~, i] = sort_nat({dlc_files.name});
        dlc_files = dlc_files(i);     

        for blocks = 1:length(reach_blocks{day})

            tmp_block = tmp_day_blocks(reach_blocks{day}(blocks));

            %%% LOAD FRAME TIMES
            Vid = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T201\neural_data\Block-' num2str(tmp_block)],'STORE','Vid0');
            [~,i] = find(Vid.scalars.Vid0.data>1e8);
            for n_i = 1:length(i)
                Vid.scalars.Vid0.data(i(n_i)) = round(mean(Vid.scalars.Vid0.data(i(n_i)-1):Vid.scalars.Vid0.data(i(n_i)+1)));
            end  
            Vid.scalars.Vid0.data = Vid.scalars.Vid0.data+1;
            first_frames = find(diff(Vid.scalars.Vid0.ts)>8);
            trial_starts_idx = [Vid.scalars.Vid0.data(1) Vid.scalars.Vid0.data(first_frames+1)];
            trial_ends_idx = [Vid.scalars.Vid0.data(first_frames) Vid.scalars.Vid0.data(end)];      

            %%% GET KINEMATICS
            kin_all = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\PROJECTS\Lemke_etal_sleep\data\T201\behavioral_data\day_' num2str(day) '\trajectories\' dlc_files(blocks).name]);

            if day==11 || day ==12 || day ==13 || day ==15 || day ==16 || day ==17
                trial_starts_idx = trial_starts_idx(1:end-1);
                trial_ends_idx = trial_ends_idx(1:end-1);
            end
                
            pellet_position = [];
            for n_trial = 1:length(trial_starts_idx)
                %%% get timesstamps (ts) for frames in trial
                trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                %%% interpolate to get missing ts
                interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));
                %%% get trial kinematics
                trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);
                %%% get pellet location from first 25 frames
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);

            tmp_beh_markers = [];
            tmp_dlc_frame = [];
            tmp_wave_time = []; 
            tmp_lock_frame = []; 
            tmp_trial_video = cell(1,length(trial_starts_idx));
            tmp_first_last_frame_time = [];
            for n_trial = 1:length(trial_starts_idx)
                %%% get timesstamps (ts) for frames in trial
                    trial_frames_idx = Vid.scalars.Vid0.data(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                    trial_frames_ts = Vid.scalars.Vid0.ts(Vid.scalars.Vid0.data>=trial_starts_idx(n_trial) & Vid.scalars.Vid0.data<=trial_ends_idx(n_trial));
                %%% interpolate to get missing ts
                    interp_trial_ts = interp1(trial_frames_idx,trial_frames_ts,trial_starts_idx(n_trial):trial_ends_idx(n_trial));
                %%% get trial kinematics
                    trial_kin = kin_all(trial_starts_idx(n_trial):trial_ends_idx(n_trial),:);
                %%% get paw trajectory
                    tmp_x = trial_kin(:,2);
                    tmp_y = trial_kin(:,3);
                    tmp_x(trial_kin(:,4)<.99) = NaN;
                    tmp_y(trial_kin(:,4)<.99) = NaN;
                %%% get paw velocity
                    diff_x = diff(tmp_x);
                    diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                    [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                    if (i<31 || i>size(trial_kin,1)-31)
                        tmp_beh_markers = [tmp_beh_markers nan];
                        tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                    else
                        tmp_beh_markers = [tmp_beh_markers interp_trial_ts(i)];
                        tmp_lock_frame = [tmp_lock_frame; trial_kin(1,1) trial_kin(i,1) trial_kin(end,1)];
                        tmp_first_last_frame_time = [tmp_first_last_frame_time; interp_trial_ts(1)-interp_trial_ts(i) interp_trial_ts(end)-interp_trial_ts(i)];
                    end
                tmp_trial_video{n_trial} = dlc_files(blocks).name;
            end
            timelock_frame{blocks} = tmp_lock_frame;
            beh_markers{blocks} = tmp_beh_markers;
            trial_video{blocks} = tmp_trial_video;
            first_last_frame_time{blocks} = tmp_first_last_frame_time;            
        end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['T201_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
   
end

%% SLDD1

clc; clear all; close all;

day_blocks = {{'SLDD1-191113-112021'}, ...
                {'SLDD1-191114-095200'}, ...
                {'SLDD1-191115-100221','SLDD1-191115-113606','SLDD1-191115-123455'}, ...
                {'SLDD1-191116-094344','SLDD1-191116-114540','SLDD1-191116-123308'}, ...
                {'SLDD1-191117-095617','SLDD1-191117-115905','SLDD1-191117-123758'}, ...
                {'SLDD1-191118-100800','SLDD1-191118-120904','SLDD1-191118-124021','SLDD1-191118-132255'}, ...
                {'SLDD1-191119-094610','SLDD1-191119-114803','SLDD1-191119-124914','SLDD1-191119-133446'}, ...
                {'SLDD1-191120-085736','SLDD1-191120-110014','SLDD1-191120-114742'}, ...
                {'SLDD1-191121-105732','SLDD1-191121-132514','SLDD1-191121-141748','SLDD1-191121-151605'}, ...
                {'SLDD1-191122-101205','SLDD1-191122-121505','SLDD1-191122-130930'}, ...
                {'SLDD1-191123-102258','SLDD1-191123-122504','SLDD1-191123-125725','SLDD1-191123-135233'}, ...
                {'SLDD1-191124-084031','SLDD1-191124-100953','SLDD1-191124-114327','SLDD1-191124-124703','SLDD1-191124-144901'}, ...
                {'SLDD1-191125-085338','SLDD1-191125-105650','SLDD1-191125-115721'}, ...
                {'SLDD1-191126-124718','SLDD1-191126-142017','SLDD1-191126-155329'}};

reach_blocks = {[],[],[2],[2],[2],[3],[3],[2],[3],[2],[3],[3],[2],[2]};      

plx_files = {'', ...
            '', ...
            'sldd1-191115-100221_mrg-01', ...
            'sldd1-191116-094344_mrg-01', ...
            'sldd1-191117-095617_mrg-01', ...
            'sldd1-191118-100800_mrg-01', ...
            'sldd1-191119-094610_mrg-01', ...
            'sldd1-191120-085736_mrg-01', ...
            'sldd1-191121-105732_mrg-01', ...
            'sldd1-191122-101205_mrg-01', ...
            'sldd1-191123-102258_mrg-01', ...
            'sldd1-191124-084031_mrg-01', ...
            'sldd1-191125-085338_mrg-01', ...
            'sldd1-191126-124718_mrg-01'};

save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\SLDD1\neural_data_reach';
make_plot = 0;

for day = 3:14
        
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\PLX_sorts\day_' num2str(day) '\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_blocks = ls(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\PLX_sorts\day_' num2str(day) '\*.txt']);
        tmp_blocks_cell = cell(size(tmp_blocks,1),1);
        for block = 1:size(tmp_blocks,1)
            tmp_blocks_cell{block} = deblank(tmp_blocks(block,:));
        end
        tmp_ind = ~contains(tmp_blocks_cell,'mrg');
        tmp_blocks_cell = tmp_blocks_cell(tmp_ind);

        tmp_ts = cell(64,length(tmp_blocks_cell));
        for block = 1:length(tmp_blocks_cell)
            disp(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\PLX_sorts\day_' num2str(day) '\' tmp_blocks_cell{block}]);
            tmp_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\PLX_sorts\day_' num2str(day) '\' tmp_blocks_cell{block}]);
            tmp_spiking = table2array(tmp_spiking);
            for chan = 1:64
                tmp_ts{chan,block} = tmp_spiking(tmp_spiking(:,1)==chan,2);
            end
        end
        
        tmp_length = zeros(64,length(tmp_blocks_cell));
        for chan = 1:64
            for block = 1:length(tmp_blocks_cell)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(64,1) tmp_length];
        
        tmp_sort = cell(64,length(tmp_blocks_cell));
        for chan = 1:64
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_blocks_cell)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end
        
        m1_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end

    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),32);
        dls_lfp = cell(length(reach_blocks{day}),32);
        for blocks = 1:length(reach_blocks{day})
            
            count = 1;
            for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
                tmp_lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\RS4_data\' day_blocks{day}{reach_blocks{day}(blocks)}],'CHANNEL',m1_chan);
                m1_lfp{blocks,count} = downsample(tmp_lfp.RSn1.data,20);
                count = count + 1;
            end
            
            count = 1;
            for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
                tmp_lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\RS4_data\' day_blocks{day}{reach_blocks{day}(blocks)}],'CHANNEL',dls_chan);
                dls_lfp{blocks,count} = downsample(tmp_lfp.RSn1.data,20);
                count = count + 1;
            end
            
        end
        lfp_samp_rate = tmp_lfp.RSn1.fs/20;        
        
	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);
        trial_video = cell(length(reach_blocks{day}),1);
        timelock_frame = cell(length(reach_blocks{day}),1);
        first_last_frame_time = cell(length(reach_blocks{day}),1);

        for blocks = 1:length(reach_blocks{day})

            tmp_block = day_blocks{day}{reach_blocks{day}(blocks)};

            frames = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\neural_data\TDT_tanks\' tmp_block]);
            frames = frames.epocs.PtC1.onset;
            trial_starts = [1; find(diff(frames)>1)+1];
            trial_ends = [find(diff(frames)>1); length(frames)];
            
            files = dir(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\reach_data\day_' num2str(day-2) '\trajectories\*SLDD1*.csv']);
            if day == 14
                files = dir(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\reach_data\day_' num2str(day-2) '\trajectories\*corrected.csv']);
            end
            [~, i] = sort_nat({files(:).name});
            files = files(i);
            
            %%% FIND PELLET POSITION
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\reach_data\day_' num2str(day-2) '\trajectories\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);
   
            tmp_tdt_length = [];
            tmp_dlc_length = [];
            tmp_beh_markers = [];
            tmp_lock_frame = []; 
            tmp_trial_video = cell(1,length(files));
            tmp_first_last_frame_time = [];
            for file = 1:length(files)
                frame_times = frames(trial_starts(file):trial_ends(file));
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\SLDD1\reach_data\day_' num2str(day-2) '\trajectories\' files(file).name]);
                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;
                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                if (i<31 || i>size(trial_kin,1)-31)
                    tmp_beh_markers = [tmp_beh_markers nan];
                    tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                    tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                else
                    tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                    tmp_lock_frame = [tmp_lock_frame; 1 i length(frame_times)];
                    tmp_first_last_frame_time = [tmp_first_last_frame_time; frame_times(1)-frame_times(i) frame_times(end)-frame_times(i)];                   
                end
                tmp_dlc_length = [tmp_dlc_length size(trial_kin,1)];
                tmp_tdt_length = [tmp_tdt_length frames(trial_ends(file))-frames(trial_starts(file))];
                tmp_trial_video{file} = files(file).name;
            end
            timelock_frame{blocks} = tmp_lock_frame;
            beh_markers{blocks} = tmp_beh_markers;
            trial_video{blocks} = tmp_trial_video;
            first_last_frame_time{blocks} = tmp_first_last_frame_time;

            figure
                hold on
                plot(zscore(tmp_dlc_length));
                plot(zscore(tmp_tdt_length));
                title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - dlc files: ' num2str(length(files)) ' - wave starts: ' num2str(length(trial_starts))])
        end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['SLDD1_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
    
 
end

%% SLDD2

clc; clear all; close all;
            
day_blocks = {{'SLDD2-191113-143329'}, ...
                {'SLDD2-191114-121834'}, ...
                {'SLDD2-191115-154323'}, ...
                {'SLDD2-191116-145119','SLDD2-191116-170045','SLDD2-191116-193647'}, ...
                {'SLDD2-191117-150115','SLDD2-191117-170253','SLDD2-191117-181931'}, ...
                {'SLDD2-191118-154100','SLDD2-191118-180452','SLDD2-191118-192527','SLDD2-191118-205741'}, ...
                {'SLDD2-191119-160043','SLDD2-191119-180437','SLDD2-191119-183709','SLDD2-191119-190631','SLDD2-191119-193444'}, ...
                {'SLDD2-191120-140404','SLDD2-191120-160657','SLDD2-191120-170507'}, ...
                {'SLDD2-191121-170233','SLDD2-191121-180559','SLDD2-191121-182818','SLDD2-191121-192200'}, ...
                {'SLDD2-191122-152109','SLDD2-191122-172516','SLDD2-191122-184340','SLDD2-191122-204428'}, ...
                {'SLDD2-191123-173250','SLDD2-191123-193545','SLDD2-191123-205347'}, ...
                {'SLDD2-191124-152841','SLDD2-191124-173146','SLDD2-191124-183235'}, ...
                {'SLDD2-191125-141910','SLDD2-191125-162007','SLDD2-191125-171539'}, ...
                {'SLDD2-191126-173517','SLDD2-191126-191048','SLDD2-191126-195858'}};
            
reach_blocks = {[], ...
                [], ...
                [], ...
                [2], ...
                [2], ...
                [2], ...
                [2 3 4], ...
                [2], ...
                [3], ...
                [2], ...
                [2], ...
                [2], ...
                [2], ...
                [2]};

plx_files = {'sldd2-191113-143329-01', ...
            'sldd2-191114-121834-01', ...
            'sldd2-191115-154323-01', ...
            'sldd2-191116-145119_mrg-01', ...
            'sldd2-191117-150115_mrg-01', ...
            'sldd2-191118-154100_mrg-01', ...
            'sldd2-191119-160043_mrg-01', ...
            'sldd2-191120-140404_mrg-01', ...
            'sldd2-191121-170233_mrg-01', ...
            'sldd2-191122-152109_mrg-01', ...
            'sldd2-191123-173250_mrg-01', ...
            'sldd2-191124-152841_mrg-01', ...
            'sldd2-191125-141910_mrg-01', ...
            'sldd2-191126-173517_mrg-01'};

save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\SLDD2\neural_data_reach';
make_plot = 0;

for day = 4:13
        
    %%% LOAD SPIKING DATA
   
        plx_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\' plx_files{day} '.txt']);
        plx_spiking = table2array(plx_spiking);
        
        tmp_blocks = ls(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\*.txt']);
        tmp_blocks_cell = cell(size(tmp_blocks,1),1);
        for block = 1:size(tmp_blocks,1)
            tmp_blocks_cell{block} = deblank(tmp_blocks(block,:));
        end
        tmp_ind = ~contains(tmp_blocks_cell,'mrg');
        tmp_blocks_cell = tmp_blocks_cell(tmp_ind);

        tmp_ts = cell(64,length(tmp_blocks_cell));
        for block = 1:length(tmp_blocks_cell)
            disp(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\' tmp_blocks_cell{block}]);
            tmp_spiking = readtable(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\PLX_sorts\day_' num2str(day) '\' tmp_blocks_cell{block}]);
            tmp_spiking = table2array(tmp_spiking);
            for chan = 1:64
                tmp_ts{chan,block} = tmp_spiking(tmp_spiking(:,1)==chan,2);
            end
        end
        
        tmp_length = zeros(64,length(tmp_blocks_cell));
        for chan = 1:64
            for block = 1:length(tmp_blocks_cell)
                tmp_length(chan,block) = length(tmp_ts{chan,block});
            end
        end  
        tmp_length = [zeros(64,1) tmp_length];
        
        tmp_sort = cell(64,length(tmp_blocks_cell));
        for chan = 1:64
            chan_sorts = plx_spiking(plx_spiking(:,1)==chan,:);
            for block = 1:length(tmp_blocks_cell)
                tmp_sort{chan,block} = chan_sorts(1+sum(tmp_length(chan,1:block)):sum(tmp_length(chan,1:block+1)),2);
            end
        end
        
        m1_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
                for m1_unit = 1:max(plx_spiking(:,2))
                    m1_spiking{blocks,count,m1_unit} = tmp_ts{m1_chan,reach_blocks{day}(blocks)}(tmp_sort{m1_chan,reach_blocks{day}(blocks)}==m1_unit);
                end
                count = count + 1;
            end
        end
        
        dls_spiking = cell(length(reach_blocks{day}),32,max(plx_spiking(:,2)));
        for blocks = 1:length(reach_blocks{day})
            count = 1;
            for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
                for dls_unit = 1:max(plx_spiking(:,2))
                    dls_spiking{blocks,dls_chan,dls_unit} = tmp_ts{dls_chan,reach_blocks{day}(blocks)}(tmp_sort{dls_chan,reach_blocks{day}(blocks)}==dls_unit);
                end
                count = count + 1;
            end
        end

    %%% LOAD LFP DATA
        
        m1_lfp = cell(length(reach_blocks{day}),32);
        dls_lfp = cell(length(reach_blocks{day}),32);
        for blocks = 1:length(reach_blocks{day})
            
            count = 1;
            for m1_chan = [8 4 12 5 22 6 21 17 13 14 28 1 30 25 29 20 52 53 54 61 49 41 38 36 62 57 60 44 37 46 45 33]
                tmp_lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' day_blocks{day}{reach_blocks{day}(blocks)}],'CHANNEL',m1_chan);
                m1_lfp{blocks,count} = downsample(tmp_lfp.RSn1.data,20);
                count = count + 1;
            end
            
            count = 1;
            for dls_chan = [2 3 9 7 11 10 16 15 19 18 24 23 27 26 32 31 35 34 40 39 43 42 48 47 51 50 56 55 59 58 64 63]
                tmp_lfp = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\RS4_data\' day_blocks{day}{reach_blocks{day}(blocks)}],'CHANNEL',dls_chan);
                dls_lfp{blocks,count} = downsample(tmp_lfp.RSn1.data,20);
                count = count + 1;
            end
            
        end
        lfp_samp_rate = tmp_lfp.RSn1.fs/20;        
        
	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(length(reach_blocks{day}),1);
        trial_video = cell(length(reach_blocks{day}),1);
        timelock_frame = cell(length(reach_blocks{day}),1);
        first_last_frame_time = cell(length(reach_blocks{day}),1);

        for blocks = 1:length(reach_blocks{day})

            tmp_block = day_blocks{day}{reach_blocks{day}(blocks)};

            frames = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\neural_data\TDT_tanks\' tmp_block]);
            frames = frames.epocs.PtC1.onset;
            trial_starts = [1; find(diff(frames)>1)+1];
            trial_ends = [find(diff(frames)>1); length(frames)];

            if day == 7
                if blocks == 1
                    files = dir(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\reach_data\day_' num2str(day-3) '\trajectories\*18h4m*corrected*']);
                    [~, i] = sort_nat({files(:).name});
                    files = files(i);
                    trial_starts = trial_starts([1:20 22:end]);
                    trial_ends = trial_ends([1:20 22:end]);
                elseif blocks == 2
                    files = dir(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\reach_data\day_' num2str(day-3) '\trajectories\*18h37m*corrected*']);
                    [~, i] = sort_nat({files(:).name});
                    files = files(i);
                    files = files(1:end-1);
                else
                    files = dir(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\reach_data\day_' num2str(day-3) '\trajectories\*19h6m*corrected*']);
                    [~, i] = sort_nat({files(:).name});
                    files = files(i);
                end
            else
                files = dir(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\reach_data\day_' num2str(day-3) '\trajectories\*corrected*']);
                [~, i] = sort_nat({files(:).name});
                files = files(i);
            end
            
            %%% FIND PELLET POSITION
            pellet_position = [];
            for file = 1:length(files)
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\reach_data\day_' num2str(day-3) '\trajectories\' files(file).name]);
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                pellet_position = [pellet_position; pellet_location];
            end
            pellet_position = nanmean(pellet_position);
   
            tmp_tdt_length = [];
            tmp_dlc_length = [];
            tmp_beh_markers = [];
            tmp_lock_frame = []; 
            tmp_trial_video = cell(1,length(files));
            tmp_first_last_frame_time = []; 
            for file = 1:length(files)
                frame_times = frames(trial_starts(file):trial_ends(file));
                trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\SLDD2\reach_data\day_' num2str(day-3) '\trajectories\' files(file).name]);
                %%% get paw trajectory
                tmp_x = trial_kin(:,2);
                tmp_y = trial_kin(:,3);
                tmp_x(trial_kin(:,4)<.99) = NaN;
                tmp_y(trial_kin(:,4)<.99) = NaN;
                %%% get paw velocity
                diff_x = diff(tmp_x);
                diff_y = diff(tmp_y);
                %%% find frame with paw closest to pellet
                [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
                %%% skip bad trials
                if (i<31 || i>size(trial_kin,1)-31)
                    tmp_beh_markers = [tmp_beh_markers nan];
                    tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                    tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
                else
                    tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                    tmp_lock_frame = [tmp_lock_frame; 1 i length(frame_times)];
                    tmp_first_last_frame_time = [tmp_first_last_frame_time; frame_times(1)-frame_times(i) frame_times(end)-frame_times(i)];                   
                end
                tmp_dlc_length = [tmp_dlc_length size(trial_kin,1)];
                tmp_tdt_length = [tmp_tdt_length frames(trial_ends(file))-frames(trial_starts(file))];
                tmp_trial_video{file} = files(file).name;
            end
            timelock_frame{blocks} = tmp_lock_frame;
            beh_markers{blocks} = tmp_beh_markers;
            trial_video{blocks} = tmp_trial_video;
            first_last_frame_time{blocks} = tmp_first_last_frame_time;
            
%             figure
%                 hold on
%                 plot(zscore(tmp_dlc_length));
%                 plot(zscore(tmp_tdt_length));
%                 title(['Day ' num2str(day) ' Block ' num2str(blocks) ' - dlc files: ' num2str(length(files)) ' - wave starts: ' num2str(length(trial_starts))])
        end

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['SLDD2_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
 
end

%% T398

clc; clear all; close all;

spike_path = {['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\2.26.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\2.27.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\2.28.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.1.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.2.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.3.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.4.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.5.19'], ...
    ['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\kilosort\3.6.19']};

% BANKS
bank_1_2 = {[5 6 7 8 9 21 39 40 41 42 43 46 51 52 53 54 56 57 58 59 60 61 100 102 105 107 108 110 111 113 114 115 116 119 120 122 123], ...
    [129 130 133 140 141 153 159 160 166 168 170 173 176 179 180 181 187 188 194 195 196 197 198 199 200 201 202 203 205 0206 207 211 212 214 223 224 231]};
bank_2_1 = {[1 2 3 5 13 19 25 31 32 38 40 42 44 45 48 51 52 53 60 67 68 69 70 71 72 73 74 75 77 78 79 83 84 86 95 96 103], ...
    [135 136 137 139 140 149 167 168 169 170 171 179 180 181 182 184 185 186 187 188 189 228 233 235 236 238 239 241 242 243 244 247 248 249 250 251 255]};

% BANK ORDER (M1 then DLS)
bank_order = {[1 2], [1 2], [1 2], [2 1], [1 2], [2 1], [1 2], [1 2], [2 1]};

pre_sleep = {'T398-190226-112203','T398-190227-113133','T398-190228-105336','T398-190301-112846','T398-190302-131102','T398-190303-140804','T398-190304-114807','T398-190305-130625','T398-190306-091252'};
reach_blocks = {'T398-190226-151524','T398-190227-144600','T398-190228-125951','T398-190301-135817','T398-190302-151642','T398-190303-161420','T398-190304-135101','T398-190305-151059','T398-190306-111957'};
post_sleep = {'T398-190226-160611','T398-190227-154618','T398-190228-141352','T398-190301-151942','T398-190302-162611','T398-190303-171207','T398-190304-145843','T398-190305-161031','T398-190306-124648'};

save_path = '\\MyCloudPR4100\Public\SL_Backups\MEGA\data\T398\neural_data_reach';
make_plot = 0;

for day = [1:5 7:9]

    %%% LOAD SPIKING DATA

        opts = delimitedTextImportOptions("NumVariables", 2);
        opts.DataLines = [2, Inf];
        opts.Delimiter = "\t";
        opts.VariableNames = ["cluster_id", "group"];
        opts.VariableTypes = ["double", "categorical"];
        opts = setvaropts(opts, 2, "EmptyFieldRule", "auto");
        opts.ExtraColumnsRule = "ignore";
        opts.EmptyLineRule = "read";
        
        cd([spike_path{day} '\M1']);
        m1_ts = readNPY('spike_times.npy');
        m1_id = readNPY('spike_clusters.npy');
        m1_unit_class = readtable('cluster_group.tsv', opts);
        m1_clust_id = zeros(1,max(m1_unit_class{:,1}));
        m1_clust_id_idx = m1_unit_class{:,1}'+1;
        idx_count = 1;
        for clust_id = m1_clust_id_idx
            tmp_qual = cellstr(m1_unit_class{idx_count,2});
            if strcmp(tmp_qual{1,1},'noise')
                m1_clust_id(clust_id) = 0;
            else
                m1_clust_id(clust_id) = 1;
            end
            idx_count = idx_count + 1;
        end
        m1_clust_id = find(m1_clust_id)-1;

        cd([spike_path{day} '\DLS']);
        dls_ts = readNPY('spike_times.npy');
        dls_id = readNPY('spike_clusters.npy');
        dls_unit_class = readtable('cluster_group.tsv', opts);
        dls_clust_id = zeros(1,max(dls_unit_class{:,1}));
        dls_clust_id_idx = dls_unit_class{:,1}'+1;
        idx_count = 1;
        for clust_id = dls_clust_id_idx
            tmp_qual = cellstr(dls_unit_class{idx_count,2});
            if strcmp(tmp_qual{1,1},'noise')
                dls_clust_id(clust_id) = 0;
            else
                dls_clust_id(clust_id) = 1;
            end
            idx_count = idx_count + 1;
        end
        dls_clust_id = find(dls_clust_id)-1;

        %%% CHECK UNIT %%%
        pre_sleep_length = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' pre_sleep{day}],'CHANNEL',1);
        samp_freq = pre_sleep_length.RSn1.fs;
        pre_sleep_length = length(pre_sleep_length.RSn1.data);

        reach_length = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' reach_blocks{day}],'CHANNEL',1);
        reach_length = length(reach_length.RSn1.data);

        post_sleep_length = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' post_sleep{day}],'CHANNEL',1);
        post_sleep_length = length(post_sleep_length.RSn1.data);

        max(m1_ts)
        max(dls_ts)
        pre_sleep_length + reach_length + post_sleep_length
        %%% CHECK UNIT %%%
        
        m1_spiking = cell(1,length(m1_clust_id),1);
        m1_count = 1;
        for m1_n = m1_clust_id
            tmp_spikes = m1_ts(m1_id==m1_n);
            m1_spiking{1,m1_count,1} = double(tmp_spikes(tmp_spikes>pre_sleep_length & tmp_spikes<pre_sleep_length+reach_length)-pre_sleep_length)/samp_freq;
            m1_count = m1_count+1;
        end
        
        dls_spiking = cell(1,length(dls_clust_id),1);
        dls_count = 1;
        for dls_n = dls_clust_id
            tmp_spikes = dls_ts(dls_id==dls_n);
            dls_spiking{1,dls_count,1} = double(tmp_spikes(tmp_spikes>pre_sleep_length & tmp_spikes<pre_sleep_length+reach_length)-pre_sleep_length)/samp_freq;
            dls_count = dls_count+1;
        end

    %%% LOAD LFP DATA

        if bank_order{day}(1)==1
            M1_bank = bank_1_2{1};
            M1_offset = 1;
            DLS_bank = bank_1_2{2};
            DLS_offset = 129;
        elseif bank_order{day}(1)==2
            M1_bank = bank_2_1{2};
            M1_offset = 129;
            DLS_bank = bank_2_1{1};
            DLS_offset = 1;
        end

        M1_mapping = {(M1_offset+[15 16 14 17 13 18 12 19]), ...%01
            (M1_offset+[11 20 10 21 9 22 8 23]), ...            %02
            [], ...                                             %03
            [], ...                                             %04
            (M1_offset+[63 32 62 33 61 34 60 35]), ...          %05
            (M1_offset+[59 36 58 37 57 38 56 39]), ...          %06
            [], ...                                             %07
            [], ...                                             %08
            (M1_offset+[79 80 78 81 77 82 76 83]), ...          %09
            (M1_offset+[75 84 74 85 73 86 72 87]), ...          %10
            [], ...                                             %11
            [], ...                                             %12
            (M1_offset+[127 96 126 97 125 98 124 99]), ...      %13
            (M1_offset+[123 100 122 101 121 102 120 103]), ...  %14
            [], ...                                             %15
            [], ...                                             %16
            (M1_offset+[111 112 110 113 109 114 108 115]), ...  %17
            (M1_offset+[107 116 106 117 105 118 104 119]), ...  %18
            [], ...                                             %19
            [], ...                                             %20
            (M1_offset+[95 64 94 65 93 66 92 67]), ...          %21
            (M1_offset+[91 68 90 69 89 70 88 71]), ...          %22
            [], ...                                             %23
            [], ...                                             %24
            (M1_offset+[47 48 46 49 45 50 44 51]), ...          %25
            (M1_offset+[43 52 42 53 41 54 40 55]), ...          %26
            [], ...                                             %27
            [], ...                                             %28
            (M1_offset+[31 0 30 1 29 2 28 3]), ...              %29
            (M1_offset+[27 4 26 5 25 6 24 7]), ...              %30
            [], ...                                             %31
            []};                                                %32

        DLS_mapping = {[], ...                                  %01
            [], ...                                             %02
            [], ...                                             %03
            [], ...                                             %04
            [], ...                                             %05
            [], ...                                             %06
            [], ...                                             %07
            [], ...                                             %08
            [], ...                                             %09
            [], ...                                             %10
            [], ...                                             %11
            [], ...                                             %12
            [], ...                                             %13
            [], ...                                             %14
            [], ...                                             %15
            [], ...                                             %16
            (DLS_offset+[15 16 14 17 13 18 12 19]), ...         %17
            (DLS_offset+[11 20 10 21 9 22 8 23]), ...           %18
            (DLS_offset+[7 24 6 25 5 26 4 27]), ...             %19
            (DLS_offset+[3 28 2 29 1 30 0 31]), ...             %20
            (DLS_offset+[63 32 62 33 61 34 60 35]), ...         %21
            (DLS_offset+[59 36 58 37 57 38 56 39]), ...         %22
            (DLS_offset+[55 40 54 41 53 42 52 43]), ...         %23
            (DLS_offset+[51 44 50 45 49 46 48 47]), ...         %24
            (DLS_offset+[79 80 78 81 77 82 76 83]), ...         %25
            (DLS_offset+[75 84 74 85 73 86 72 87]), ...         %26
            (DLS_offset+[71 88 70 89 69 90 68 91]), ...         %27
            (DLS_offset+[67 92 66 93 65 94 64 95]), ...         %28
            (DLS_offset+[127 96 126 97 125 96 124 95]), ...     %29
            (DLS_offset+[123 100 122 101 121 102 120 103]), ... %30
            (DLS_offset+[119 104 118 105 117 106 116 107]), ... %31
            (DLS_offset+[115 108 114 109 113 110 112 111])};    %32

        m1_lfp = cell(1,37);
        count = 1;
        for shank = [1 2 5 6 9 10 13 14 17 18 21 22 25 26 29 30]
            disp(shank);
            [tmp_chans, ~] = intersect(M1_mapping{shank},M1_bank);
            for electrode = 1:length(tmp_chans)
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' reach_blocks{day}],'CHANNEL',tmp_chans(electrode));
                m1_lfp{count} = downsample(data.RSn1.data,20);
                count = count+1;
            end
        end
        
        dls_lfp = cell(1,38);
        count = 1;
        for shank = 17:32
            disp(shank);
            [tmp_chans, ~] = intersect(DLS_mapping{shank},DLS_bank);
            for electrode = 1:length(tmp_chans)
                data = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\RS4_data\M1-DLS_256Probe-190222-144047\' reach_blocks{day}],'CHANNEL',tmp_chans(electrode));
                dls_lfp{count} = downsample(data.RSn1.data,20);
                count = count+1;
            end
        end
        
        lfp_samp_rate = data.RSn1.fs/20;
        
	%%% LOAD BEHAVIORAL MARKERS
    
        beh_markers = cell(1,1);
        trial_video = cell(1,1);
        timelock_frame = cell(1,1);
        first_last_frame_time = cell(1,1);
        
        files = dir(['\\MyCloudPR4100\Public\SL_Backups\T398\reach_data\day_' num2str(day) '\trajectories\*T398*']);
        [~, i] = sort_nat({files(:).name});
        files = files(i);
                
        %%% FIND PELLET POSITION
        pellet_position = [];
        for file = 1:length(files)
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T398\reach_data\day_' num2str(day) '\trajectories\' files(file).name]);
            tmp_pellet_x = trial_kin(1:25,8);
            tmp_pellet_y = trial_kin(1:25,9);
            tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
            tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
            pellet_location = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
            pellet_position = [pellet_position; pellet_location];
        end
        pellet_position = nanmean(pellet_position);

        frames = TDTbin2mat(['\\MyCloudPR4100\Public\SL_Backups\T398\neural_data\TDT_tanks\M1-DLS_256Probe-190222-144047\' reach_blocks{day}]);
        frames = frames.epocs.PC2_.onset;
        trial_starts = [1; find(diff(frames)>1)+1];
        trial_ends = [find(diff(frames)>1); length(frames)];
        
        if day == 1
            files = files(1:73);
        end
        
        if day == 3
            trial_starts = trial_starts(2:end-1);
            trial_ends = trial_ends(2:end-1);
        end
        
        tmp_tdt_length = [];
        tmp_dlc_length = [];
        tmp_beh_markers = [];
        tmp_lock_frame = [];
        tmp_trial_video = cell(1,length(files));
        tmp_first_last_frame_time = [];
        for file = 1:length(files)
            frame_times = frames(trial_starts(file):trial_ends(file));
            trial_kin = readmatrix(['\\MyCloudPR4100\Public\SL_Backups\T398\reach_data\day_' num2str(day) '\trajectories\' files(file).name]);
            %%% get paw trajectory
            tmp_x = trial_kin(:,2);
            tmp_y = trial_kin(:,3);
            tmp_x(trial_kin(:,4)<.99) = NaN;
            tmp_y(trial_kin(:,4)<.99) = NaN;
            %%% get paw velocity
            diff_x = diff(tmp_x);
            diff_y = diff(tmp_y);
            %%% find frame with paw closest to pellet
            [~,i] = min(sqrt(((tmp_x-pellet_position(1)).^2)+((tmp_y-pellet_position(2)).^2)));
            %%% skip bad trials
            if (i<31 || i>size(trial_kin,1)-31)
                tmp_beh_markers = [tmp_beh_markers nan];
                tmp_lock_frame = [tmp_lock_frame; nan nan nan];
                tmp_first_last_frame_time = [tmp_first_last_frame_time; nan nan];
            else
                tmp_beh_markers = [tmp_beh_markers frame_times(i)];
                tmp_lock_frame = [tmp_lock_frame; 1 i length(frame_times)];
                tmp_first_last_frame_time = [tmp_first_last_frame_time; frame_times(1)-frame_times(i) frame_times(end)-frame_times(i)];
            end
            tmp_dlc_length = [tmp_dlc_length size(trial_kin,1)];
            tmp_tdt_length = [tmp_tdt_length frames(trial_ends(file))-frames(trial_starts(file))];
            tmp_trial_video{file} = files(file).name;
        end
        timelock_frame{1} = tmp_lock_frame;
        beh_markers{1} = tmp_beh_markers;
        trial_video{1} = tmp_trial_video;
        first_last_frame_time{1} = tmp_first_last_frame_time;
            
%         %%% CHECK %%%
%         figure
%         hold on
%         plot(zscore(tmp_dlc_length));
%         plot(zscore(tmp_tdt_length));
%         title(['Day ' num2str(day) ' - dlc files: ' num2str(length(files)) ' - wave starts: ' num2str(length(trial_starts))])
%         %%% CHECK %%%

    [m1_reach_spiking, m1_reach_lfp, dls_reach_spiking, dls_reach_lfp, trial_information] = reach_activity(m1_spiking,dls_spiking,m1_lfp,dls_lfp,lfp_samp_rate,beh_markers,trial_video,timelock_frame,first_last_frame_time,day,make_plot);
    
    %%% SAVE
    cd(save_path);
    save(['T398_day_' num2str(day) '_reach_spiking_LFP.mat'],'m1_reach_spiking','dls_reach_spiking','m1_reach_lfp','dls_reach_lfp','trial_information','-v7.3');
    pause(0.1)
 
end