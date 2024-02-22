function [reach_data] = get_reaching_data_withSuccess(animal,day,path,make_plot,video_direction)
% function to pull reaching data for specific day/animal
%
% outputs:
%
%   reach_data = 1 x trial # struct with following fields: 
%
%       reach_data.position = 2 x 10000 matrix representing X (dim. 1) 
%       and Y (dimension 2) paw position over time (bin = 1ms, locked 
%       to "pellet touch" (same as reaching neural data)
%
%       reach_data.velocity = same size and organization as paw 
%       position, representing paw velocity
%
%       reach_data.acceleration = same size and organization as paw 
%       position, representing paw acceleration

%%% LOAD DATA 

    load([path '\' animal '\neural_data_reach\' animal '_day_' num2str(day) '_reach_spiking_LFP.mat'],'trial_information');
    trial_success = load([path '\' animal '\success_rate\success_rate_day_' num2str(day) '.mat']);
    
%%% GET REACH DATA

    for trial = 1:length(trial_information)

        % load trial kinematics
        
            trial_kin = readmatrix([path '\' animal '\reach_trajectories\day_' num2str(day) '\' trial_information(trial).trial_video]);
            if size(trial_kin,1)<length(trial_information(trial).timelock_frame(1):trial_information(trial).timelock_frame(3))
            else
                trial_kin = trial_kin(trial_information(trial).timelock_frame(1):trial_information(trial).timelock_frame(3),:);
            end
            
        % load trial success rate
        

            
        % set kinematics relative to pellet position and intuitive sign (positive change = towards pellet (x) /up (y)

            if video_direction == 'l'
                tmp_pellet_x = trial_kin(1:25,5);
                tmp_pellet_y = trial_kin(1:25,6);
                tmp_pellet_x(trial_kin(1:25,7)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,7)<.99) = NaN;
                pellet_position = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                trial_kin(:,2) = (trial_kin(:,2)-pellet_position(1))*-1;
                trial_kin(:,3) = (trial_kin(:,3)-pellet_position(2))*-1;
            elseif video_direction == 'r'
                tmp_pellet_x = trial_kin(1:25,8);
                tmp_pellet_y = trial_kin(1:25,9);
                tmp_pellet_x(trial_kin(1:25,10)<.99) = NaN;
                tmp_pellet_y(trial_kin(1:25,10)<.99) = NaN;
                pellet_position = [nanmean(tmp_pellet_x) nanmean(tmp_pellet_y)];
                trial_kin(:,2) = (trial_kin(:,2)-pellet_position(1));
                trial_kin(:,3) = (trial_kin(:,3)-pellet_position(2))*-1;
            end

        % remove low probability tracking
        
            trial_kin(trial_kin(:,4)<.99,2)= NaN;
            trial_kin(trial_kin(:,4)<.99,3)= NaN;

        % get time (relative to pellet touch) for every kinematic bin
        
            if size(trial_kin,1)<length(trial_information(trial).timelock_frame(1):trial_information(trial).timelock_frame(3))
                kin_times = interp1([1 size(trial_kin,1)], ...
                            [trial_information(trial).first_last_frame_time(1) trial_information(trial).first_last_frame_time(2)], ...
                            [1:size(trial_kin,1)]);
            else
                kin_times = interp1([trial_information(trial).timelock_frame(1) trial_information(trial).timelock_frame(3)], ...
                            [trial_information(trial).first_last_frame_time(1) trial_information(trial).first_last_frame_time(2)], ...
                            [trial_information(trial).timelock_frame(1):trial_information(trial).timelock_frame(3)]);
            end

        % interpolate kinematic bin TIMES to 1ms bins
        
            kin_times_1ms = [trial_information(trial).first_last_frame_time(1):0.001:trial_information(trial).first_last_frame_time(2)];    

        % interpolate kinematics to 1ms bins
        
            trial_kin_1ms_x = interp1(kin_times, trial_kin(:,2), kin_times_1ms);
            trial_kin_1ms_y = interp1(kin_times, trial_kin(:,3), kin_times_1ms);

        % "fit" kinematics (binned at 1ms) to -5s to +5s pellet touch window

            % generate NaN arrays
            kin_timelock_x = nan(1,10000);
            kin_timelock_y = nan(1,10000);

            % find index to fit kinematics to -5s to +5s window
            index = round(kin_times_1ms*1000)+5000;

            % remove out-of-window kinematics
            trial_kin_1ms_x(index<1) = [];
            trial_kin_1ms_y(index<1) = [];
            index(index<1) = [];
            
            trial_kin_1ms_x(index>10000) = [];
            trial_kin_1ms_y(index>10000) = [];
            index(index>10000) = [];

            kin_timelock_x(index) = trial_kin_1ms_x;
            kin_timelock_y(index) = trial_kin_1ms_y;

        % define output
        
            reach_data(trial).position = [kin_timelock_x; kin_timelock_y];
            reach_data(trial).velocity = [nan(1,1) diff(kin_timelock_x); nan(1,1) diff(kin_timelock_y)];
            reach_data(trial).acceleration = [nan(1,2) diff(diff(kin_timelock_x)); nan(1,2) diff(diff(kin_timelock_y))];
            reach_data(trial).success = trial_success.trial_success(trial);
            
    end

%%% PLOT

    if make_plot==1

        x_pos = [];
        y_pos = [];
        x_vel = [];
        y_vel = [];
        x_acc = [];
        y_acc = [];
        
        figure
        
            for trial = 1:length(reach_data)
                
                subplot(3,2,1); hold on;
                    plot(reach_data(trial).position(1,:),'color',[.5 .5 .5])
                    x_pos = [x_pos; reach_data(trial).position(1,:)];
                
                subplot(3,2,2); hold on;
                    plot(reach_data(trial).position(2,:),'color',[.5 .5 .5])
                    y_pos = [y_pos; reach_data(trial).position(2,:)];

                subplot(3,2,3); hold on;
                    plot(reach_data(trial).velocity(1,:),'color',[.5 .5 .5])
                    x_vel = [x_vel; reach_data(trial).velocity(1,:)];
                    
                subplot(3,2,4); hold on;
                    plot(reach_data(trial).velocity(2,:),'color',[.5 .5 .5])
                    y_vel = [y_vel; reach_data(trial).velocity(2,:)];
                    
                subplot(3,2,5); hold on;
                    plot(reach_data(trial).acceleration(1,:),'color',[.5 .5 .5])
                    x_acc = [x_acc; reach_data(trial).acceleration(1,:)];
                    
                subplot(3,2,6); hold on;
                    plot(reach_data(trial).acceleration(2,:),'color',[.5 .5 .5])
                    y_acc = [y_acc; reach_data(trial).acceleration(2,:)];
                    
            end    

            subplot(3,2,1); hold on;
                plot(nanmean(x_pos),'k','LineWidth',2)
                title('x position')

            subplot(3,2,2); hold on;
                plot(nanmean(y_pos),'k','LineWidth',2)
                title('y position')

            subplot(3,2,3); hold on;
                plot(nanmean(x_vel),'k','LineWidth',2)
                title('x velocity')

            subplot(3,2,4); hold on;
                plot(nanmean(y_vel),'k','LineWidth',2)
                title('y velocity')

            subplot(3,2,5); hold on;
                plot(nanmean(x_acc),'k','LineWidth',2)
                title('x acceleration')

            subplot(3,2,6); hold on;
                plot(nanmean(y_acc),'k','LineWidth',2)
                title('y acceleration')

    end

end
    