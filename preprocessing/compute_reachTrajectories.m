%% COMPUTE REACH TRAJECTORIES 
% Pre-compute and generate one .mat file for all animals and days
% containing reach trajectories

    clear all; clc; close all;

%%% INITIALIZE PARAMETERS

    % data paramaters
    animals = {'T102','T107','T110','T200','T201','SLDD1','SLDD2','T398'};
%     days = {[1:8],[2:6],[2:7],[2:10],[4:17],[3:14],[4:13],[3:9]};
    days = {[1:8],[2:6],[2:7],[2:8],[4:17],[3:12],[4:13],[3:9]};
    days{3}(3) = []; days{8}(4) = []; % remove missing day  
    video_direction = {'l','l','l','l','l','r','r','r'}; % in the video, is the rat reaching right (r) or left (l)
    data_path = 'G:\final_IIT\data\animal_data';
    save_path = 'G:\final_IIT\data\preprocessed_data';

	% general paramaters     
    make_plot = 1;
    
    % temporal parameters        
    tMin = 4001;
    tMax = 5000;

%%% COMPUTE

    for animal = 1:numel(animals)
        disp(['Animal ', num2str(animal)])
        dayCount = 1;
        for day = days{animal}

            [reach_data] = get_reaching_data_withSuccess(animals{animal},day,data_path,make_plot,video_direction{animal});
            nTrials(animal,dayCount) = numel(reach_data);

            trajectories{animal}{dayCount} = zeros(6,nTrials(animal,dayCount),tMax-tMin+1); % 1 = Xpos; 2 = Ypos; 3 = Xvel; 4 = Yvel;  5 = Xacc; 6 = Yacc;
            for tidx = 1:nTrials(animal,dayCount)
                trajectories{animal}{dayCount}(1,tidx,:) = reach_data(tidx).position(1,tMin:tMax);
                trajectories{animal}{dayCount}(2,tidx,:) = reach_data(tidx).position(2,tMin:tMax);
                trajectories{animal}{dayCount}(3,tidx,:) = reach_data(tidx).velocity(1,tMin:tMax);
                trajectories{animal}{dayCount}(4,tidx,:) = reach_data(tidx).velocity(2,tMin:tMax);
                trajectories{animal}{dayCount}(5,tidx,:) = reach_data(tidx).acceleration(1,tMin:tMax);
                trajectories{animal}{dayCount}(6,tidx,:) = reach_data(tidx).acceleration(2,tMin:tMax);            
            end
            success_rate{animal}{dayCount} =  [reach_data.success];
            
            dayCount = dayCount + 1;
        end
    end

%     fname=['all_animals_trajectories_' num2str(tMin) '_to_' num2str(tMax) '_' date];
%     save([save_path '\' fname],'trajectories','success_rate','tMin','tMax')
