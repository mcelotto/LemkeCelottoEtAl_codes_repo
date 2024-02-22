function [] = plot_instPID_Pooled_Sh_multi_timePoints(clusterParams,params,feat2plot,reach_bins,plot_params,fixedKin)
% same as compute_PID_Pooled but with new cluster stat as in Combrisson et
% al. (2022) used for channels selection

animal_colors = distinguishable_colors(numel(params.animals));
cols = distinguishable_colors(4);

n_lfp = 1;
PIDtime = ((26:225)-175)*10; % time points reference used for PID analyses (0 is time of pellet touch)
early_pairs = 295; late_pairs = 510;
infoMaps_early = nan(size(fixedKin,1),early_pairs,39,51); 
infoMaps_late = nan(size(fixedKin,1),late_pairs,39,51); 

for i = 1:size(fixedKin,1) 
    PIDpath = ['/media/DATA/slemke/InfoFlowProject/processed_data/instPID/100Shuff_recRef/fixedKin_',num2str(fixedKin(i,1)),'_',num2str(fixedKin(i,2)),'/'];
    centerTime(i) = 0.5*(PIDtime(fixedKin(i,1)) + PIDtime(fixedKin(i,2)));
    centerTime(i)
    for fidx = feat2plot
        
        disp(params.instReachFeatures{fidx})
        
        early_info_all = [];
        early_info_Orig_all = [];
        early_info_Sh_all = [];
        early_pairs_per_animal = 0;
        late_info_all = [];
        late_info_Orig_all = [];
        late_info_Sh_all = [];
        late_pairs_per_animal = 0;

        tmpFiles = dir([PIDpath 'PID*']);        
        tmpFiles = tmpFiles([3 4 5 6 7 1 2 8],:);
        
        for animal = 1:size(tmpFiles,1)
            
            load([PIDpath tmpFiles(animal).name])
            
            %%% EARLY
            early_info = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51); % for submitted PID analysis last dim = 51
            early_info_Orig = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            early_info_Sh = zeros(size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            for pairs = 1:size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)
                tmp_info = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared(pairs,:,:));
                tmp_infoSh = squeeze(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).sharedSh(pairs,:,:,:));
                early_info(pairs,:,:) = tmp_info - mean(tmp_infoSh,3);
            end
            early_info_all = cat(1,early_info_all,early_info);
            early_info_Orig_all = cat(1,early_info_Orig_all,early_info_Orig);
            early_info_Sh_all = cat(1,early_info_Sh_all,early_info_Sh);
            early_pairs_per_animal = [early_pairs_per_animal size(PIDout.LFP_early_late{1}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)];

            %%% LATE
            late_info = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            late_info_Orig = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            late_info_Sh = zeros(size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1),39,51);
            for pairs = 1:size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)
                tmp_info = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared(pairs,:,:));
                tmp_infoSh = squeeze(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).sharedSh(pairs,:,:,:));
                late_info(pairs,:,:) = tmp_info - mean(tmp_infoSh,3);
            end
            late_info_all = cat(1,late_info_all,late_info);
            late_info_Orig_all = cat(1,late_info_Orig_all,late_info_Orig);
            late_info_Sh_all = cat(1,late_info_Sh_all,late_info_Sh);
            late_pairs_per_animal = [late_pairs_per_animal size(PIDout.LFP_early_late{2}.(params.lfpFeatures{n_lfp}).(params.instReachFeatures{fidx}).shared,1)];

        end
        
        %% PLOT MEAN INFO OVER TIME AND DELAYS (IMAGESC)
            infoMaps_early(i,:,:,:) = early_info_all;
        	infoMaps_late(i,:,:,:) = late_info_all;
    end
end


% info_M1_DLS_early(i,:,:) = squeeze(mean(early_info_all(:,:,27:51),3));
% info_DLS_M1_early(i,:,:) = squeeze(mean(early_info_all(:,:,1:26),3));
% info_M1_DLS_late(i,:,:) = squeeze(mean(late_info_all(:,:,27:51),3));
% info_DLS_M1_late(i,:,:) = squeeze(mean(late_info_all(:,:,1:26),3));
% z_slices = 11:37; % previsouly 2:12
z_slices = 1:16; 

info_M1_DLS_early = squeeze(max(infoMaps_early(:,:,:,27:51),[],4));
info_DLS_M1_early = squeeze(max(infoMaps_early(:,:,:,1:25),[],4));
info_M1_DLS_late = squeeze(max(infoMaps_late(:,:,:,27:51),[],4));
info_DLS_M1_late = squeeze(max(infoMaps_late(:,:,:,1:25),[],4));

% take maximum over time per pair
peaks_M1_DLS_early = max(info_M1_DLS_early(z_slices,:,reach_bins),[],3);
peaks_DLS_M1_early = max(info_DLS_M1_early(z_slices,:,reach_bins),[],3);
peaks_M1_DLS_late = max(info_M1_DLS_late(z_slices,:,reach_bins),[],3);
peaks_DLS_M1_late = max(info_DLS_M1_late(z_slices,:,reach_bins),[],3);


[~,delays_early] = max(squeeze(max(infoMaps_early,[],3)),[],3); 
[~,delays_late] = max(squeeze(max(infoMaps_late,[],3)),[],3); 
change_pos=[]; change_neg=[];
for i = 1:size(delays_late,1)
    change_pos(i) = sum(delays_late(i,:)>=27 & delays_late(i,:)<=41)/size(delays_late,2)-sum(delays_early(i,:)>=27 & delays_early(i,:)<=41)/size(delays_early,2);
    change_neg(i) = sum(delays_late(i,:)>=6 & delays_late(i,:)<=21)/size(delays_late,2)-sum(delays_early(i,:)>=6 & delays_early(i,:)<=21)/size(delays_early,2);
end
% change_pos is the difference in the fraction of channels having an
% M1->DLS delay early and those having an M1->DLS delay late
figure()
hold on
h(1) = plot(centerTime(z_slices),change_pos(z_slices),'color',cols(1,:));
h(2) = plot(centerTime(z_slices),change_neg(z_slices),'color',cols(2,:));
ylabel('change in fraction of channels pairs from naive to skilled')
legend({'M1->DLS','DLS->M1'})
xlim([centerTime(z_slices(1)),centerTime(z_slices(end))])

% Plot info in the four directions
figure(); 

clear h
subplot(3,1,1)
hold on
h(1) = plot(centerTime(z_slices),mean(peaks_M1_DLS_early,2),'color',cols(1,:));
shadedErrorBar(centerTime(z_slices),mean(peaks_M1_DLS_early,2),std(peaks_M1_DLS_early,[],2)/sqrt(size(peaks_M1_DLS_early,2)),'LineProps',{'color',cols(1,:)},'patchSaturation',0.2)
h(2) = plot(centerTime(z_slices),mean(peaks_DLS_M1_early,2),'color',cols(2,:));
shadedErrorBar(centerTime(z_slices),mean(peaks_DLS_M1_early,2),std(peaks_DLS_M1_early,[],2)/sqrt(size(peaks_DLS_M1_early,2)),'LineProps',{'color',cols(2,:)},'patchSaturation',0.2)
legend([h(1) h(2)],{'M1->DLS','DLS->M1'})
xlim([centerTime(z_slices(1)),centerTime(z_slices(end))])
ylabel('[bits]')
xlabel('time from pellet touch [ms]')
title('Naive')

clear h
subplot(3,1,2)
hold on
h(1) = plot(centerTime(z_slices),mean(peaks_M1_DLS_late,2),'color',cols(1,:));
shadedErrorBar(centerTime(z_slices),mean(peaks_M1_DLS_late,2),std(peaks_M1_DLS_late,[],2)/sqrt(size(peaks_M1_DLS_late,2)),'LineProps',{'color',cols(1,:)},'patchSaturation',0.2)
h(2) = plot(centerTime(z_slices),mean(peaks_DLS_M1_late,2),'color',cols(2,:));
shadedErrorBar(centerTime(z_slices),mean(peaks_DLS_M1_late,2),std(peaks_DLS_M1_late,[],2)/sqrt(size(peaks_DLS_M1_late,2)),'LineProps',{'color',cols(2,:)},'patchSaturation',0.2)
xlim([centerTime(z_slices(1)),centerTime(z_slices(end))])
legend([h(1),h(2)],{'M1->DLS','DLS->M1'})
ylabel('[bits]')
xlabel('[time from pellet touch]')
title('Skilled')

subplot(3,1,3)
hold on
h(1) = plot(centerTime(z_slices),mean(peaks_M1_DLS_early-peaks_DLS_M1_early,2),'color',cols(1,:));
shadedErrorBar(centerTime(z_slices),mean(peaks_M1_DLS_early-peaks_DLS_M1_early,2),std(peaks_M1_DLS_early-peaks_DLS_M1_early,[],2)/sqrt(size(peaks_M1_DLS_late,2)),'LineProps',{'color',cols(1,:)},'patchSaturation',0.2)
h(2) = plot(centerTime(z_slices),mean(peaks_M1_DLS_late-peaks_DLS_M1_late,2),'color',cols(2,:));
shadedErrorBar(centerTime(z_slices),mean(peaks_M1_DLS_late-peaks_DLS_M1_late,2),std(peaks_M1_DLS_late-peaks_DLS_M1_late,[],2)/sqrt(size(peaks_DLS_M1_late,2)),'LineProps',{'color',cols(2,:)},'patchSaturation',0.2)
xlim([centerTime(z_slices(1)),centerTime(z_slices(end))])
legend([h(1),h(2)],{'Naive','Skilled'},'AutoUpdate','off')
yline(0,'k--')
ylabel('[bits]')
xlabel('[time from pellet touch]')
title('M1->DLS - DLS->M1')

sgtitle('time-delay peaks differences')

figure(); 

early_diff = peaks_M1_DLS_early-peaks_DLS_M1_early;
late_diff = peaks_M1_DLS_late-peaks_DLS_M1_late;

h(1) = plot(centerTime(z_slices),mean(early_diff,2),'color',cols(1,:));
hold on
shadedErrorBar(centerTime(z_slices),mean(early_diff,2),std(early_diff,[],2)/sqrt(size(early_diff,2)),'LineProps',{'color',cols(1,:)},'patchSaturation',0.2)
h(2) = plot(centerTime(z_slices),mean(late_diff,2),'color',cols(2,:));
shadedErrorBar(centerTime(z_slices),mean(late_diff,2),std(late_diff,[],2)/sqrt(size(late_diff,2)),'LineProps',{'color',cols(2,:)},'patchSaturation',0.2)
ylabel('M1->DLS - DLS->M1')
yline(0,'k--')
legend([h(1) h(2)], 'Naive','Skilled')
xlim([centerTime(z_slices(1)),centerTime(z_slices(end))])

% figure(); 
% plot(centerTime,mean(mean(info_M1_DLS_early,2),3))
% hold on
% plot(centerTime,mean(mean(info_DLS_M1_early,2),3))
% plot(centerTime,mean(mean(info_M1_DLS_late,2),3))
% plot(centerTime,mean(mean(info_DLS_M1_late,2),3))
% legend('M1->DLS Naive','DLS->M1 Naive','M1->DLS Skilled','DLS->M1 Skilled')

figure(); 
plot(centerTime,max(mean(info_M1_DLS_early,2),[],3)-max(mean(info_M1_DLS_late,2),[],3),'color',cols(3,:))
hold on
plot(centerTime,max(mean(info_DLS_M1_early,2),[],3)-max(mean(info_DLS_M1_late,2),[],3),'color',cols(4,:))
yline(0,'k--')
legend('M1->DLS','DLS->M1')
ylabel('Naive - Skilled')

%% Slice plot 
x_neural_time = linspace(-1,0.5,size(infoMaps_early,3));
y_delay = -250:10:250;
z_trajectory_time = centerTime;
z_slices = 5:2:15; % previsouly 2:2:12

tmp_early = squeeze(mean(infoMaps_early,2));
tmp_early = permute(tmp_early,[2,3,1]);
tmp_late = squeeze(mean(infoMaps_late,2));
tmp_late = permute(tmp_late,[2,3,1]);


figure('Position',[2.7,357.7,1274.3,260.3]); 

% Slice plot Naive

ax(1) = subplot(1,2,1);
hs = slice(tmp_early,[],[],z_slices);
colorbar()
caxis([0, 0.012]);

xlim([1,51])
xticks(1:10:51)
xticklabels(-250:100:250) 
xlabel('inter-area delays [ms]')

ylim([11 29])
yticks([11 16 21 26])
yticklabels([-.75 -.5 -.25 0])
ylabel('M1 time [s]')

zticks(z_slices)
zticklabels(centerTime(z_slices))
zlabel('kinematic time [ms]')

hold on
for zi = z_slices
    plot3([25.5 25.5],[0.5 39.5],[zi zi],'r--')
    plot3([26.5 26.5],[0.5 39.5],[zi zi],'r--')
    plot3([.5 51.5],[26 26],[zi zi],'k')
end

shading interp
set(hs,'FaceAlpha',0.8);
view(-52,9.4)

title('Naive')

% Slice plot Skilled

ax(2) = subplot(1,2,2)
hs = slice(tmp_late,[],[],z_slices);
colorbar()
caxis([0, 0.012]);

xlim([1,51])
xticks(1:10:51)
xticklabels(-250:100:250) 
xlabel('inter-area delays [ms]')

ylim([11 29])
yticks([11 16 21 26])
yticklabels([-.75 -.5 -.25 0])
ylabel('M1 time [s]')

zticks(z_slices)
zticklabels(centerTime(z_slices))
zlabel('kinematic time [ms]') 

hold on
for zi = z_slices
    plot3([25.5 25.5],[0.5 39.5],[zi zi],'r--')
    plot3([26.5 26.5],[0.5 39.5],[zi zi],'r--')
    plot3([.5 51.5],[26 26],[zi zi],'k')
end

shading interp
set(hs,'FaceAlpha',0.8);
view(-52,9.4)
title('Skilled')

% Slice plot diff
% 
% ax(3) = subplot(1,3,3);
% hs = slice(tmp_late-tmp_early,[],[],z_slices);
% colormap(ax(3),redblue)
% colorbar()
% caxis([-0.006, 0.006])
% 
% xlim([1,51])
% xticks(1:10:51)
% xticklabels(-250:100:250) 
% xlabel('inter-area delays [ms]')
% 
% ylim([11 29])
% yticks([11 16 21 26])
% yticklabels([-.75 -.5 -.25 0])
% % yticks([6 16 26 36])
% % yticklabels([-1 -.5 0 .5])
% ylabel('M1 time [s]')
% 
% zticks(z_slices)
% zticklabels(centerTime(z_slices))
% zlabel('kinematic time [ms]')
% 
% hold on
% for zi = z_slices
%     plot3([25.5 25.5],[0.5 39.5],[zi zi],'r--')
%     plot3([26.5 26.5],[0.5 39.5],[zi zi],'r--')
%     plot3([.5 51.5],[26 26],[zi zi],'k')
% end
% 
% shading interp
% set(hs,'FaceAlpha',0.8);
% view(-52,9.4)
% 
% title('Skilled-Naive')

end

function c = redblue(m)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009

if nargin < 1, m = size(get(gcf,'colormap'),1); end

if (mod(m,2) == 0)
    % From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
    m1 = m*0.5;
    r = (0:m1-1)'/max(m1-1,1);
    g = r;
    r = [r; ones(m1,1)];
    g = [g; flipud(g)];
    b = flipud(r);
else
    % From [0 0 1] to [1 1 1] to [1 0 0];
    m1 = floor(m*0.5);
    r = (0:m1-1)'/max(m1,1);
    g = r;
    r = [r; ones(m1+1,1)];
    g = [g; 1; flipud(g)];
    b = flipud(r);
end

c = [r g b]; 

end
