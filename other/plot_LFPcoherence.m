function [] = plot_LFPcoherence(params,coh_params,plot_params)

n_lfp = 1;
fidx = 1;

all_animal_coh = zeros(2,41,19,8);
all_channel_early = [];
all_channel_late = [];

for animal = 1:numel(params.animals)

    animal_neuralActivity = load(['G:\CurrBiol_revision\data\preprocessed_data\neuralActivity_' params.animals{animal} '_10ms_04-Jun-2021.mat']);

    m1_chans = size(animal_neuralActivity.LFP{1}.reach_related_lfp_ref.M1,1);
    dls_chans = size(animal_neuralActivity.LFP{1}.reach_related_lfp_ref.DLS,1);
    
    % early
    early_animal_coh = zeros(2,41,19,m1_chans*dls_chans);            
    pairIdx = 1;
    for m1_chan = 1:m1_chans
        for dls_chan = 1:dls_chans

            tmp_C = zeros(41,19,params.num_earlylate_days{animal});
            dayIdx = 1;
            for day = 1:params.num_earlylate_days{animal}
                early_m1_lfp = squeeze(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1(m1_chan,:,:));
                early_dls_lfp = squeeze(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.DLS(dls_chan,:,:));
                [C,phi,S12,S1,S2,t,f]=cohgramc(early_m1_lfp,early_dls_lfp,coh_params.movingwin,coh_params);
                tmp_C(:,:,dayIdx) = C;
                dayIdx = dayIdx + 1;
            end

            all_channel_early = cat(3,all_channel_early,squeeze(mean(tmp_C,3)));
            early_animal_coh(1,:,:,pairIdx) = squeeze(mean(tmp_C,3));
            pairIdx = pairIdx + 1;
        end
    end

    % late
    late_animal_coh = zeros(2,41,19,m1_chans*dls_chans);            
    pairIdx = 1;
    for m1_chan = 1:m1_chans
        for dls_chan = 1:dls_chans

            tmp_C = zeros(41,19,params.num_earlylate_days{animal});
            dayIdx = 1;
            for day = length(animal_neuralActivity.LFP)-params.num_earlylate_days{animal}+1:length(animal_neuralActivity.LFP)
                late_m1_lfp = squeeze(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.M1(m1_chan,:,:));
                late_dls_lfp = squeeze(animal_neuralActivity.LFP{day}.reach_related_lfp_ref.DLS(dls_chan,:,:));
                [C,phi,S12,S1,S2,t,f]=cohgramc(late_m1_lfp,late_dls_lfp,coh_params.movingwin,coh_params);
                tmp_C(:,:,dayIdx) = C;
                dayIdx = dayIdx + 1;
            end

            all_channel_late = cat(3,all_channel_late,squeeze(mean(tmp_C,3)));
            late_animal_coh(2,:,:,pairIdx) = squeeze(mean(tmp_C,3));
            pairIdx = pairIdx + 1;
        end
    end

    % plot
    figure;
        subplot(1,3,1);
            imagesc(squeeze(mean(early_animal_coh(1,:,:,:),4))')
            set(gca, 'YDir','normal')
            colorbar;
            title(num2str(animal));
        subplot(1,3,2);
            imagesc(squeeze(mean(late_animal_coh(2,:,:,:),4))')
            set(gca, 'YDir','normal')
            colorbar;
        subplot(1,3,3);
            imagesc(squeeze(mean(late_animal_coh(2,:,:,:),4))'-squeeze(mean(early_animal_coh(1,:,:,:),4))')
            set(gca, 'YDir','normal')
            colorbar;

    all_animal_coh(1,:,:,animal) = squeeze(mean(early_animal_coh(1,:,:,:),4));
    all_animal_coh(2,:,:,animal) = squeeze(mean(late_animal_coh(2,:,:,:),4));
    pause(0.1);

end

max_early_LFPCoh = [];
for pair = 1:size(all_channel_early,3)
    max_early_LFPCoh = [max_early_LFPCoh squeeze(mean(mean(all_channel_early(26:36,1:4,pair),1),2))];
end
max_late_LFPCoh = [];
for pair = 1:size(all_channel_late,3)
    max_late_LFPCoh = [max_late_LFPCoh squeeze(mean(mean(all_channel_late(26:36,1:4,pair),1),2))];
end

figure;
    subplot(1,3,1); hold on;
        imagesc(squeeze(mean(all_channel_early(:,:,max_early_LFPCoh>0.2),3))',[.2 .5])
        plot([26 26],[1 19],'color','r');
        set(gca, 'YDir','normal')
        colorbar;
        xlim([5 41])
        ylim([1 19])
        xticks([5:2:36])
        xticklabels([t(5:2:36)])
        yticks([1:2:19])
        yticklabels([f(1:2:19)])
    subplot(1,3,2); hold on;
        imagesc(squeeze(mean(all_channel_late(:,:,max_late_LFPCoh>0.2),3))',[.2 .5])
        plot([26 26],[1 19],'color','r');
        set(gca, 'YDir','normal')
        colorbar;
        xlim([5 41])
        ylim([1 19])
        xticks([5:2:36])
        xticklabels([t(5:2:36)])
        yticks([1:2:19])
        yticklabels([f(1:2:19)])
    subplot(1,3,3); hold on;
        imagesc(squeeze(mean(all_channel_late(:,:,max_late_LFPCoh>0.2),3))'-squeeze(mean(all_channel_early(:,:,max_early_LFPCoh>0.2),3))')
        plot([26 26],[1 19],'color','r');
        set(gca, 'YDir','normal')
        colorbar;
        xlim([5 41])
        ylim([1 19])
        xticks([5:2:36])
        xticklabels([t(5:2:36)])
        yticks([1:2:19])
        yticklabels([f(1:2:19)])
    if plot_params.save == 1
        print([plot_params.save_path '\lfp_coh_imagesc.eps'],'-painters','-depsc');
    end

figure;
    clearvars g
    x = [1:41];
    y = [squeeze(mean(all_channel_early(:,1:4,max_early_LFPCoh>0.2),2))'; squeeze(mean(all_channel_late(:,1:4,max_late_LFPCoh>0.2),2))'];
    c = [ones(1,size(all_channel_early(:,:,max_early_LFPCoh>0.2),3)) 2*ones(1,size(all_channel_late(:,:,max_late_LFPCoh>0.2),3))];
    g(1,1)=gramm('x',x,'y',y,'color',c);
    g(1,1).stat_summary('type','sem');
    g(1,1).axe_property('XLim',[5 41]);
    g(1,1).axe_property('YLim',[.35 .5]);
    g.draw();
    if plot_params.save == 1
        print([plot_params.save_path '\line_plot_imagesc.eps'],'-painters','-depsc');
    end

figure;
    subplot(2,1,1); hold on;
        histogram(max_early_LFPCoh(max_early_LFPCoh>0.2),'Normalization','Probability','DisplayStyle','Stairs');
        histogram(max_late_LFPCoh(max_late_LFPCoh>0.2),'Normalization','Probability','DisplayStyle','Stairs');
        [p,~,stats] = ranksum(max_early_LFPCoh(max_early_LFPCoh>0.2),max_late_LFPCoh(max_late_LFPCoh>0.2));
        title(p);
        xlim([0 1]);
    subplot(2,1,2);
        hold on;
        cdfplot(max_early_LFPCoh(max_early_LFPCoh>0.2));
        cdfplot(max_late_LFPCoh(max_late_LFPCoh>0.2));
        xlim([.2 1]);
    if plot_params.save == 1
        print([plot_params.save_path '\cdf_imagesc.eps'],'-painters','-depsc');
    end
    

        

end