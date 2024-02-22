function [] = plot_mutualInformation_alignComparison(plot_params)

mo_lfp_info = load('MO_LFP_info.mat');
pt_lfp_info = load('PT_LFP_info.mat');

mo_spikes_info = load('MO_Spiking_info.mat');
pt_spikes_info = load('PT_spikes_info.mat');

figure;
    subplot(2,2,1); hold on;
        v_mo = max(mo_lfp_info.m1_early_info');
        v_pt = max(pt_lfp_info.m1_early_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('m1 early lfp info')
    subplot(2,2,2); hold on;
        v_mo = max(mo_lfp_info.dls_early_info');
        v_pt = max(pt_lfp_info.dls_early_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('dls early lfp info')
    subplot(2,2,3); hold on;
        v_mo = max(mo_lfp_info.m1_late_info');
        v_pt = max(pt_lfp_info.m1_late_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('m1 late lfp info')
    subplot(2,2,4); hold on;
        v_mo = max(mo_lfp_info.dls_late_info');
        v_pt = max(pt_lfp_info.dls_late_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('dls late lfp info')
    if plot_params.save == 1
        print([plot_params.save_path '\LFP_mo_pt_align_comparison.eps'],'-painters','-depsc');
    end

figure;
    subplot(2,2,1); hold on;
        v_mo = max(mo_spikes_info.m1_early_info');
        v_pt = max(pt_spikes_info.m1_early_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('m1 early spikes info')
    subplot(2,2,2); hold on;
        v_mo = max(mo_spikes_info.dls_early_info');
        v_pt = max(pt_spikes_info.dls_early_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('dls early spikes info')
    subplot(2,2,3); hold on;
        v_mo = max(mo_spikes_info.m1_late_info');
        v_pt = max(pt_spikes_info.m1_late_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('m1 late spikes info')
    subplot(2,2,4); hold on;
        v_mo = max(mo_spikes_info.dls_late_info');
        v_pt = max(pt_spikes_info.dls_late_info');
        cdfplot(v_mo);
        cdfplot(v_pt);
        [p,~,stats] = ranksum(v_mo,v_pt)
        xlim([0 0.15]);
        grid off;
        title('dls late spikes info')
    if plot_params.save == 1
        print([plot_params.save_path '\SPIKES_mo_pt_align_comparison.eps'],'-painters','-depsc');
    end

end
