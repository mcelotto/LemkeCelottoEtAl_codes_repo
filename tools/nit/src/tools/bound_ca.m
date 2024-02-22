function [bound_ca] = bound_ca(signal, indicator, moving_average_window, ca_rest)
%%% *function [bound_ca] = bound_ca(signal, indicator, moving_average_window)*
%%%
%%% ### Description
%%% bound_ca calculates the moles of bound Ca from the fluorescence signal by inverting Hill's eqn.
%%% The function operates on a filtered version of the input signal, calculated by performing a moving average over *moving_average_window*.
%%%
%%% ### Inputs:
%%% - *signal*: time array for the signal
%%% - *indicator*: fluorescence indicator, must be one of `"GCaMP6f"` or `"GCaMP6s"`
%%% - *moving_average_window*: number of frames used to smooth signal
%%% - *ca_rest*: resting state calcium concentration
%%%
%%% ### Outputs:
%%% - *bound_ca*: time trace of bound Ca
plot_traces = false;

% check inputs
allowed_indicators = ["GCaMP6f", "GCaMP6s"];
assert(any(contains(allowed_indicators,indicator)),...
    "Unknown Ca indicator. Allowed indicators are: "...
    + join(allowed_indicators));
% get Hill's eq params
[dF_F_max, K_d] = get_ca_sat_params(indicator);

% calculate moving average signal
filtered_signal = movmean(signal, moving_average_window);

% cap the signal to max value of dF_F_max as inversion of saturation equation fails
% otherwise. The reason for fluorescence to be, sometimes, higher than dF_F_max
% is because of noise added after generation of dF_F
filtered_signal(filtered_signal >= dF_F_max-0.05) = dF_F_max-0.05;

% remove negative values (at this point, if any, they should only be due to
% noise so it seems plausible to just take the absolute value
filtered_signal = abs(filtered_signal);
bound_ca = invert_sat_eq(filtered_signal,dF_F_max, K_d, ca_rest);

if plot_traces
    cols = brewermap(8, 'Dark2');
    figure()
    hold on
    subplot(3,1,1)
    plot(signal,'DisplayName','F','LineWidth',2,'Color',cols(7,:))
    legend()
    subplot(3,1,2)
    plot(filtered_signal,'DisplayName','Filtered F','LineWidth',2,'Color',cols(4,:))
    legend()
    subplot(3,1,3)
    hold on
    plot(bound_ca,'DisplayName','Ca_b','LineWidth',2,'Color',cols(2,:))
    legend()
    box off
    set(gcf,'color',[1,1,1])  
end
end

function [bound_ca] = invert_sat_eq(filtered_signal,dF_F_max, K_d, ca_rest)
bound_ca = (ca_rest + K_d*filtered_signal/dF_F_max)./(1 - filtered_signal/dF_F_max);
end
