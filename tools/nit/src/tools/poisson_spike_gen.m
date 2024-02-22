function [spikes] = poisson_spike_gen(time, rate, noise_prob)
%%% *function [spikes] = poisson_spike_gen(time, rate, noise_prob)*
%%% 
%%% ### Description
%%% poisson_spike_gen generates a train of Poisson spikes with defined rate.
%%%
%%% ### Inputs:
%%% - *time*: time array for the spike generation of size *n_timesteps*.
%%% - *rate*: rate for the Poisson process, if a double it is assumed a constant rate, if a vector a time varying rate can be specified. In the latter case *length(rate) = n_timesteps*.
%%% - *noise_prob*: noise probability of Poisson process (probability of a spike to be randomly generated/suppressed.
%%%
%%% ### Outputs:
%%% - *spikes*: spike train corresponding to given inputs.

dt = diff(time); dt = dt(1);
if length(rate)==1
    lambda = rate;
end
spikes = zeros(size(time));
for i=1:length(time)
    t = time(i);
    if length(rate)>1
        lambda = rate(i);
    end
    % generate spikes
    if unifrnd(0,1) <= dt*lambda
        spikes(i) = 1;
    end
    % add noise
    if unifrnd(0,1) < noise_prob
        if spikes(i) == 1
            spikes(i) = 0;
        else
            spikes(i) = 1;
        end
    end
end
end

