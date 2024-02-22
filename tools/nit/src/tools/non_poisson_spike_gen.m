function [spikes] = non_poisson_spike_gen(time, rate, noise_prob)
%%% *function [spikes] = non_poisson_spike_gen(time, rate, k)*
%%% 
%%% ### Description
%%% non_poisson_spike_gen generates a train of non-Poisson spikes with
%%% defined rate using a gamma distribution for the ISI
%%%
%%% ### Inputs:
%%% - *time*: time array for the spike generation of size *n_timesteps*.
%%% - *rate*: rate for the Poisson process, must be a single value.
%%% - *k*: shape parameter of the gamma distribution (k<1 super-Poisson k=1 Poisson-like k>1 sub-Poisson).
%%%
%%% ### Outputs:
%%% - *spikes*: spike train corresponding to given inputs.

dt = diff(time); dt = dt(1);
assert(length(rate)==1);
spikes = zeros(1,length(time));

last_spike_time = 0;
mean_isi = 1/rate;
theta = mean_isi/k;
spike_times = [];
while last_spike_time < time(end)
    isi = gamrnd(k,theta);
    last_spike_time = last_spike_time + isi;
    if last_spike_time <= time(end)
        spike_times = [spike_times last_spike_time];
    end
end
% create spike train
for j=1:length(spike_times)
    spikes(round(spike_times(j)/dt)+1) = 1;
end

end

