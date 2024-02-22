%
%  This source code is part of:
%  NIT - Neuroscience Information Toolbox
%  Copyright (C) 2020  Roberto Maffulli, Miguel Angel Casal Santiago
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%%% # TE binned methods tutorials
%%% ## Tutorial 1
%%% Tutorial 1 shows how to use NIT to calculate Transfer Entropy with 1D rate coding
%%
clear all
close all
clc
rng('default')

% generate rate time traces of two shifted Poisson spiking neurons and plot
t = linspace(0,5,500);
f = 3;
time_lag = 0.1;     % time lag between X and Y in seconds
rate_x = max(0,200*sin(2*pi*f*t));
rate_y = max(0,200*sin(2*pi*f*(t - time_lag)));
frames_lag = round(time_lag/t(2));

figure()
subplot(2,1,1)
title("Firing rates")
hold on
plot(t,rate_x,'Color',[0.6350, 0.0780, 0.1840],'Linewidth',2,'DisplayName','Firing rate X');
plot(t,rate_y,'Color',[0, 0.75, 0.75],'Linewidth',2,'DisplayName','Firing rate Y');
xlabel('t [s]')
ylabel('Firing rate [Hz]')
set(gca,'TickDir','out')
legend()

% generate spike trains for trials and plot
n_trials = 200;
for i=1:n_trials
    spikes_x(:,i) = poisson_spike_gen(t,rate_x,0);
    spikes_y(:,i) = poisson_spike_gen(t,rate_y,0);
end

subplot(2,1,2)
title("Raster Plot")
hold on
for i=1:n_trials
    plot(t(spikes_x(:,i) == 1), i*ones(size(t(spikes_x(:,i) == 1))), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize',3,'MarkerFaceColor',[0.6350, 0.0780, 0.1840],'MarkerEdgeColor',[0.6350, 0.0780, 0.1840]);
    plot(t(spikes_y(:,i) == 1), n_trials + i*ones(size(t(spikes_y(:,i) == 1))), 'LineStyle', 'none', 'Marker', '.', 'MarkerSize',3,'MarkerFaceColor',[0, 0.75, 0.75],'MarkerEdgeColor',[0, 0.75, 0.75]);
end
xlabel('t [s]')
set(gca,'ytick',[])
set(gca,'TickDir','out')

% set transfer entropy calculation
opts.method = 'dr';                                             % direct plug-in probability estimate
opts.bias   = 'naive';                                          % quadratic extrapolation bias correction
taus = linspace(1e-2,0.2,20);                                   % time lag between X and Y in seconds
taus_frames = round(taus/t(2));

% calculate both Transfer Entropy and Normalized Transfer Entropy
for i=1:length(taus_frames)
    opts.taux   = [-taus_frames(i)];        % time lags for x, dimensions are number of time-steps (should always be NEGATIVE!)
    opts.tauy   = [-taus_frames(i)];        % time lags for y, dimensions are number of time-steps (should always be NEGATIVE!)
    outputs = transferentropy(spikes_x, spikes_y, opts, {'TE', 'NTE'});
    TE(i) = outputs{1};
    NTE(i) = outputs{2};
end

% plot TE as a function of the time lags considered
figure()
title("Normalized Transfer Entropy values Vs time lags")
plot(taus,NTE,'Linewidth',2,'Color',[0.6350, 0.0780, 0.1840]);
xlabel('Time lags [s]')
ylabel('Normalized Transfer Entropy [bits]')
set(gca,'TickDir','out')


