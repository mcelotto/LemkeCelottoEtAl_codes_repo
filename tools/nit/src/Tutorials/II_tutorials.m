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
%%% # II tutorials
%%% ## Tutorial 1
%%% Tutorial 1 shows how to use NIT to calculate Intersection Information with 1D rate coding
clear all
close all
clc

dt = 1/50;
trialendtime = 0.4;
t_trial = 0:dt:trialendtime;
nStimuli = 2;
nTrials = 100;

% generate stimuli
mu = .2;  % mu for the gaussian
sigmaT = [.05,.05];  % sigma for each stimulus
rate = [50 20];  % peak rate for each stimulus
stimuli = nan(nStimuli,length(t_trial));
for i = 1:nStimuli
    signal = normpdf(t_trial,mu,sigmaT(i));
    stimuli(i,:) =  rate(i) * signal / max(signal);  %normalize for peak
end

% generate responses and choice (spike count)
R = [];
S = [];
C = [];
for i = 1:nStimuli
    for j = 1:nTrials
        R = [R sum(poisson_spike_gen(t_trial, stimuli(i,:), 0))];
        S = [S; i];
        C = [C; 0]; % II is expected to be 0 given this choice vector
    end
end

% define options
opts.bin_method = 'eqspace';
opts.bias = 'qe';
opts.draws_per_split_number = 100;
opts.th_bias_corr_convergence_slope = 0.5E-4;
opts.n_bins = 2;
[II_biased, II_unbiased] = II(S, R, C, opts);
