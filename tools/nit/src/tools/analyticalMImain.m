%%% ### Description
%%% *script used to compile analyticalMI in C++ with Matlab Coder*

clear all;
clc;

addpath(genpath(pwd)); 

dt = 0.02; t = .2;
time = dt:dt:t;


maxRate = 50; %Hz
stimulus(1,:) = maxRate * ones(size(time));
stimulus(2,:) = 0 * ones(size(time));

windows = 5;
windowVec = windowing(time,windows);


info = analyticalMI(stimulus*dt,windowVec,0,0);
