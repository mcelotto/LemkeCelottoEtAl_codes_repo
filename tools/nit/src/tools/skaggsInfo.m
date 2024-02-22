function [I_bits,I_time,I_spike] = skaggsInfo(data,sVec,time)
%%% *function [MI_bits,MI_time,MI_spike] = skaggsInfo(data,sVec,time)*
%%% 
%%% ### Description
%%% poisson_spike_gen generates a train of Poisson spikes with defined rate.
%%%
%%% ### Inputs:
%%% - *data*: *nTimesteps x nTrials* array including if spikes were produced for each timestep in each trial.
%%% - *sVec*: 1D array including which stimulus is applied in each of the trials *length(rate) = nTrials*.
%%% - *time*: 1 x nTimeSteps time array. Array contains the actual time for each timeStep.
%%%
%%% ### Outputs:
%%% - *MI_bits*: information in bits following Skaggs approach.
%%% - *MI_time*: information in bits/second following Skaggs approach.
%%% - *MI_spike*: information in bits/spike following Skaggs approach.


deltat = time(end);
uniqueS = unique(sVec);
nStimuli = length(uniqueS);

for s = 1:nStimuli
    pS(s) = length(find(sVec==uniqueS(s))) / length(sVec);
end

spikeCount = sum(data,2);
lambda_avg = mean(spikeCount(:))/deltat;

lambda = zeros(nStimuli,1);

for s = 1:nStimuli
    idx = find(sVec==uniqueS(s));
    dataTemp = spikeCount(idx);
    lambda(s) = mean(dataTemp(:))/deltat;
end
clearvars dataTemp

addedInfo = log2(lambda./lambda_avg);
addedInfo(isnan(addedInfo)) = 0;
addedInfo(isinf(addedInfo)) = 0;

I_bits = sum(lambda .* pS' * deltat .* addedInfo,1);
I_time = sum(lambda .* pS' .* addedInfo,1);
I_spike = sum(lambda / lambda_avg .* pS' .* addedInfo,1);

end

