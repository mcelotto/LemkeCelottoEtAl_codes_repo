function MI = analyticalMI(lambda,windowVec, probTolerance1, probTolerance2)
%%% *function MI = analyticalMI(lambda,windowVec,probTolerance)*
%%% 
%%% ### Description
%%% analyticalMI computes analytically the mutual information.
%%%
%%% ### Inputs:
%%% - *lambda*: *nStimuli x timeSteps* matrix including the probability of spike for every time step.
%%% - *windowVec*: *1 x nWindows* array including the length of the windows (in timeSteps) in which the responses are binned.
%%% - *probTolerance1*: tolerance determining which spike words are not significantly contributing to the MI calculation. If > 0 it causes the function to go faster but to underestimate the real value of MI.
%%% - *probTolerance2*: tolerance determining which values of the MI should not be computed. If > 0 it causes the function to go faster but to underestimate the real value of MI.
%%%
%%% ### Outputs:
%%% - *MI*: MI value.

    if any(lambda(:)>1)
        warning('Probabilities should not be higher than 1.')
    end

    nStimuli = size(lambda,1);
    nWindows = length(windowVec);
    
    pS = 1/nStimuli * ones(1,nStimuli);
    binnedProb = zeros([max(windowVec)+1,nStimuli,nWindows]);
    
    idx = zeros(1,max(windowVec)); % revise C compilation with idx = zeros(1,max(windowVec));
    for win = 1:nWindows%-1
        idx = idx(end) + (1:windowVec(win));
        idx2 = idx(1:1:windowVec(win));
        lambda_idx = lambda(:,idx2);
        probMat = zeros(length(idx2)+1,nStimuli);
        for s = 1:nStimuli
            idxTol = find(lambda_idx(s,:)>probTolerance1); 
            probMat(1:length(idxTol)+1,s) = ...
                probsCalculation(lambda_idx(s,idxTol)); %getting the probability of having a given amount of spikes in a window (given the stimulus)
        end
        binnedProb(1:size(probMat,1),1:nStimuli,win) = probMat; %gathering all the data
    end

    sequences = permut(windowVec);
    if size(sequences,1) == 1  % if working with spikeRates (aka univariate data)
        sequences = sequences';
    end
    
    MI = 0;
    parfor i = 1:size(sequences,1)
        prod = ones(1,nStimuli);
        for s = 1:nStimuli
            prod(s) = 1;
            for win = 1:nWindows
                prod(s) = prod(s) * binnedProb(sequences(i,win),s,win); %multiply the probability of the sequence by the probability og having n spikes given stimulus s
                if prod(s) <= probTolerance2
                    prod(s) = 0;
                    continue;
                end
            end
        end
        pR = sum(pS.*prod);
        for s = 1:nStimuli
            addedInfo = log2(prod(s)/pR);                                       %Information in bits 
            if isinf(addedInfo) || isnan(addedInfo)                             %If NaN, convert to Zero
                addedInfo=0;
            end
            MI = MI + pS(s)*prod(s)*addedInfo;
        end
    end
end


function prob = probsCalculation(lambda)
% function prob = probsCalculation(lambda)
% 
% ### Description
% probsCalculation computes analytically the probability of summing a
% given amount of spikes when the stimuli lambda is presented
%
% ### Inputs:
% - lambda: nStimuli x timeSteps matrix including the probability of
% spike for every time step
%
% ### Outputs:
% - prob: size(lambda,2)+1 x size(lambda,1) or (possible values of
% spikes x nStimuli) including the probability of having n spikes (in row
% n) when stimulus s (column s) is presented 

    timeSteps = length(lambda);
    prob = zeros(1,timeSteps+1);
    combinations = 2^timeSteps;
    for i = 1:combinations
        sequence = de2bi(2^(timeSteps)+i-1);
        sequence = logical(sequence(end-1:-1:1));
        idx1 = sum(sequence) + 1;
        temp = zeros(1,length(lambda));
        temp(sequence) = lambda(sequence);
        temp(~sequence) = 1 - lambda(~sequence);
        prob(idx1) = prob(idx1) + prod(temp);
    end

end


function sequences = permut(windows)
% function sequences = permut(windows)
% 
% ### Description
% permut computes the possible combinations of sequences that we can have
% in a given windowed trace
%
% ### Inputs:
% - windows: 1 x nWindows matrix including the length (in timesteps) of
% each window
%
% ### Outputs:
% - sequences: nSequences x nWindows including all the possible
% combinations of spikes for each window

    windows = windows+1;
    nperms = prod(windows);
    sequences = zeros(nperms+1,length(windows));

    for j = length(windows):-1:1 % notice reverse order
        tile_size = prod(windows(j:end));
        if j == length(windows)

            tile = (1:tile_size)';

        else
            tile = [];
            for i = windows(j):-1:1
               tile = [tile; i*ones(prod(windows(j+1:end)),1)];
            end
        end
        sequences(2:end,j) = repmat(tile,nperms/tile_size,1);

    end
    sequences = sequences(2:end,:);
end