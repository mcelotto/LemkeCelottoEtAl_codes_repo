% temporal rebinning

function outArray = temporal_rebinning(inArray,newBin,method,timeJump)

% k = movmean range

if nargin < 3
    method = 'sum';
elseif nargin < 4
    timeJump = newBin;
end
% Inputs:
% inArray = nUnits x timePoints x nTrials array
% newBin = new time bin in units of the input time bin
%%
nUnits = size(inArray,1);
oldTimePoints = size(inArray,2);
nTrials = size(inArray,3);

oldTcount = 1;

if strcmp(method,'sum')
    newTimePoints = oldTimePoints/newBin;
    outArray = zeros(nUnits,newTimePoints,nTrials);
    for t = 1:newTimePoints
        outArray(:,t,:) = sum(inArray(:,oldTcount:oldTcount+newBin-1,:),2);
        oldTcount = oldTcount + newBin;
    end
elseif strcmp(method,'mean')
    newTimePoints = oldTimePoints/newBin;
    outArray = zeros(nUnits,newTimePoints,nTrials);
    for t = 1:newTimePoints
        outArray(:,t,:) = mean(inArray(:,oldTcount:oldTcount+newBin-1,:),2);
        oldTcount = oldTcount + newBin;
    end
elseif strcmp(method,'movmean')
    newTimePoints = floor((oldTimePoints-(newBin-timeJump))/timeJump);
    outArray = zeros(nUnits,newTimePoints,nTrials);
    for t = 1:newTimePoints
        outArray(:,t,:) = mean(inArray(:,oldTcount:oldTcount+newBin-1,:),2);
        oldTcount = oldTcount + timeJump;
    end
end

end
