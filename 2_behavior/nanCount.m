function [totNan, maxConsNan] = nanCount(v)
% v is a 1 x N array

% totNan is the number of nan values in the array
% maxConsNan is the maximum number of consecutively Nan values

if sum(isnan(v)) > 0
    totNan = sum(isnan(v));
    [nanClus, L] = bwlabel(isnan(v));
    maxConsNan = 0;
    for cidx = 1:L
        cSize = sum(nanClus == cidx);
        maxConsNan = max([cSize,maxConsNan]);
    end
    
else
    totNan = 0;
    maxConsNan = 0;
end