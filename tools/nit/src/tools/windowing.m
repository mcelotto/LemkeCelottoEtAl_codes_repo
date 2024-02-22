function windowVec = windowing(time,nWindows)
%%% *function windowVec = windowing(time,nWindows)*
%%% 
%%% ### Description
%%% windowing chooses the best way of splitting the data given the amount of windows provided by the user. The amount of windows can be slightly modified from the input.
%%%
%%% ### Inputs:
%%% - *time*: vector including time steps.
%%% - *nWindows*: integer including how many windows we want to split the data into.
%%%
%%% ### Outputs:
%%% - *windowVec*: vector inluding the length of each window.

assert(nWindows<=length(time),'nWindows should not be higher than the amount of timesteps');

windowSize = length(time)/nWindows;
windowSizeRounded = floor(windowSize);
windowVec = windowSizeRounded*ones(1,nWindows);

if windowSize - windowSizeRounded > 0.5 %if the remaining time steps are more than half a window, add a new window
    re = (windowSize - windowSizeRounded) * nWindows;
    nWindows = nWindows + 1;
    windowVec = [windowVec,double(uint8(re))];
elseif windowSize - windowSizeRounded <=0.5 %if the remaining time steps are less than half a window, pack them into the last window
    re = (windowSize - windowSizeRounded) * nWindows;
    windowVec(end) = windowVec(end) + uint8(re);%rem(length(time),nWindows);
end