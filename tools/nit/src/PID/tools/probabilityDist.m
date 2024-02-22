function result = probabilityDist(pastX, Y, pastY, stimulus)

%Computes probability distribution among all
%combinations of quadruples pastX(i), Y(i), pastY(i), stimulus(i) 

% All inputs have to be 1 X nTrials arrays of discretized variables
% containing only positive values.
if min(stimulus) <= 0
    stimulus = stimulus - min(stimulus) + 1;
end
if min(pastX) <= 0
    pastX = pastX - min(pastX) + 1;
end
if min(pastY) <= 0
    pastY = pastY - min(pastY) + 1;
end
if min(Y) <= 0
    Y = Y - min(Y) + 1;
end

assert(min(pastX) > 0, 'The pastX has invalid values');
assert(min(Y) > 0, 'The Y has invalid values');
assert(min(pastY) > 0, 'The pastY has invalid values');
assert(min(stimulus) > 0, 'The stimulus has invalid values');

count = length(pastX);

result = accumarray([pastX; Y; pastY; stimulus]', 1, [max(pastX) max(Y) max(pastY) max(stimulus)]) / count;

end

