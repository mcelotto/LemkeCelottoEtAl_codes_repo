function PiXY = computePiXY( pxyhys )

% This function computes the partial information on the {X}{Y} node of a
% trivariate lattice with target variable T and source variables (X,Y,hY).
% The I_min definition of redundant information provided by Williams
% and Beer 'Nonnegative decomposition of multivariate information' (2010) 
% is used here since, to the best of our knowledge, it is the only one to 
% be nonnegative on the trivariate lattice.

% Inputs:
% pxyhys = discrete probability distribution of 4 variables. X has to be 
% stored in the first dimension, Y in the second, hY in the third and S in 
% the fourth.

%%%
% In the coming code abbreviations in the names of variables are as
% follows:

% r -> denotes one set of synergy. Basically stands for {} but they
% cannot be in a name of a variable. rX -> {X}, rXrYhY -> {X}{Y, hY}
% X -> past of the sender
% Y -> current of the receiver
% hY -> past of the receiver
% ps -> probability of S - P(S)
% psa -> P(S|A)
% pas -> P(A|S)
% PiR -> greek letter Pi with index R as used in the paper


% Check whether the matrix is a valid probability distribution matrix (sum = 1)
% assert(abs(sum(pxyhys(:)) - 1) > eps, 'The matrix is not a valid joint probability distribution');

% Compute specific surprises fro all subsets of {X, Y, hY}

specSurpriseX = zeros(size(pxyhys, 4), 1);
specSurpriseY = zeros(size(pxyhys, 4), 1);
specSurprisehY = zeros(size(pxyhys, 4), 1);

% For all possible stimuli compute all these matrices
for s = 1:size(pxyhys, 4)
    
    ps = sum(sum(sum(pxyhys(:, :, :, s)))); % probability of S
            
    % {x}
    for x = 1:size(pxyhys, 1)
        psa = sum(sum(sum(sum(pxyhys(x, :, :, s))))) / (sum(sum(sum(sum(pxyhys(x, :, :, :))))) + eps);
        pas = sum(sum(sum(sum(pxyhys(x, :, :, s))))) / (sum(sum(sum(sum(pxyhys(:, :, :, s))))) + eps);
        specSurpriseX(s) = specSurpriseX(s) + pas * (log2(1/(ps + eps)) - log2(1/(psa + eps)));
    end
            
    % {y}
    for y = 1:size(pxyhys, 2)
        psa = sum(sum(sum(sum(pxyhys(:, y, :, s))))) / (sum(sum(sum(sum(pxyhys(:, y, :, :))))) + eps);
        pas = sum(sum(sum(sum(pxyhys(:, y, :, s))))) / (sum(sum(sum(sum(pxyhys(:, :, :, s))))) + eps);
        specSurpriseY(s) = specSurpriseY(s) + pas * (log2(1/(ps + eps)) - log2(1/(psa + eps)));
    end
            
    % {hy}
    for hy = 1:size(pxyhys, 3)
        psa = sum(sum(sum(sum(pxyhys(:, :, hy, s))))) / (sum(sum(sum(sum(pxyhys(:, :, hy, :))))) + eps);
        pas = sum(sum(sum(sum(pxyhys(:, :, hy, s))))) / (sum(sum(sum(sum(pxyhys(:, :, :, s))))) + eps);
        specSurprisehY(s) = specSurprisehY(s) + pas * (log2(1/(ps + eps)) - log2(1/(psa + eps)));
    end
end

% Computation of iMin as defined in the WB paper for all nodes in the lattice
iMinrXrYrhY = 0;
iMinrXrY = 0;


for s = 1:size(pxyhys, 4)
    % {X}{Y}{hY}
    iMinrXrYrhY = iMinrXrYrhY + sum(sum(sum(pxyhys(:, :, :, s)))) * min([specSurpriseX(s) specSurpriseY(s) specSurprisehY(s)]);
    
    % {X}{Y}
    iMinrXrY = iMinrXrY + sum(sum(sum(pxyhys(:, :, :, s)))) * min([specSurpriseX(s) specSurpriseY(s)]);
end

% PiR term for {X}{Y}{hY} is just its iMin
% The next PiR terms are computed as iMin - sum of all preceding PiR terms
% in the lattice. It has to be computed from below up, of course.

PiRrXrYrhY = iMinrXrYrhY;
PiXY = iMinrXrY - PiRrXrYrhY;

end