function [I_second,I_spike] = skaggs(data,sVec)

stimuli = unique(sVec);

for i = 1:length(stimuli)
    idx = find(sVec == stimuli(i));
    p_i(i) = length(idx) / length(sVec);
    lambda_i(i) = mean(data(idx));
end
lambda = mean(data);

addedInfo = log2(lambda_i./lambda);
addedInfo(isnan(addedInfo)) = 0;
addedInfo(isinf(addedInfo)) = 0;

I_second = sum(p_i .* lambda_i .* addedInfo);
I_spike = I_second / lambda;
if isnan(I_spike), I_spike = 0; end

end
    
    



