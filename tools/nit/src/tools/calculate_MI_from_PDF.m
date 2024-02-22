function [MI] = calculate_MI_from_PDF(jointProb)
%%% *function [MI] = calculate_MI_from_PDF(jointProb)*
%%%
%%% ### Description
%%% Computes mutual information given a estimated joint probability distribution.
%%%
%%% ### Inputs:
%%% - *jointProb*: n x m matrix. n and m are the possible states first and second variable can take, respectively.
%%%
%%% ### Outputs:
%%% - *[MI]*: mutual information.

pstimulus = sum(jointProb,2);                                               %probability of X
pspike = sum(jointProb,1);                                                  %probability of Y

MI = 0;                                                                     %Set Normalized Mutual Information
for i = 1:size(jointProb,1)                                                 %Sum over Y
    for j = 1:size(jointProb,2)                                             %Sum over X
        pst = pstimulus(i);
        psp = pspike(j);
        addedInfo = log2(jointProb(i,j)/(pst*psp));                         %Information in bits 
        if isinf(addedInfo) || isnan(addedInfo)                             %If NaN, convert to Zero
            addedInfo=0;
        end
        MI = MI + jointProb(i,j)*addedInfo;                                 %Sum MI in bits 
    end
end

end

