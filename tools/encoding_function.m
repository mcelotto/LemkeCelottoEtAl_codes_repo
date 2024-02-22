function X = encoding_function(S, vars_order, delta, plt_function)

valsS = unique(S);
if plt_function
    stairs(0.5:1:numel(valsS)+1,[vars_order,vars_order(end)]*delta)
    ylim([0,max(vars_order*delta)+1])
    xlabel('Stimulus')
    ylabel('Response')
end
    
X = nan(size(S));

count = 0;
for i = valsS
    count = count + 1;
    idxs{i} = find(S == i);
    
    X(idxs{i}) = vars_order(count)*delta;
end

end