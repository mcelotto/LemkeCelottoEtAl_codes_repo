function X_cont = NPC_discrete_to_cont(X)
%%% *function X_cont = NPC_discrete_to_cont(X)*
%%%
%%% ### Description
%%% Turns a vector of discrete values *X* to continuous *X_cont* by adding
%%% uniform random noise to each value as specified in Safaai, et al.
%%% "Information estimation using nonparametric copulas." *Physical Review
%%% E* 98.5 (2018): 053302. Sec. B.1
%%%
%%% ### Inputs:
%%% - *X*: vector of discrete values.
%%%
%%% ### Outputs:
%%% - *X_cont*: continuous version of input *X*, transformed by adding
%%% uniformly distributed (between suitable bounds) random noise to 
%%% each of its components
%%%

assert(~any(logical(rem(X,1))),"Input vector 'X' should be made of discrete values.");
X_cont = nan(size(X));
X_vals = sort(unique(X));
X_diff = diff(X_vals);

for i = 1:length(X_vals)
    idxs = find(X == X_vals(i));
    if i < length(X_vals)
        %X_cont(idxs) = X(idxs) + unifrnd(X_vals(i),X_vals(i+1),size(idxs));
        X_cont(idxs) = X(idxs) + unifrnd(0,X_diff(i),size(idxs));
    else
        %X_cont(idxs) = X(idxs) + unifrnd(X_vals(i),X_vals(i)+1,size(idxs));
        X_cont(idxs) = X(idxs) + unifrnd(0,1,size(idxs));
    end
end
end

