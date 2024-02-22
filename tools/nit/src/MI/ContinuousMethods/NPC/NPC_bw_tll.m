function bw=NPC_bw_tll(data,deg)

n = size(data,1);
bw = 5 * n^(-1 / (4 * deg + 2)) * (chol(cov(data)))';
