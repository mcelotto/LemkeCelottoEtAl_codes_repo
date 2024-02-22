
function ura = NPC_kerncoprnd(copula,data,cases,vine,tim)

%%%%%%% this is from 5.15 and Algorithm 5.3 of the following book
% % % Simulating Copulas: Stochastic Models, Sampling Algorithms and Applications

                                
rng('shuffle', 'v5uniform');

[pts,GRID_u]= NPC_mk_grid(size(copula{1,2}.C_grid,1),vine.METH{1,2});
mag=max(GRID_u(:));
mig=min(GRID_u(:));

d = length(copula);
w = rand(cases,d);
w=(mag-mig)*((w-min(w))./(max(w)-min(w)))+mig;


v = zeros(cases,d,d);
v(:,1,1) = reshape(w(:,1),[cases 1 1]);
for i = 2:d
    v(:,i,i) = reshape(w(:,i),[cases 1 1]);
    for k = (i-1):-1:1
        
        tr=k;
        str=i-k+1;
        method=vine.METH{tr,str};
        
        %workerID = get(getCurrentTask(),'ID');
        %save(['crashData/',num2str(workerID),'kerncopccdfinv.mat'])
        
        v(:,k,i)=NPC_kerncopccdfinv(copula{tr,str},[v(:,k,k) v(:,k+1,i)],method);%copula{tr,str}.fit.F.ccdf([v(:,tr,1) v(:,tr,str)]);
    end
end

u=squeeze(v(:,1,:));


%%%%%% from CDF to samples

for i=1:size(u,2)
    
u(u(:,i)>=max(copula{i,1}.MarginG.p),i)=max(copula{i,1}.MarginG.p)-1e-10.*abs(rand(sum(u(:,i)>=max(copula{i,1}.MarginG.p)),1));
u(u(:,i)<=min(copula{i,1}.MarginG.p),i)=min(copula{i,1}.MarginG.p)+1e-10.*abs(rand(sum(u(:,i)<=min(copula{i,1}.MarginG.p)),1));

if strcmp(vine.families{i,1},'discrete')
    ura(:,i)=interp1(copula{i,1}.MarginG.p,copula{i,1}.MarginG.s,u(:,i),'nearest');
else
    ura(:,i)=interp1(copula{i,1}.MarginG.p,copula{i,1}.MarginG.s,u(:,i),'linear');
end
    
end


end



