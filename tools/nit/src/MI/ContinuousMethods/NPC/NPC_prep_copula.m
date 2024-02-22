function [vine,X_new,ord]=NPC_prep_copula(X,families,range)

% X  t x d matrix, with the first column being the neuron
% margins    1 x d     in this version can be only 'kernel' but easily we
% add other parmeteric

% families 'cont' for continous data or discrete


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d= length(families);
vine.margins = cell(d,1);

vine.families = cell(d);
vine.theta = cell(d);
ord=1:size(X,2);
X_new=X(:,ord);


for i=1:d
    for j=i:d
        if i==j
            vine.families{i,j} = families{1,j};
            vine.families{j,i} = families{1,j};
        end
        vine.METH{i,j}=[1 1];
    end
end



for dd=1:d
    
    RA=rand(size(X_new(:,1)));
    
    vine.margins{dd}.dist = 'kernel';
    vine.margins{dd}.theta = [];
    vine.margins{dd}.iscont = 1;
    
    if strcmp(families{1,dd},'cont')
        
        [nn xx]=unique(X_new(:,dd));
        a=sort(unique(nn));ad=diff(a);ad(numel(ad)+1)=ad(end);Di=X_new(:,dd);
        for i=a'
            Di(X_new(:,dd)==i)=ad(find(a==i));
        end
        vine.margins{dd}.ker = X_new(:,dd)+Di.*rand(size(X_new(:,dd)))*1e-10;
        
    elseif strcmp(families{1,dd},'discrete')
        
        [nn xx]=unique(X_new(:,dd));
        a=sort(unique(nn));ad=diff(a);ad(numel(ad)+1)=1;Di=X_new(:,dd);
        for i=a'
            Di(X_new(:,dd)==i)=ad(find(a==i));
        end
        vine.margins{dd}.ker = X_new(:,dd)+Di.* RA;
    end
    vine.theta{dd,dd} = vine.margins{dd}.ker;  
end
vine.range=range;

