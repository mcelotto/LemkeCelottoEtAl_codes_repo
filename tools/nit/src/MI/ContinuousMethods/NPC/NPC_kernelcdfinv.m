function x=kernelcdfinv(C,par)



if par.fit~=0

    error('fit should be zero')

else
    
s=par.s;
pp=par.p;
Mar=NaN; 
end

% C=interp1(s,pp,y,'nearest');
x=interp1(pp,s,C,'linear');




return

% 
% x=sort(x);
% xx=cat(2,x,x);
% TYPEE=[TYPE TYPE];
% [bandwidth,density,X,Y]=kde2d(xx,N,TYPEE);
% density=density'/sum(sum(density));
% p1=sum(density);
% p1=p1/sum(p1);
% p2=cumsum(p1);
% 
% 
% U=unique(y);
% p=zeros(size(y));
% for i=1:numel(U)
% [~,j]=min(abs(U(i)-p2));    
% p(y==U(i))=X(1,j);
% end