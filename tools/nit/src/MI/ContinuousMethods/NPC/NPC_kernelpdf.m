function [p,Mar]=NPC_kernelpdf(x,TYPE,y,par)

if par.fit~=0
if TYPE==2
[~,pp,s,~]=kde(x,128,min([x]),max([x])+eps);
else
[~,pp,s,~]=kde(x,128,min([x])-eps,max([x])+eps);
end
dX=mean(abs(diff(s)));
pp=pp/sum(pp*dX);
Mar.s=s;
Mar.p=pp;
else
s=par.s;
pp=par.p;
Mar=NaN;
end
y(abs(y-max(s))<1e-9)=max(s)-eps;
y(abs(y-min(s))<1e-9)=min(s)+eps;
p=interp1(s,pp,y,'nearest');
p(y>max(s) | y<min(s))=0;

return 