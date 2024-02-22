function [C,Mar]=NPC_kernelcdf(x,y,par)

if par.fit~=0
    u=x*0;
    for j=1:numel(x)
        u(j)=sum(x<=x(j))/(numel(x)+1);
    end
    [~,nn]=unique(u);
    Mar.s=x(nn);
    Mar.p=u(nn);
    s=x(nn);
    pp=u(nn);
else
    
    s=par.s;
    pp=par.p;
    Mar=NaN;
end

C=interp1(s,pp,y,'linear');

C(isnan(C) & y>max(s))=par.max*(1-1e-10*abs(rand(size(C(isnan(C) & y>max(s))))));
C(isnan(C) & y<min(s))=par.min*(1+1e-10*abs(rand(size(C(isnan(C) & y<min(s))))));
C(C>=par.max)=par.max*(1-1e-10*abs(rand(size(C(C>=par.max)))));
C(C<=par.min)=par.min*(1+1e-10*abs(rand(size(C(C<=par.min)))));




