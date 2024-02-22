
function [ccdf_points,pd_points,ccdf_grid,pd_grid,ccdf_data,pd_data,COPULA]=NPC_KernelCop(knots,dat,poi,method,method_fit,Lfit,fitt,parallel)

%%%%%    _u --> uv space
%%%%%    _S --> \Phi space
%%%%%    _X --> PC spacec

n=size(dat,1);
d=size(dat,2);
%%%%%%%%%
data.u=dat;
data.S=norminv(data.u,0,1);    
[COEFF,data.X,~,~,~,mu] = pca(data.S);
data.X=data.S * COEFF - repmat(mu,size(data.S,1),1) * COEFF;

%%%%%%%%
% knots
if fitt==1
    bw=NPC_bw_tll(data.X,2);
    bw(bw<1e-3)=1e-3;
    bw=[bw(1,1) bw(2,2)];
    bw=bw/10;
else
    bw=Lfit.bw;
end
%%%%%%% fit the model
points.u=poi;
points.S=norminv(points.u,0,1);
points.X=points.S * COEFF - repmat(mu,size(points.S,1),1) * COEFF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ischar(method)
    Grid=NPC_GRID_Bands(dat,knots,method);
else
    Grid=NPC_GRID_Bands(dat,knots,method);
end
GRID_X=Grid.X;
GRID_S=Grid.S;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.u=double(data.u);
data.S=double(data.S);
data.X=double(data.X);

points.u=double(points.u);
points.S=double(points.S);
points.X=double(points.X);

[x1,~,~]=unique(sort(Grid.S(:,1)));
xd11=[diff(x1)];xd11=[xd11;xd11(1)];
[x2,~,~]=unique(sort(Grid.S(:,2)));
xd22=[diff(x2)];xd22=[xd22;xd22(1)];
xd1=nan(size(Grid.u,1),1);
xd2=nan(size(Grid.u,1),1);
xd=xd1;
for k=1:numel(x1)
    xd1(Grid.S(:,1)==x1(k))=xd11(k);
    xd2(Grid.S(:,2)==x2(k))=xd22(k);
end

for k1=1:numel(x1)
    for k2=1:numel(x2)
        xd(Grid.S(:,1)==x1(k1) & Grid.S(:,2)==x2(k2))=xd11(k1)*xd22(k2);
    end
end

if numel(uniquetol(xd(:),1e-3))==1
    xd=uniquetol(xd(:),1e-3);
end

P1=normpdf(x1);
P2=normpdf(x2);
NORM=(P1*P2')';

samp=1:size(data.X,1);

afk = cvpartition(numel(samp),'KFold',5);  % Stratified cross-validation
bw0=double(bw);

if fitt==1
    
    if strcmp(method_fit,'LL1') | strcmp(method_fit,'LL2')
        for i=1:numel(afk.TestSize)
            dd{i}.u=data.u(samp(~afk.test(i)==1),:);
            dd{i}.S=data.S(samp(~afk.test(i)==1),:);
            dd{i}.X=data.X(samp(~afk.test(i)==1),:);
            ddT{i}.u=data.u(samp(afk.test(i)==1),:);
            ddT{i}.S=data.S(samp(afk.test(i)==1),:);
            ddT{i}.X=data.X(samp(afk.test(i)==1),:);
        end
        
        options=optimset('MaxIter',100,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',100);
        
        if strcmp(method_fit,'LL2')
            options=optimset('MaxIter',50,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',50);
            bw1=fminsearchbnd(@(x) NPC_MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,1,parallel),bw0,[1e-4 1e-4],[2 2],options);
            options=optimset('MaxIter',50,'Display','off','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50);
            bw1=fminsearchbnd(@(x) NPC_MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,0,parallel),bw0,bw1*0.8,bw1*1.2,options);
            bw=abs(bw1);
            %disp('Opt=LL2')
        end
        
        if strcmp(method_fit,'LL1')
            options=optimset('MaxIter',50,'Display','off','TolX',1e-1,'TolFun',1e-1,'MaxFunEvals',50);
            bw1=fminbnd(@(x) NPC_MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,1,parallel),[1e-4],[2],options);
            options=optimset('MaxIter',50,'Display','off','TolX',1e-2,'TolFun',1e-2,'MaxFunEvals',50);
            bw1=fminbnd(@(x) NPC_MISE(x,data,Grid,points.u(1,:),xd,afk,dd,ddT,bw0,0,NORM,0,parallel),[bw1*0.8],[bw1*1.2],options);
             bw=abs(bw0*bw1);
            %disp('Opt=LL1')
        end
        
    else
        error('Bandwidth method unknown')
    end
    
    lfit = NPC_loclik_fit(bw,data,Grid);
    [pd_data,ccdf_data,pd_grid,ccdf_grid,pd_points,ccdf_points,F,norma]=NPC_func_tll(lfit,Grid,points.S,data,1,0,NORM);

    fit.F=F;
    fit.bw=bw;
    fit.norm=norma;
    COPULA.data=data;
    COPULA.fit=fit;

else
    
    if fitt==-1
        lfit = NPC_loclik_fit(bw,data,Grid);
        fit.lfit=[];
    end
    [pd_data,ccdf_data,pd_grid,ccdf_grid,pd_points,ccdf_points,Fn,~]=NPC_func_tll(Lfit,Grid,points.S,data,fitt,0,NORM);
    
    if fitt==0
        COPULA.fit=0;
        COPULA.Grid=0;
        COPULA.data=0;
        COPULA.points=0;
    elseif fitt==-1
        COPULA=Lfit;
        fit.F=[];
        fit.F=Fn;
        fit.bw=Lfit.bw;
        COPULA.data=data;
        COPULA.fit=fit;
    end
    
end

return
