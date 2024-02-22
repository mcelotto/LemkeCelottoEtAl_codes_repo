% #' Evaluate the density of the transformation log likelihood estimator


function [pd_data,ccdf_data,pd_grid,ccdf_grid,pd_points,ccdf_points,F,norma] = NPC_func_tll(lfit,Grid,points,data,fit,short,NORM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Grid.X,2)==1
    kn=round((size(Grid.u,1)));
end
if size(Grid.X,2)==2
    kn=round(sqrt(size(Grid.u,1)));
end
if size(Grid.X,2)==3
    kn=round((size(Grid.u,1)^(1/3)));
end

if size(Grid.X,2)==1
    x1=unique(Grid.u(:,1));
    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
end
if size(Grid.X,2)==2
    x1=unique(Grid.u(:,1));
    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
    x2=unique(Grid.u(:,2));
    xd2=[diff(x2)];xd2=[xd2;xd2(1)];
end
if size(Grid.X,2)==3
    x1=unique(Grid.u(:,1));
    xd1=[diff(x1)];xd1=[xd1;xd1(1)];
    x2=unique(Grid.u(:,2));
    xd2=[diff(x2)];xd2=[xd2;xd2(1)];
    x3=unique(Grid.u(:,3));
    xd3=[diff(x3)];xd3=[xd3;xd3(1)];
end

if isfield(lfit,'Kergrid')
    
    if size(Grid.X,2)==1
        t1=(lfit.Kergrid);
        F1=griddedInterpolant(Grid.S,double(t1),'linear','none');
    end
    
    if size(Grid.X,2)==2
        t1=reshape(lfit.Kergrid,kn,kn);
        F1=scatteredInterpolant(Grid.S,double(t1(:)),'linear','none');
    end   
else
    LF = NPC_loclik_fit(lfit.bw,data,Grid);
    if size(Grid.X,2)==1
        t1=reshape(LF.Kergrid,kn);
    end
    if size(Grid.X,2)==2
        t1=reshape(LF.Kergrid,kn,kn);
    end
    F1=scatteredInterpolant(Grid.S,double(t1(:)),'linear','none');
end


if size(Grid.X,2)==1
    if short==0
        pd_grid=t1(:);
        tnorm=pd_grid;
        
        for n=1
            I2=sum(xd1.*tnorm,1);
            K=I2;
            tnorm=tnorm./K;
        end
        pd_grid=tnorm;
    end
end


%%%%% Normalisation
%%%%% f(S,T)/phi(S)phi(T)--> c(u,v)  


if short==1
    NUM=50;
else
    NUM=1000;
end

if size(Grid.X,2)==2
    if  short~=1
        pd_grid=t1(:);
        t1=reshape(pd_grid,kn,kn);
        tnorm=t1./NORM;
        for n=1:NUM
            I1=sum(xd2'.*tnorm,2);
            I2=sum(xd1.*tnorm,1);
            
            K=I1.*I2;
            tnorm=tnorm./K;
            
        end
            II=sum(xd1.*sum(xd2'.*tnorm,2));
            tnorm=tnorm/II;           
        pd_grid=tnorm;
    end
end


if short==2 & fit==0
        t1=pd_grid.*NORM;
        F1=scatteredInterpolant(Grid.S,double(t1(:)),'linear','none');
        short=1;
end

norma=1;

if short==0 | fit==-1
    
    if fit~=0
        if size(Grid.X,2)==1
            F.pdf=griddedInterpolant(Grid.S,double(pd_grid(:)),'linear','none');
        end
        if size(Grid.X,2)==2
            F.pdf=scatteredInterpolant(Grid.S,double(pd_grid(:)),'nearest','none');
        end        
    else
        
        if size(Grid.X,2)==2 || size(Grid.X,2)==1
            F.pdf=lfit.F.pdf;
        end
        
    end
    
    if size(Grid.X,2)==2 || size(Grid.X,2)==1
        pd_points=F.pdf(points);
        pd_data=F.pdf(data.S);
        pd_data(pd_data<0)=0;
        pd_points(pd_points<0)=0;
        pd_points(isnan(pd_points))=0;
        pd_data(isnan(pd_data))=0;
    end
    
else
    
    if size(Grid.X,2)==2 | size(Grid.X,2)==1
        pd_grid=t1;
        pd_points=NaN;
        pd_data=F1(data.S);
        pd_data(pd_data<0)=eps;
        pd_data(isnan(pd_data))=eps;
    end
        
end

if short==0
    t1=pd_grid;    

    if size(Grid.X,2)==1
        h=t1*NaN;        
    end

    if size(Grid.X,2)==2
        h=t1*NaN;
        for i2=1:size(t1,2)
            for i1=1:size(t1,1)
                if i2==1
                    h(i1,1)=t1(i1,1)*xd2(1)/sum(t1(i1,1:end)'.*xd2(1:end));                    
                else
                    if sum(t1(i1,1:end)'.*xd2(1:end))~=0
                        h(i1,i2)=sum(t1(i1,1:i2)'.*xd2(1:i2))/sum(t1(i1,1:end)'.*xd2(1:end));
                    else
                        h(i1,i2)=0;
                    end
                end
            end
        end        
    end
    
    [~,GRID_u]= NPC_mk_grid(kn,'');

        h(h<=min(Grid.u(:)))=min(GRID_u(:));
        h(h>=max(Grid.u(:)))=max(GRID_u(:));    
    
    if fit~=0
        if size(Grid.X,2)==2
            F.ccdf=scatteredInterpolant(Grid.S,double(h(:)),'linear','none');
        end
    else
        if size(Grid.X,2)==2
            F.ccdf=lfit.F.ccdf;
        end
                
    end
    
    
    if size(Grid.X,2)==1
        ccdf_points=NaN;
        ccdf_data=NaN;
        ccdf_grid=NaN;        
    end
    
    if size(Grid.X,2)==2
        ccdf_points=F.ccdf(points);
        ccdf_data=F.ccdf(data.S);
        ccdf_grid=h;
        ccdf_points(isnan(ccdf_points))=0;
        ccdf_data(isnan(ccdf_data))=0;
        ccdf_data(ccdf_data<=min(GRID_u(:)))=min(GRID_u(:));
        ccdf_data(ccdf_data>=max(GRID_u(:)))=max(GRID_u(:));
        ccdf_points(ccdf_points<=min(GRID_u(:)))=min(GRID_u(:));
        ccdf_points(ccdf_points>=max(GRID_u(:)))=max(GRID_u(:));
        
    end
    
else
    ccdf_grid=NaN;
    ccdf_data=NaN;
    ccdf_points=NaN;
    F=NaN;
end





