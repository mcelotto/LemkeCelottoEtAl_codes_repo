

function err=NPC_MISE(a,data,Grid,points,xd,afk,dd,ddT,bw0,nn,NORM,uspace,parallel)

if numel(a)==1 && nn==0
    bw=abs(a*bw0);
else
    if nn==0
        bw=abs(a);
    end
end

if uspace==0
    x1=unique(Grid.S(:,1));
    x2=unique(Grid.S(:,2));
    
    lfit = NPC_loclik_fit(bw,data,Grid);
    [~,~,pd_grid,~,~,~,~]=NPC_func_tll(lfit,Grid,points,data,0,2,NORM);
    
    if parallel
        parfor i=1:numel(afk.TestSize)
            if nn==0
                lfitCV = NPC_loclik_fit(bw(:),dd{i},Grid);
            else
                lfitCV = NPC_loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
            end
            [kkk{i},~,~,~,~,~,~]=NPC_func_tll(lfitCV,Grid,points,ddT{i},0,2,NORM);
        end
    else
        for i=1:numel(afk.TestSize)
            if nn==0
                lfitCV = NPC_loclik_fit(bw(:),dd{i},Grid);
            else
                lfitCV = NPC_loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
            end
            [kkk{i},~,~,~,~,~,~]=NPC_func_tll(lfitCV,Grid,points,ddT{i},0,2,NORM);
        end
    end
    
    
    for i=1:numel(afk.TestSize)
        pd_dataCV(afk.test(i)==1)=kkk{i}/sum(pd_grid(:)*xd);
    end
    
    e=pd_grid.^2;
    err=sum(e(:).*xd)-2*mean(pd_dataCV);
    
elseif uspace==1
    if parallel
        parfor i=1:numel(afk.TestSize)
            if nn==0
                lfitCV = NPC_loclik_fit(bw(:),dd{i},Grid);
            else
                lfitCV = NPC_loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
            end
            [kkk{i},~,~,~,~,~,~]=NPC_func_tll(lfitCV,Grid,points,ddT{i},0,1,NORM);
            kkk{i}=kkk{i}/(sum(lfitCV.Kergrid(:).*xd));
        end
    else
        for i=1:numel(afk.TestSize)
            if nn==0
                lfitCV = NPC_loclik_fit(bw(:),dd{i},Grid);
            else
                lfitCV = NPC_loclik_fit(bw(:,~afk.test(i)),dd{i},Grid);
            end
            [kkk{i},~,~,~,~,~,~]=NPC_func_tll(lfitCV,Grid,points,ddT{i},0,1,NORM);
            kkk{i}=kkk{i}/(sum(lfitCV.Kergrid(:).*xd));
        end
    end
    for i=1:numel(afk.TestSize)
        pd_dataCV(afk.test(i)==1)=kkk{i};
    end
    
    lfit = NPC_loclik_fit(bw,data,Grid);
    pd_grid=lfit.Kergrid/(sum(lfit.Kergrid(:).*xd));
    e=pd_grid.^2;
    err=sum(e(:).*xd)-2*mean(pd_dataCV);
    
end