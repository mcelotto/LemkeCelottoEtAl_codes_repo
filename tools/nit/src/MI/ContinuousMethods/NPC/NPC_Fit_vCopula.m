%%%%%% Fitting kernel vine copula and computing the joint density funciton
function [p,pdata,Copula,PDF,p_copula] = NPC_Fit_vCopula(vine,x,METH_fit,fit,lfit,knots,parallel)

if nargin < 6
    knots = lfit{1,2}.fit.knots;
end

d = length(vine.margins);

if ~isstruct(vine)
    error('NPC_Fit_vCopula: Argument "vine" must be a struct');
end
if ~ismatrix(x)
    error('NPC_Fit_vCopula: Argument "x" must be a matrix');
end
if size(x,2) ~= d
    error('NPC_Fit_vCopula: Second dimension of u must match number of margins');
end
for i=1:d
    PDF{i}=[];
end
cases = size(x,1);
newnode = true(d);
Fp = zeros(cases,d,d);
logf = zeros(cases,d);
logfdata=zeros(size(vine.theta{1,1},1),d);

no=0;
if fit==3
    fit=0;
    no=1;
end

[~,GRID_u]= NPC_mk_grid(knots,'');

for i = 1:d
    if fit~=0
        parG{i}.fit=fit;
        parS{i}.fit=fit;
        parT{i}.fit=fit;
        parP{i}.fit=fit;
        parG{i}.max=max(GRID_u(:));
        parG{i}.min=min(GRID_u(:));
        parT{i}.max=max(GRID_u(:));
        parT{i}.min=min(GRID_u(:));
    else
        parG{i}.fit=0;
        parG{i}.p=lfit{i,1}.MarginG.p;
        parG{i}.s=lfit{i,1}.MarginG.s;
        parS{i}.fit=0;
        parS{i}.p=lfit{i,1}.MarginS.p;
        parS{i}.s=lfit{i,1}.MarginS.s;
        parT{i}.fit=0;
        parT{i}.p=lfit{i,1}.MarginG.p;
        parT{i}.s=lfit{i,1}.MarginG.s;
        parP{i}.fit=0;
        parP{i}.p=lfit{i,1}.MarginS.p;
        parP{i}.s=lfit{i,1}.MarginS.s;
        
        parG{i}.max=max(GRID_u(:));
        parG{i}.min=min(GRID_u(:));
        parT{i}.max=max(GRID_u(:));
        parT{i}.min=min(GRID_u(:));
    end
end

for i = 1:d  
%     disp(['pdf ',num2str(i)])
    [G0{i},Mar_G{i}]=NPC_kernelcdf(vine.margins{i}.ker,x(:,i),parG{i});
    [S0{i},Mar_S{i}]=NPC_kernelpdf(vine.margins{i}.ker,vine.margins{i}.iscont,x(:,i),parS{i});
    [T0{i},Mar_T{i}]=NPC_kernelcdf(vine.margins{i}.ker,vine.margins{i}.ker,parT{i});
    [P0{i},Mar_P{i}]=NPC_kernelpdf(vine.margins{i}.ker,vine.margins{i}.iscont,vine.margins{i}.ker,parP{i});
    if fit==0 | fit==-1
        PDF{i}.margins=S0{i};
    end
end

for i = 1:d
    Fp(:,1,i)=G0{i};
    logf(:,1)=logf(:,1)+log(S0{i});
    logfdata(:,1)=logfdata(:,1)+log(P0{i});
    vine.theta{1,i}=T0{i};
end

clear G0 T0 S0

DD=1;

for tr=1
    
    clear co logftr G1p G2p G3p G4p CopC Copf Copul CCC
    
    if fit==1
        
        for str=2:d- (tr-1)   
            tic
            [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},Copul{str}]=NPC_KernelCop(knots,[vine.theta{tr,1} vine.theta{tr,str}],[vine.theta{tr,1}(1) vine.theta{tr,str}(1)],vine.METH{tr,str},METH_fit,0,1,parallel);
            logftr{str}=log(G2p{str});
            logfda{str}=log(G4p{str});
            %disp(['Fit copula -> ',num2str(tr),',',num2str(str),' , Time= ',num2str(toc)])
        end
        
    else
        
        for str=2:d- (tr-1)  
            Lfit=lfit{tr,str}.fit;
            [G1p{str},G2p{str},CopC{str},Copf{str},G3p{str},G4p{str},LF{str}]=NPC_KernelCop(knots,[vine.theta{tr,1} vine.theta{tr,str}],[Fp(:,tr,1) Fp(:,tr,str)],vine.METH{tr,str},METH_fit,Lfit,fit,parallel);
            logftr{str}=log(G2p{str});
            logfda{str}=log(G4p{str});
            if no~=1
                %disp(['Eval copula -> ',num2str(tr),',',num2str(str)])
            end
        end
    end
    
    for str=2:d- (tr-1)
        if fit~=0
            par1{str}.fit=1;
            par2{str}.fit=1;
            par1{str}.max=max(GRID_u(:));
            par1{str}.min=min(GRID_u(:));
            par2{str}.max=max(GRID_u(:));
            par2{str}.min=min(GRID_u(:));
        else
            par1{str}.fit=0;
            par1{str}.p=lfit{tr,str}.Margin1.p;
            par1{str}.s=lfit{tr,str}.Margin1.s;
            par2{str}.fit=0;
            par2{str}.p=lfit{tr,str}.Margin2.p;
            par2{str}.s=lfit{tr,str}.Margin2.s;
            par1{str}.max=max(GRID_u(:));
            par1{str}.min=min(GRID_u(:));
            par2{str}.max=max(GRID_u(:));
            par2{str}.min=min(GRID_u(:));
        end
    end
        
    clear Y1 Y2 Y3 Y4
    for str=2:d- (tr-1)
        Y1(:,str)=NaN*zeros(size(G1p{str}));
        Y2(:,str)=NaN*zeros(size(G3p{str}));
        Y3(:,str)=NaN*zeros(size(G4p{str}));
        [Y1(:,str),Mar_1{str}] = NPC_kernelcdf(G3p{str},G1p{str},par1{str});
        [Y2(:,str),Mar_2{str}] = NPC_kernelcdf(G3p{str},G3p{str},par2{str});
        Y3(:,str) = G4p{str};
    end
    
    ord=1:size(Y2,2)-1;
    
    for str=2:d- (tr-1)
        
        if fit==1
            Copula{tr,str}=Copul{str};
            Copula{tr,str}.fit.F='not saved here';
            Copula{tr,str}.fit.knots=knots;
            Copula{tr,str}.C_grid='not saved here';
            Copula{tr,str}.f_grid='not saved here';
            Copula{tr,str}.ord=ord;
        elseif fit==0
            Copula=NaN;
        elseif fit==-1
            Copula{tr,str}=LF{str};
            Copula{tr,str}.fit=lfit{tr,str}.fit;
            Copula{tr,str}.fit.F='not saved here';
            Copula{tr,str}.C_grid=CopC{str};
            Copula{tr,str}.f_grid=Copf{str};
            Copula{tr,str}.ord=ord;
            Copula{tr,str}.Margin1=[];
            Copula{tr,str}.Margin2=[];
            Copula{tr,str}.knots=size(Copula{tr,str}.C_grid,1);
        end
        Fp(:,tr+1,str-1) = Y1(:,str);
        vine.theta{tr+1,str-1} = Y2(:,str);
        vine.c{tr+1,str-1} = Y3(:,str);
        logf(:,tr+1) = logf(:,tr+1) + logftr{str};
        logfdata(:,tr+1) = logfdata(:,tr+1) + logfda{str};
    end
    
end

logp = logf(:,1);
logp_copula = 0;
for i = 2:d
    logp = logp + logf(:,i);
    logp_copula = logp_copula + logf(:,i);
end

logdata = logfdata(:,1);
for i = 2:d
    logdata = logdata + logfdata(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if fit==-1
    for i=1:d
        Copula{i,1}.MarginG=Mar_G{i};
        Copula{i,1}.MarginS=Mar_S{i};
    end
end

if fit==1
    [nor]=1;
    Copula{1,1}.norm=nor;
elseif fit==0
    if no==1
        nor=1;
    else
        nor=lfit{1,1}.norm;
    end
elseif fit==-1
    Copula{1,1}.norm=1;
    nor=Copula{1,1}.norm;
end

logp = reshape(logp,[cases 1]);
% Correct numerical inaccuracies
logp = real(logp);
% Massless intervals can result in NaNs; set to 0
logp(isnan(logp)) = -inf;
p = exp(logp)/nor;

logp_copula = reshape(logp_copula,[cases 1]);
% Correct numerical inaccuracies
logp_copula = real(logp_copula);
% Massless intervals can result in NaNs; set to 0
logp_copula(isnan(logp_copula)) = -inf;
p_copula = exp(logp_copula)/nor;

logdata = reshape(logdata,[size(vine.theta{tr,1},1) 1]);
% Correct numerical inaccuracies
logdata = real(logdata);
% Massless intervals can result in NaNs; set to 0
logdata(isnan(logdata)) = -inf;
pdata = exp(logdata)/nor;


end





