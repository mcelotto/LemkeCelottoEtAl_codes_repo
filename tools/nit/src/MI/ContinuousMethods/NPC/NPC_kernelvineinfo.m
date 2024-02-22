function [info,stderr_tot,in,En] = NPC_kernelvineinfo(vines,copula,opts,pcond)
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
% code to compute mutual information given two vine copulas or between two vine copulas
% copyright Houman Safaai March 2018
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%

if isfield(opts,'parallel')
    parallel=opts.parallel;
else
    parallel=0;
end    
if isfield(opts,'alpha')
    alpha=opts.alpha;
else
    alpha=0.05;
end    
if isfield(opts,'erreps')
    erreps=opts.erreps;
else
    erreps=5e-3;
end
if isfield(opts,'cases')
    cases=opts.cases;
else
    cases=1e4;
end
if isfield(opts,'iter')
    kmin=opts.iter;
else
    kmin=50;
end
if isfield(opts,'plot')
    plo=opts.plot;
else
    plo=0;
end

if nargin>3
    pcond0=pcond;
    % Number of conditions
    ncond = length(cell2mat(pcond));
    % Probabilities of conditions
else
    ncond=1;
end

% Gaussian confidence interval for erreps and level alpha
conf = norminv(1 - alpha,0,1);
h=0;
hcond=0;
k = 0;

for w=1:15
    En(w,1) = 0;
    varsum(w,1) = 0;
    stder(w,1) = 0;
end

stderr1 = inf;
stderr2 = inf;
stderr_tot = inf;
varsum1 = 0;
varsum_tot = 0;
infoc1 = 0;


while (stderr1 >= erreps | stderr2 >= erreps | stderr_tot >= erreps) & k<kmin
    
        if length(vines)==1           
            
            infok = zeros(cases,1);
            infokc1 = zeros(cases,1);
            infokc = zeros(cases,1);
            
            k = k + 1;
            % Generate samples
            
            clear X pp qq pp_par qq_par
            for hh=1:numel(vines.margins)
                X(:,hh)=vines.margins{hh}.ker;
            end
                        
            D=NPC_kerncoprnd(copula,X,cases,vines);            
            
            [~,~,~,~,p_cop] = NPC_Fit_vCopula(vines,D,[],-1,copula,opts.knots_fit,parallel);

            
            log2pp=log2(p_cop);
            log2pp(p_cop==0)=0;
            I= log2pp  ;
            infokc1 =  I;
            
            infoc1 = infoc1 + ( nanmean(infokc1(~isinf(infokc1))) - infoc1) / k;
            
            varsum1 = varsum1 + nansum((infokc1(~isinf(infokc1)) - infoc1) .^ 2);
            stderr1 = conf * sqrt(varsum1 / (k * cases * (k * cases - 1)));
            info = (infoc1);
            
            in.info(k)=info;
            in.stderr(k)=stderr1;
            
            if plo
            figure(100)
            hold on
            subplot(1,2,1);plot(k,info,'O')
            title('Info')
            hold on
            subplot(1,2,2);plot(k,stderr1,'*');hold on;plot(k,stderr_tot,'O')
            title('std error')
            hold on
            drawnow
            end
            
        else
            
            
            clear infok infokc1 infokc2 infokc infokt pp qq pp_par qq_par D cc
            
            SET_CON=cell2mat(SET_CO);
            pcon=cell2mat(pcond);
            
            SETS=[];
            for i=1:numel(SET_CO)
            SETS=[SETS i*ones(1,numel(SET_CO{i}))];
            end
            conditions=unique(SETS);
            
            
            
            NSAM=200;
            k = k + 1;
            % Generate samples
            
            clear Dt
            
            if k==1
                for i=1:ncond
                    D{i}=[];
                    clear DA
                    for j=1:numel(vines{SET_CON(i)}.margins)
                        DA(:,j)=vines{SET_CON(i)}.margins{j}.ker;
                    end
                    [~,~,copS{SET_CON(i)},~,~] = Fit_vCopula(vines{SET_CON(i)},DA(1,:),size(copula{1},2),vines{SET_CON(i)}.METH,-1,copula{SET_CON(i)},'rand',[],150); %200
                    %                     [~,~,copSB{SET_CON(i)},~,~] = Fit_vCopula(vines{SET_CON(i)}.B,DA(1,2:end),size(copula{1},2),vines{SET_CON(i)}.B.METH,-1,vines{SET_CON(i)}.copulaB,'rand',[],100); %200
                    
                    %                     copula{SET_CON(i)}=copS{SET_CON(i)};
                    %                     vines{SET_CON(i)}.copulaB=copSB{SET_CON(i)};
                end
            end
                     
            
            parfor i=1:ncond
                [D{i} TV{i}]=vinecopula_sample(vines{SET_CON(i)},copS{SET_CON(i)},tvec,NSAM,cases,nn);    %%%% nn is time
            end
            
            
            for i=1:ncond
                if i==1
                    tind=TV{1};
                else
                    tind= intersect(tind,TV{i});
                end
            end

                        
            tvec=tvec(tind);
            
           
            for i=1:ncond
                if i==1
                    Dt=D{1};
                else
                    Dt=cat(1,Dt,D{i});
                end
            end
            
            if size(D{1},1)>10000
                for i=1:ncond
                    if i==1
                        Dt=D{1}(randsample(1:size(D{1},1),10000),:);
                    else
                        Dt=cat(1,Dt,D{i}(randsample(1:size(D{i},1),10000),:));
                    end
                end
            end
                      
            
            clear D Q
            D=Dt;
            
            [~,h2]=histc(D(:,nn),tvec);
            
            if min(h2)==0
                h2=h2+1;
            end
            
            
            DS=D(:,2:end);
            
            QB=[];
            clear Q Qs
            for j=1:ncond
                
                if j==1
                    %                     [~,~,cop{SET_CON(j)},~,~] = Fit_vCopula(vines{SET_CON(j)},D(1:round(size(D,1)/2),:),size(copula{1},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                else
                    %                     [~,~,cop{SET_CON(i)},~,~] = Fit_vCopula(vines{SET_CON(j)},D(round(size(D,1)/2)+1:end,:),size(copula{1},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                end
                
                [Q(:,j),~,~,~,~] = Fit_vCopula(vines{SET_CON(j)},D,size(copula{SET_CON(1)},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                %                 [Qs(:,j),~,~,~,~] = Fit_vCopula(vines{SET_CON(j)},D,size(copula{1},2),vines{SET_CON(j)}.METH,-1,cop{SET_CON(j)},'rand',[]);
                
                if isfield(vines{SET_CON(j)},'B')
                    if j==1
                        %                     [~,~,copB{SET_CON(j)},~,~] = Fit_vCopula(vines{SET_CON(j)}.B,DS(1:round(size(D,1)/2),:),size(vines{SET_CON(1)}.copulaB,2),vines{SET_CON(j)}.B.METH,-1,vines{SET_CON(j)}.copulaB,'rand',[]);
                    else
                        %                     [~,~,copB{SET_CON(j)},~,~] = Fit_vCopula(vines{SET_CON(j)}.B,DS(round(size(D,1)/2)+1:end,:),size(vines{SET_CON(1)}.copulaB,2),vines{SET_CON(j)}.B.METH,-1,vines{SET_CON(j)}.copulaB,'rand',[]);
                    end
                    
                    %                 QB=Q;
                    
                    [QB(:,j),~,~,~,~] = Fit_vCopula(vines{SET_CON(j)}.B,DS,size(vines{SET_CON(j)}.copulaB,2),vines{SET_CON(j)}.B.METH,-1,vines{SET_CON(j)}.copulaB,'rand',[]);
                    %                     [QsB(:,j),~,~,~,~] = Fit_vCopula(vines{SET_CON(j)}.B,DS,size(vines{SET_CON(j)}.copulaB,2),vines{SET_CON(j)}.B.METH,-1,copB{SET_CON(j)},'rand',[]);
                end
            end
            
            
            %             copula=cop;
            %             Q=Qs;
            %             QB=QsB;
            
            for i=1:ncond
                if i==1
                    LAB=ones(round(size(D,1)/ncond),1);
                else
                    LAB=cat(1,LAB,SETS(i)*ones(round(size(D,1)/ncond),1));
                end
            end
            
                        
            NNN=250;
            dims=[1,2];
                        
            y_vectorT=linspace(vines{1}.range(mm,1),vines{1}.range(mm,2),500);
            
            for j=1:numel(SET_CON)
                if j==1
                    yy=linspace(min(vines{SET_CON(j)}.margins{mm}.ker),max(vines{SET_CON(j)}.margins{mm}.ker),250)';
                else
                    yy=cat(1,yy,linspace(min(vines{SET_CON(j)}.margins{mm}.ker),max(vines{SET_CON(j)}.margins{mm}.ker),250)');
                end
            end
            
            y_vectorT=sort(unique(yy));
            
            clear pMar pdf pdfTV pdfn pp qq qqc qqcc pq pmmm fm pmmmT fmT pppp
            
            for j = 1:ncond
                
                if size(vines{SET_CON(1)}.range,1)>2
                    [vineMar{j},copulaMar{j}]=kernelmarginalize(dims,vines{SET_CON(j)},copula{SET_CON(j)});
                    [pMar(:,j),~,~,~,~] = Fit_vCopula(vineMar{(j)},D(:,dims),size(copulaMar{1},2),vineMar{(j)}.METH,-1,copulaMar{(j)},'rand',[]); %% f(n,t|C)
                else
                    vineMar{j}=vines{SET_CON(j)};
                    copulaMar{j}=copula{SET_CON(j)};
                    pMar(:,j)=Q(:,j);
                end
                
                %%%%%% LATER I SHOULD ADD I(n;b_i| B_{-i} ) for
                %%%%%% each component of B
                %                             for x=1:size(D,2)-2
                %                             [pSingle{x}(:,j),~,~,~,~] = Fit_vCopula(vineMar{j},D(:,[1 x+2]),size(copulaMar{1},2),[],-1,copulaMar{j},'rand',[]); %% f(n,t|C)cc
                %                             end
                
                p=Q(:,j);
                
                par.fit=0;par.s=copula{SET_CON(j)}{nn,1}.MarginS.s;par.p=copula{SET_CON(j)}{nn,1}.MarginS.p;
                [pdf(:,j),~]=kernelpdf(D(:,nn),2,D(:,nn),par);
                
                par.fit=0;par.s=copula{SET_CON(j)}{mm,1}.MarginS.s;par.p=copula{SET_CON(j)}{mm,1}.MarginS.p;
                [pdfn(:,j),~]=kernelpdf(D(:,mm),2,D(:,mm),par);
                
                X=p./pdf(:,j);
                qq(:,j)=X;   %%% f(B,n|t,C)
                
                X=p./pMar(:,j); %%% f(B|t,n,C)
                qqc(:,j)=X;
                
                X=pMar(:,j)./pdf(:,j);
                pp(:,j)=X;   %%% f(n|t,C)
                
                X=pdf(:,j);
                pppp(:,j)=X;   %%% f(t|C)
                
                % Marginalizing the n
                if nn~=1 & isempty(QB)
                    
                    for sor=1:round(numel(y_vectorT)/NNN)
                        
                        %disp(['condition = ', num2str(j),' , sor = ',num2str(sor)])
                        
                        if sor*NNN<numel(y_vectorT)
                            y_vector=y_vectorT((sor-1)*NNN+1:sor*NNN);
                        else
                            y_vector=y_vectorT((sor-1)*NNN+1:end);
                        end
                        ne=0;
                        points=zeros(numel(y_vector)*size(D,1),size(D,2));
                        mmm=setdiff(1:size(D,2),mm);
                        
                        for i2=1:size(D,1)
                            for i1=1:numel(y_vector)
                                ne=ne+1;
                                Z(mm)=y_vector(i1);
                                Z(mmm)=D(i2,mmm);
                                points(ne,:)=Z;
                            end
                        end
                        
                        [pT,~,~,~,~] = Fit_vCopula(vines{SET_CON(j)},points,size(copula{1},2),vines{SET_CON(j)}.METH,-1,copula{SET_CON(j)},'rand',[]);
                        
                        ppp=reshape(pT,numel(y_vector),[]);
                        
                        if sor==1
                            pm=ppp;
                        else
                            pm=cat(1,pm,ppp);
                        end
                    end
                    
                    %                     pmm=sum(pm,1)'*mean(diff(y_vectorT));
                    pmm=sum(pm(1:end-1,:).*diff(y_vectorT),1)';
                    X=pmm./pdf(:,j);
                    qqcc(:,j)=X;
                    
                else
                    
                    if ~isempty(QB)
                        qqcc(:,j)=QB(:,j)./pdf(:,j);
                    else
                        qqcc=qq;
                    end
                    
                end
                
            end
            
            
            
            Vnbt=prod(diff(vines{1}.range(setdiff(1:size(Dt,2),nn),:)'));
            Vbt=prod(diff(vines{1}.range(setdiff(1:size(Dt,2),[nn mm]),:)'));
            Vn=prod(diff(vines{1}.range(mm,:)'));
            
            
            %                         qqcc(qqcc<1e-1/Vbt)=0;
            %                         qq(qq<1e-1/Vnbt)=0;
            %                         pp(pp<1e-1/Vn)=0;
            
%%%% convert to condtions, summing over sub_conditions
            qq0=qq;
            pp0=pp;
            qqcc0=qqcc;
            Q0=Q;
            clear qq pp qqcc Q
            
            PO=repmat(pcon,size(qq0,1),1);
            for i=1:numel(conditions)
                Q(:,i)=sum(PO(:,find(SETS==i)).*Q0(:,find(SETS==i)),2);
                qq(:,i)=sum(PO(:,find(SETS==i)).*qq0(:,find(SETS==i)),2);
                qqcc(:,i)=sum(PO(:,find(SETS==i)).*qqcc0(:,find(SETS==i)),2);
                pp(:,i)=sum(PO(:,find(SETS==i)).*pp0(:,find(SETS==i)),2);
            end
            
            for tt=1:numel(tvec)-1
                
                
                ff=find(h2==tt);
                
                qw=Q(ff,:);
                yy = quantile(qw(:),0.05);
                
                
                %                 ff=find(h2==tt & (Q(:,1)>yy | Q(:,2)>yy)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% this threshhold SHOULD BE CHECKED
                
                clear PR
                for i=conditions
                    PR(1,i)=sum(LAB(ff)==i)/numel(ff);
                end
                
%                                 PR=[sum(LAB(ff)==1) sum(LAB(ff)==2)]/numel(ff);%[0.5 0.5];%
                
                PR=repmat(PR,numel(ff),1);
                %                                             PR=(PR.*pdf(ff,:))./repmat(sum(PR.*pdf(ff,:),2),1,2);
                
                pdfnc=sum(PR .* qqcc(ff,:),2);
                pdfncMeasure0=sum(PR .* qq(ff,:),2);
                pdfncMeasure00=sum(PR .* pp(ff,:),2);
                
                
                for co=conditions
                    G=-(qqcc(ff,co).*log2(qqcc(ff,co))./pdfnc);
                    In= PR(:,co).*G;
                    
                    %                                                                 G=-log2(qqcc(ff(LAB(ff)==co),co));%G(qqcc(ff(LAB(ff)==co),co)==0)=0;
                    %                                                                 In= PR(LAB(ff)==co,co).*G;
                    
                    In(isinf(abs(In)) | isnan(In))=0;
                    H{co+2,tt}=In;
                    H1{co+2,tt}=In;
                end
                
                G=-log2(pdfnc);
                In= G;
                In(isinf(abs(In)) | isnan(In))=0;
                H{15,tt}=In;                   %%H(B|C)
                
                %%%%%%%%%%%
                
                G=-log2(pdfncMeasure0);
                I= G;
                I(isinf(abs(I)) | isnan(I))=0;
                H{5,tt}=I;
                
                G=-log2(pdfnc.*pdfncMeasure00);
                I= G;
                I(isinf(abs(I)) | isnan(I))=0;
                H{6,tt}=I;
                
                for co=conditions
                    G=log2(qq(ff(LAB(ff)==co),co)./(pp(ff(LAB(ff)==co),co).*qqcc(ff(LAB(ff)==co),co)));
                    In= PR(LAB(ff)==co,co).*G;
                    In(isinf(abs(In)) | isnan(In))=0;
                    H{12+co,tt}=In;
                    H2{12+co,tt}=In;
                end
                
                %%%%%%%%%%%
                
                
                for co=conditions
                    G=-pp(ff,co).*log2(pp(ff,co))./pdfncMeasure00;
                    I= PR(:,co).*G;
                    I(isinf(abs(I)) | isnan(I))=0;
                    H{6+co,tt}=I;
                    H3{6+co,tt}=I;
                end
                
                G=-log2(pdfncMeasure00);
                In= G;
                In(isinf(abs(In)) | isnan(In))=0;
                H{9,tt}=In;
                
                %%%%%%%%%%
                
                for co=conditions
                    
                    G=-(qq(ff,co).*log2(qq(ff,co))./pdfncMeasure0);
                    I= PR(:,co).*G;
                    
                    %                                                                 G=-log2(qq(ff(LAB(ff)==co),co));%G(qq(ff(LAB(ff)==co),co)<1e-6)=NaN;
                    %                                                                 I= PR(LAB(ff)==co,co).*G;
                    
                    I(isinf(abs(I)) | isnan(I))=0;
                    H{9+co,tt}=I;
                    H4{9+co,tt}=I;
                end
                
                G=-log2(pdfncMeasure0);
                G(isinf(abs(G)) | isnan(G))=0;
                H{12,tt}=G;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ff~=0
                    for w=1:15
                        En(w,tt) = En(w,tt)+  (nanmean(H{w,tt}(~isinf(H{w,tt}))) - En(w,tt)) /k;
                        varsum(w,tt) = varsum(w,tt) + nansum((H{w,tt}(~isinf(H{w,tt})) - En(w,tt)).^ 2);
                        NN=numel(H{w,tt});
                        stder(w,tt) = conf * sqrt(varsum(w,tt) / (NN*k * (NN*k - 1)));
                    end
                    
                    
                    for w=conditions
                        En1(w+2,tt) = En1(w+2,tt)+  (nanmean(H1{w+2,tt}(~isinf(H1{w+2,tt}))) - En1(w+2,tt)) /k;
                        En2(w+12,tt) = En2(w+12,tt)+  (nanmean(H2{w+12,tt}(~isinf(H2{w+12,tt}))) - En2(w+12,tt)) /k;
                        En3(w+6,tt) = En3(w+6,tt)+  (nanmean(H3{w+6,tt}(~isinf(H3{w+6,tt}))) - En3(w+6,tt)) /k;
                        En4(w+9,tt) = En4(w+9,tt)+  (nanmean(H4{w+9,tt}(~isinf(H4{w+9,tt}))) - En4(w+9,tt)) /k;
                    end
                end
                
            end
            
            
            clear in stde
            
            for tt=1:numel(tvec)-1
                
                %%%   I(n;C)
                                                
%                 in.info(1,tt)=En(9,tt)-En(8,tt)-En(7,tt);
                in.info(1,tt)=En(9,tt)-sum(En3(7:6+ncond,tt));
                
                %%%   I(n;B)
                
                in.info(2,tt)=En(6,tt)-En(5,tt);
                
                %%%   I(B;C)
                                
%                 in.info(3,tt)=En(15,tt)-En(3,tt)-En(4,tt);
                in.info(3,tt)=En(15,tt)-sum(En1(3:2+ncond,tt));
                
                %%% I(n,B;C)
                
%                 in.info(4,tt)=En(12,tt)-En(10,tt)-En(11,tt);
                in.info(4,tt)=En(12,tt)-sum(En4(10:9+ncond,tt));
                
                %%% I(n;C | B)
                
                %                             in.info(in.info<0)=0;
                
                in.info(5,tt)=in.info(4,tt) - in.info(3,tt);
                
                %%% I(n;B | C)
                                                
%                 in.info(6,tt)=En(13,tt)+En(14,tt);
                in.info(6,tt)=sum(En2(13:12+ncond,tt));
                
                %%% I(n;B , C)
                
                in.info(7,tt)=in.info(5,tt) + in.info(2,tt);
                
                
                stde(1,tt)=stder(3,tt);
                stde(2,tt)=stder(4,tt);
                stde(3,tt)=stder(15,tt);
                stde(4,tt)=stder(10,tt);
                stde(5,tt)=stder(11,tt);
                stde(6,tt)=stder(12,tt);
                
            end
            
            
            %             in.info(in.info<0)=0;
            
            
            stderr1=0;%norm(stderr1T(~isnan(stderr1T) & ~isinf(stderr1T)));
            stderr2=0;%norm(stderr2T(~isnan(stderr2T) & ~isinf(stderr2T)));
            stderr_tot=max(stder(:));
            
            in.stderr(k)=stderr_tot;
            
            
            
            figure(100)
            subplot(2,3,1)
            plot(tvec(1:end-1),in.info(1,:))
            hold on
            %                         xlim([5 290])
            title('I(n;C)')
            
            subplot(2,3,2)
            plot(tvec(1:end-1),in.info(4,:))
            hold on
            %                                                 plot(tvec(1:end-1),squeeze(mean(HH(:,12,:)-HH(:,10,:)-HH(:,11,:),1)),'.-')
            title('I(n,B;C)')
            
            subplot(2,3,3)
            plot(tvec(1:end-1),in.info(5,:))
            hold on
            %                                                 xlim([5 290])
            title('I(n;C|B)')
            
            subplot(2,3,4)
            plot(tvec(1:end-1),in.info(3,:))
            hold on
            %                         plot(tvec(1:end-1),inSH.info(3,:),'r')
            %                         xlim([5 290])
            title('I(B;C)')
            
            subplot(2,3,5)
            plot(tvec(1:end-1),in.info(2,:))
            hold on
            %                         xlim([5 290])
            title('I(n;B)')
            
            subplot(2,3,6)
            plot(tvec(1:end-1),in.info(6,:))
            hold on
            %                         plot(tvec(1:end-1),in.info(4,:)-in.info(1,:))
            %                         xlim([5 290])
            title('I(n;B|C)')
            drawnow
            
            
            
            
        end

        
    
end

