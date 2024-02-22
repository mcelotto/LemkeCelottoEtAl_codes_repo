function [ker_grid]=NPC_Kern_LL(KerN_grid,B)


switch numel(B)
    case 1
        
        e1=B(1)*sqrt(abs(KerN_grid{3}./KerN_grid{1}-(KerN_grid{2}./KerN_grid{1}).^2));
        
        e1(isnan(e1))=0;
        
        C=-e1.^2.*(KerN_grid{2}./KerN_grid{1}).^2/(2*B(1)^2);
        
        C(isnan(C))=0;
        
        ker_grid=KerN_grid{1}.*e1.*exp(C);
        
    case 2
        
        e1=B(1)*sqrt(abs(KerN_grid{4}./KerN_grid{1}-(KerN_grid{2}./KerN_grid{1}).^2));
        e2=B(2)*sqrt(abs(KerN_grid{5}./KerN_grid{1}-(KerN_grid{3}./KerN_grid{1}).^2));
        
        e1(isnan(e1))=0;
        e2(isnan(e2))=0;
        
        C=-e1.^2.*(KerN_grid{2}./KerN_grid{1}).^2/(2*B(1)^2)-e2.^2.*(KerN_grid{3}./KerN_grid{1}).^2/(2*B(2)^2);
        
        C(isnan(C))=0;
        
        ker_grid=KerN_grid{1}.*e1.*e2.*exp(C);
        
    case 3
        
        e1=B(1)*sqrt(abs(KerN_grid{5}./KerN_grid{1}-(KerN_grid{2}./KerN_grid{1}).^2));
        e2=B(2)*sqrt(abs(KerN_grid{6}./KerN_grid{1}-(KerN_grid{3}./KerN_grid{1}).^2));
        e3=B(3)*sqrt(abs(KerN_grid{7}./KerN_grid{1}-(KerN_grid{4}./KerN_grid{1}).^2));
        
        e1(isnan(e1))=0;
        e2(isnan(e2))=0;
        e3(isnan(e3))=0;
        
        C=-e1.^2.*(KerN_grid{2}./KerN_grid{1}).^2/(2*B(1)^2)-e2.^2.*(KerN_grid{3}./KerN_grid{1}).^2/(2*B(2)^2)-e3.^2.*(KerN_grid{4}./KerN_grid{1}).^2/(2*B(3)^2);
        
        C(isnan(C))=0;
        
        ker_grid=KerN_grid{1}.*e1.*e2.*e3.*exp(C);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

