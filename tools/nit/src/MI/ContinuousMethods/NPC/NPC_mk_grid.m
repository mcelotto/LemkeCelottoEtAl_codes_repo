
function [pnts,expanded]=NPC_mk_grid(knots,method)

if nargin<2
    method=[1 1];
end
if isempty(method)
    method=[1 1];
end

if ~ischar(method)
    pnts=zeros(knots,numel(method));
    if numel(method)>1
        for m=1:numel(method)
            if method(m)==1
                pnts(:,m) = normcdf(linspace(-3.2,3.2,knots));
            else
                [~,GG]= NPC_mk_grid(knots,1);
                pnts(:,m) = linspace(min(GG(:)),max(GG(:)),knots);
            end
        end
        
        [~,ex2] = meshgrid(pnts(:,1));
        [ex1,~] = meshgrid(pnts(:,2));
        expanded=[ex2(:) ex1(:)];
    else
        if method==1
            pnts = normcdf(linspace(-3.2,3.2,knots));
        else
            [~,GG]= NPC_mk_grid(knots,'');
            pnts = linspace(min(GG(:)),max(GG(:)),knots);
        end
        expanded=pnts(:);
    end
else
    pnts = normcdf(linspace(-3.2,3.2,knots));
    [ex1,ex2] = meshgrid(pnts);
    expanded=[ex2(:) ex1(:)];
end


