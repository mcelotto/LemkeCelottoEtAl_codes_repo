function u=NPC_kerncopccdfinv(copula,v,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(copula,'knots')
knots = copula.knots;
else
knots = copula.fit.knots;
end
[~,GRID_u]= NPC_mk_grid(knots,method);
copula.Grid.u=double(GRID_u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cc=reshape(copula.C_grid,knots,knots);

S=copula.Grid.u(:,1);
SS=reshape(S,sqrt(size(S,1)),sqrt(size(S,1)));
X=SS(:,1);
S=copula.Grid.u(:,2);
SS=reshape(S,sqrt(size(S,1)),sqrt(size(S,1)));
Y=SS(1,:);

G=copula.C_grid;

ww=[];
for i=1:size(v,1)

[~,m1]=min(abs(X-v(i,1))); 
g=G(m1,:);

gg=find(g>=v(i,2));
ggg=g(gg);
[~,m22]=min(abs(v(i,2)-ggg));
m2=gg(m22);
u(i)=Y(m2);
% u(i)=interp1(g,Y,v(i,2),'linear','extrap');
end

