
function LL=LocLik(PAR,B,data,grid)

Ker0=zeros(size(grid,2),size(data,2));
KerI0=zeros(size(grid,2),size(data,2));
P0=zeros(size(grid,2),size(data,2));

Ker=zeros(size(grid,2),size(grid,2));
KerI=zeros(size(grid,2),size(grid,2));
P=zeros(size(grid,2),size(grid,2));

for i=1:size(data,2)
   d=data(:,i);
   dis=grid-d;

   P0(:,i)=PAR(1)+PAR(2)*dis(1,:)+PAR(3)*dis(2,:);
   
   Ker0(:,i)=exp(-dis(1,:).* dis(1,:))/(2*pi).*exp(-dis(2,:).* dis(2,:))/(2*pi);   
   KerI0(:,i)=Ker0(:,i).*(P0(:,i));   


end


for i=1:size(grid,2)
   d=grid(:,i);
   dis=grid-d;

   P(:,i)=PAR(1)+PAR(2)*dis(1,:)+PAR(3)*dis(2,:);
   
   Ker(:,i)=exp(-dis(1,:).* dis(1,:))/(2*pi).*exp(-dis(2,:).* dis(2,:))/(2*pi);   
   KerI(:,i)=Ker(:,i).*exp(P(:,i));   


end



LL=-sum(KerI0(:))+size(data,2)*sum(KerI(:));


