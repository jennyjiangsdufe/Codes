function [AR,RI,MI,HI]=RandIndex(c1,c2)


if nargin < 2 || min(size(c1)) > 1 || min(size(c2)) > 1
   error('RandIndex: Requires two vector arguments')
   return
end

C=Contingency(c1,c2);
n=sum(sum(C));
nis=sum(sum(C,2).^2);		
njs=sum(sum(C,1).^2);		

t1=nchoosek(n,2);		
t2=sum(sum(C.^2));	
t3=.5*(nis+njs);

nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;	
D=  -t2+t3;		

if t1==nc
   AR=0;			
else
   AR=(A-nc)/(t1-nc);		
end

RI=A/t1;			
MI=D/t1;			
HI=(A-D)/t1;	

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 || min(size(Mem1)) > 1 || min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1)
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end
