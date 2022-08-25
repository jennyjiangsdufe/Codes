function CE=Statistic(V,k,W,K)
for i=1:K
    cri=find(W==i);
    NK=zeros(1,K);
    for j=1:k
        a=find(V(cri)==j);
        NK(j)=length(a);
    end 
    t=find(NK==max(NK));
    NS(i)=t(1);
    NR(i)=max(NK); 
end
for i=1:k
    [A B]=find(NS==i);
    if length(B)>=2
        t=max(NR(B));
        NR(B)=0;
        NR(B(1))=t;
    end
end
CE=sum(NR)/length(V);