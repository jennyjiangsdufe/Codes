function [Entropy Purity]=EnAndPur(PC,Ci)
[rn cn]=size(PC);
%% º∆À„Ïÿ
for i=1:rn
    for j=1:cn
     precision(i,j)=PC(i,j)/Ci(1,i);    
    end
end

for i=1:rn
    for j=1:cn
     ei(i,j)=precision(i,j)*log2(precision(i,j));    
    end
end

for i=1:rn
    ei_sum(i)=-nansum(ei(i,:));
end

for j=1:cn
    mmi(j)=Ci(1,j)*ei_sum(j);
end

Entropy=nansum(mmi)/nansum(Ci);

for i=1:rn
     pr_max(i)=max(precision(i,:));    
end

for j=1:cn
    nni(j)=Ci(1,j)*pr_max(j);
end
Purity=nansum(nni)/nansum(Ci);
end

