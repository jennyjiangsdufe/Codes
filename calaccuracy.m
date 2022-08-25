
function [precision,recall,F,FMeasure,Accuracy] = calaccuracy(P,C)
 % P为人工标记簇
 % C为聚类算法计算结果
 N = length(C);
 p = unique(P);
 c = unique(C);
 P_size = length(p);
 C_size = length(c);
 Pid = double(ones(P_size,1)*P' == p*ones(1,N) );
 Cid = double(ones(C_size,1)*C' == c*ones(1,N) );
 CP = Cid*Pid';
 Pj = sum(CP,1);
Ci = sum(CP,2);

precision = CP./( Ci*ones(1,P_size) );
recall = CP./( ones(C_size,1)*Pj );
F = 2*precision.*recall./(precision+recall);
FMeasure = sum( (Pj./sum(Pj)).*max(F) );
Accuracy = sum(max(CP,[],2))/N;

end
