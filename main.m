clear;
clc
%�������ݼ�
Dataset = importdata('Soybean_small.txt');
N = size(Dataset, 2)-1;  
data=Dataset(:,1:N);
V=Dataset(:,N+1);
K=numel(unique(V))
times=30;      
%% Ԥ�ȼ������sigma
data_n = size(Dataset, 1); 
N = size(Dataset, 2)-1;   
data=Dataset(:,1:N);
dist=zeros(data_n,data_n);
for i=1:data_n
    for j=1:data_n
        dist(i,j)=sum(data(i,:)~=data(j,:));
        dist(i,i)=1000;
    end
end
[min_dist,index]=min(dist,[],2);
sigma1=max(min_dist);

% ������
for t=1:1:times
    fprintf('ʵ���ظ����е�%d��\n',t);
 [center, U, obj_fcn,Accuracy, RI, NMI,FMeasure] = KIWFKMDP(Dataset, K, sigma1);
 result_RI(t)=RI;
 result_NMI(t)=NMI;
 result_Accuracy(t)=Accuracy;
 result_FMeasure(t)=FMeasure;
end
Avg_RI=sum(result_RI)/times
Avg_NMI=sum(result_NMI)/times
Avg_Accuracy=sum(result_Accuracy)/times
Avg_FMeasure=sum(result_FMeasure)/times



