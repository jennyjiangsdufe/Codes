function [center, U, obj_fcn,Accuracy, RI, NMI,FMeasure] = KIWFKMDP(Dataset, K, sigma1, options)
tic
if nargin ~= 2 & nargin ~= 3,
	error('Too many or too few input arguments!');
end

default_options = [2;	
    10;              
    1e-5;             
    1                   
    2                 
    ];
if nargin == 3,
	options = default_options;
 else      
	if length(options) < 4,
		tmp = default_options;
		tmp(1:length(options)) = options;
		options = tmp;
    end
  	nan_index = find(isnan(options)==1);
  	options(nan_index) = default_options(nan_index);
	if options(1) <= 1, 
		error('The exponent should be greater than 1!');
	end
end
data_n = size(Dataset, 1);
N = size(Dataset, 2)-1;  
data=Dataset(:,1:N);
V=Dataset(:,N+1);
k=K;
cluster_n=K;
expo = options(1);          
max_iter = options(2);	
min_impro = options(3);		
display = options(4);	
expo2=options(5);         

pop_size=100;  
obj_fcn = zeros(pop_size, 1);	
fitness_optimum=zeros(max_iter, 1);
U_chrom=zeros(pop_size,data_n*K); 
CX=zeros(pop_size,data_n);
for n=1:1:pop_size              
U = initfcm(cluster_n, data_n);    
U_chrom(n,:)=reshape(U,1,[]);      
end
center = data(randperm(size(data,1),K),:);
GA_soap=U_chrom;       
NEW_GA_soap=zeros(pop_size,data_n*K+1);
Pm=0.1;
%% Main loop  主要循环
for i = 1:max_iter
for s=1:1:pop_size
U1=GA_soap(s,1:data_n*K);
U=reshape(U1,K,[]);
[cx,dist] = distfcm(center, U, data, K, N, expo, expo2, sigma1);
CX(s,:)=cx'; 
[U, center, obj_fcn(s)] = stepfcm(data, U, cluster_n, expo, expo2, sigma1, K, cx, N);
NEW_GA_soap(s,data_n*K+1)=obj_fcn(s);
Fitness=NEW_GA_soap(:,data_n*K+1); 
U2=reshape(U,1,[]);
NEW_GA_soap(s,1:data_n*K)=U2;
end

[order,Index]=sort(Fitness,'descend');
W=CX(Index(1),:)';
fitness_optimum(i)=Fitness(Index(1));
if i > 1,
		if abs(fitness_optimum(i) - fitness_optimum(i-1)) < min_impro, 
            break;
        end,
end
NEW_GA_soap=selection(NEW_GA_soap,Fitness,pop_size);
NEW_GA_soap=mutation(NEW_GA_soap,Pm,pop_size,K,data_n);
GA_soap=NEW_GA_soap; 
end
[Accuracy, RI, NMI,FMeasure]=performance_index(V, W);
fprintf('KIWFKM-DP:Iteration count = %d, RI=%f,NMI=%f,Accuracy=%f,FMeasure=%f\n', i,RI,NMI,Accuracy,FMeasure);  %第2种计算     
toc


%% 子函数1 
function U = initfcm(cluster_n, data_n)
U = rand(cluster_n, data_n);
col_sum = sum(U);
U = U./col_sum(ones(cluster_n, 1), :);

%% 子函数2 
function [U_new, center, obj_fcn] =stepfcm(data, U, cluster_n, expo, expo2, sigma1, K,cx, N)
U_a=U; %U_a表示隶属度
V_a=(1-(U_a).^expo).^(1/expo);  
Pi_a=1-U_a-V_a;   
U=Pi_a+U_a;

mf = U.^expo;      
center=compuecenter(data,cx,K);
[cx,dist] = distfcm(center, U, data, K, N, expo, expo2, sigma1);     
dist=dist';
obj_fcn = sum(sum((dist.^2).*mf)); 
tmp = dist.^(-2/(expo-1));     
U_new = tmp./(ones(cluster_n, 1)*sum(tmp)); 



%% 子函数3 
function [cx,dist] = distfcm(center, U, data, K, N, expo, expo2, sigma1)

n=size(data,1);
dist = zeros(n,size(center, 1));
cx=zeros(n,1);
mf = U.^expo;      
delta_1=zeros(n,N); 
delta=zeros(1,N);   
delta_w=zeros(1,N);
 for p=1:N
       for i = 1:n,   
        Kernel_dist=2*(1-exp(-(repmat(data(i,p),K,1)~=center(:,p))./sigma1));      
        delta_1(i,p)=U(:,i)'*(Kernel_dist);
        delta_1(isnan(delta_1)) = 0;
       end
 end
   delta=sum(delta_1,1);
   delta_w=1./((delta/sum(delta)).^(1/expo2-1));
    for i=1:n
    dist(i, :)=sum(((repmat(data(i,:),[K,1])~=center).*repmat(delta_w,K,1))',1);
    end
    [M,I]=min(dist,[],2);
    cx=I;  

%% 子函数4 
function center = compuecenter(data,cx,K)
    for i = 1:K
        center(i,:) = mode(data(cx==i,:));
    end
    
%% 子函数5 
function New_G=selection(G,F,G_Num)
for i=1:G_Num 
    r=rand*sum(F); 
    add_temp=0; 
    j=1; 
    while (add_temp<r)&(j<(size(G,1))+1) 
        add_temp=add_temp+F(j); 
        j=j+1; 
    end 
    if j==1 
        j=1; 
    else j=j-1; 
    end 
    New_G(i,:)=G(j,:); 
end 

%% 子函数6 
function G=mutation(G,Pm,pop_size,K,data_n) 
for i=1:1:pop_size         
    individual=G(i,1:K*data_n);
    individual=reshape(individual,K,[]);
    for j=1:1:size(individual,2) 
    Mr=rand;
    if Pm>Mr
      U_one = rand(K, 1); 
      col_sum = sum(U_one);
      U_one = U_one./col_sum;  
      individual(:,j)=U_one;
    end
    end
    individal=reshape(individual,1,[]);
    G(i,1:K*data_n)=individal;
end


