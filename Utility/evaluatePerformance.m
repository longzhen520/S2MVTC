function [res_cluster,time] = evaluatePerformance(beta, gamma, lambda, alpha_1, alpha_2)
 
% addpath('./Caltech101-Datasets')
addpath('./Utility')
load Caltech101-all
%%  Nonlinear anchor feature embedding---------------------简单非线性RBF映射
V = size(X,2);
load Anchor_Caltech101
%fprintf('Nonlinear Anchor Embedding...\n');
for it = 1:V
%fprintf('The %d-th view Nonlinear Anchor Embeeding...\n',it);
    dist = EuDist2(X{it},Anchor{it},0); %计算X与Anchor之间欧式距离,bSqrt 参数被设置为 0，表示不进行平方根处理
    sigma = mean(min(dist,[],2).^0.5)*2;%计算sigma，它是距离矩阵中每行最小值的平方根的均值的两倍
    feaVec = exp(-dist/(2*sigma*sigma));%计算 dist 距离矩阵元素的指数函数值
    X{it} = bsxfun(@minus, feaVec', mean(feaVec',2));% Centered data 将 feaVec' 中的每一列减去对应列的平均值  @min减法  % 中心化数据
end
clear feaVec dist sigma dist Anchor it


%% ----------------------Initializing parameters------------------------    初始化参数
MaxIter = 10;       % 5 iterations are okay, but better results for 10
innerMax = 10;      % 内存最大值
K = 128;            % Hashing code length



N = size(X{1},2);%sample number

N = size(X{1},2);
%初始化B_bar
rand('seed',100);
sel_sample = X{4}(:,randsample(N, 1000),:);
[pcaW, ~] = eigs(cov(sel_sample'), K);  %sel_sample转置的协方差矩阵 前 L 个特征值和对应的特征向量并存储
B_bar = sign(pcaW'*X{4});                   %用二值化函数 sign() 对矩阵 pcaW' * X{4} 进行二值化处理

%初始化B
for v =1:V
    rand('seed',100);
    sel_sample = X{v}(:,randsample(N, 1000),:);
    [pcaW, ~] = eigs(cov(sel_sample'), K);  %sel_sample转置的协方差矩阵 前 K 个特征值和对应的特征向量并存储  ?????????XXT需不需要倒
    B{v}= sign(pcaW'*X{v});          %用二值化函数 sign() 对矩阵 pcaW' * X{v} 进行二值化处理
end

% 初始化Lambda_1矩阵
Lambda_1 = cell(1, V);
for v =1:V
    Lambda_1{v}= zeros(K, N);          %用二值化函数 sign() 对矩阵 pcaW' * X{v} 进行二值化处理
end
% 初始化p张量————b是张量（Bv更新以后在合并在一起）
B_tensor = cat(3, B{:,:});
P_tensor = B_tensor;
% Lambda_2是张量
Lambda_2_tensor = zeros(K, N, V);  
dim1 = K;dim2 = N;dim3 = V;
sX = [K, N, V];
        
C = numel(unique(Y));           %聚类的类别数目

%初始化D、G
rand('seed',500);
D = B_bar(:,randsample(N, C));
HamDist = 0.5*(K - B_bar'*D);
[~,ind] = min(HamDist,[],2);
G = sparse(ind,1:N,1,C,N,N);
G = full(G);
DG = C*G;

%%初始化输入以后
%  A  %  |矩阵 A\X\B          |张量             |标量   α_1     ————E F H矩阵
%  B  %  |矩阵 B_bar\Lambda_1 |张量  Lambda_2\P_tensor |标量   α_2\β\λ ————
%B_bar%  |矩阵                |张量  B__tensor  |标量   γ       ————
% DG  %  |矩阵 D\G            |张量             |标量   γ       ————




clear HamDist ind initInd n_randm pcaW sel_sample view
%% ----------------------End Initialization------------------------


tic;
%disp('----------The proposed method (multi-view)----------');
for iter = 1:MaxIter
    %fprintf('The %d-th iteration...\n',iter);
     %% Update A_v新---------------------------------------------------------  
     %svd_EFH 
        for v = 1:V
            M{v} = X{v} * Lambda_1{v}' + alpha_1 * X{v} * B{v}';
            [E{v}, F{v}, H{v}] = svd(M{v});
            A{v} = H{v} * E{v}(:,1:K)';
        end

    %% Update B_v新----------------------------------------------------------
        for v= 1:V
            numerator = beta * B_bar - Lambda_1{v} + alpha_1 * A{v} * X{v} + (Lambda_2_tensor(:, :, v) + alpha_2* P_tensor(:, :, v));%+ 这有个逆转换
            denominator = beta + alpha_1 + alpha_2 - 2 * lambda;
            B{v} = sign(numerator ./ denominator);
            B{v}(B{v}==0) = -1;
        end
           
    %% Update P改改改改改------------------------------------------------------------
        B_tensor = cat(3, B{:,:});
        b = B_tensor(:);
        l = Lambda_2_tensor(:);
        [p, objV] = wshrinkObj_weight(b + 1/alpha_2*l,alpha_2,sX,0,3);
        P_tensor = reshape(p, sX);

        history.objval(iter+1) = objV;


    %% Update B_bar------------------------------------------------------------
        % 初始化变量来存储 B 的和
        BX = zeros(size(B{1}));
        for v = 1:V
            BX = BX + B{v};
        end
        B_bar = (BX + gamma * D * G)  / (beta * V + gamma);



    %% Update D and G-------------------------------------------------------
       for iterInner = 1:innerMax
            % For simplicity, directly using DPLM here
            D = sign(B_bar*G'); 
            D(D==0) = 1;    %0值赋1 修正sign中的0
            theta = .001; tau = .01; % Preferred for this dataset
                for iterIn = 1:3
                    grad = -B_bar*G' + theta*repmat(sum(D),K,1);  %梯度
                    D = sign(D-1/tau*grad); D(D==0) = 1; 
                end
              HamDist = 0.5*(K - B_bar'*D); % 计算b二进制码 和c聚类簇心 之间的汉明距离Hamming distance referring to "Supervised Hashing with Kernels"
              [~,indx] = min(HamDist,[],2);%indx是i元素是i行最小值的列索引
              G = sparse(indx,1:N,1,C,N,N);
        end
        DG = D*G;

    % clear iterIn grad HamDist indx tau theta
    %% Update Lambda_1 Lambda_2------------------------------------------------------------
        for v= 1:V
         Lambda_1{v} = Lambda_1{v} + alpha_1 *(B{v} - A{v} * X{v});
        end     
        Lambda_2_tensor = Lambda_2_tensor + alpha_2 * (P_tensor - B_tensor);
end

elapsed_time = toc;
%disp('----------Main Iteration Completed----------');

[~,pred_label] = max(G,[],1);
res_cluster = ClusteringMeasure(Y, pred_label);
time_char = num2str(elapsed_time);
time = str2double(time_char);

end
