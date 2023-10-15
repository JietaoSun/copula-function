

%%t3的权重版本！

% 1、提取各个监测点的权重残差数据
Monitor=readmatrix('监测点残差数据-正式3.xlsx');%从文件监测点残差数据.xlsx读取数据
X1_weighted=Monitor(:,1);%提取矩阵Monitor的第1列残差数据
X2_weighted=Monitor(:,2);%提取矩阵Monitor的第2列残差数据
X3_weighted=Monitor(:,3);%提取矩阵Monitor的第3列残差数据
X4_weighted=Monitor(:,4);%提取矩阵Monitor的第4列残差数据


% 权重核密度估计函数
[f_kx1_weighted,~,u1_weighted]=ksdensity(X1_weighted);%x1c是在第一列残差数据X1的取值范围内等间隔选取的的100个点构成的向量,f_kx1_weighted是与x1c_weighted相应的核密度估计值向量,ui是Gaussian核函数由样本数量、标准差计算的的最佳窗宽
[f_kx2_weighted,~,u2_weighted]=ksdensity(X2_weighted);%核密度估计值向量f_kx2_weighted
[f_kx3_weighted,~,u3_weighted]=ksdensity(X3_weighted);%核密度估计值向量f_kx3_weighted
[f_kx4_weighted,~,u4_weighted]=ksdensity(X4_weighted);%核密度估计值向量f_kx4_weighted

% 权重核分布估计函数
[F_kx1_weighted,x1c_weighted]=ksdensity(X1_weighted,'function','cdf');%x1c_weighted是在第一列残差数据X1的取值范围内等间隔选取的的100个点构成的向量,F_kx1_weighted是与x1c_weighted相应的核分布估计值向量
[F_kx2_weighted,x2c_weighted]=ksdensity(X2_weighted,'function','cdf');%核分布估计值向量F_kx2_weighted
[F_kx3_weighted,x3c_weighted]=ksdensity(X3_weighted,'function','cdf');%核分布估计值向量F_kx3_weighted
[F_kx4_weighted,x4c_weighted]=ksdensity(X4_weighted,'function','cdf');%核分布估计值向量F_kx4_weighted

%求残差数据Xi在权重核密度函数上取任一点xiq处的核密度函数值（缩略版的）
f_kx1q_weighted=interp1(x1c_weighted,f_kx1_weighted,X1);%interp1函数根据已有f_kx1_weighted核密度函数、已知数据点的横坐标x1c_weighted，求出横坐标X1上插值得到对应的函数值 f_kx1q_weighted；%示例vq = interp1(x,v,xq) 使用线性插值返回一维函数在特定查询点的插入值。向量 x 包含样本点，v 包含对应值 v(x)。向量 xq 包含查询点的坐标。
f_kx2q_weighted=interp1(x2c_weighted,f_kx2_weighted,X2);
f_kx3q_weighted=interp1(x3c_weighted,f_kx3_weighted,X3);
f_kx4q_weighted=interp1(x4c_weighted,f_kx4_weighted,X4);

%求残差数据Xi在权重核分布函数上取任一点xiq处的核分布函数值

F_kx1q_weighted=interp1(x1c_weighted,F_kx1_weighted,X1,'previous');%interp1函数通过线性或其他插值方式插值求出变量X1取任一点x1q处的核分布函数值，其中的'previous'参数表示取离x1q最近的一个数。
F_kx2q_weighted=interp1(x2c_weighted,F_kx2_weighted,X2,'previous');
F_kx3q_weighted=interp1(x3c_weighted,F_kx3_weighted,X3,'previous');
F_kx4q_weighted=interp1(x4c_weighted,F_kx4_weighted,X4,'previous');


%%
% t1产生的各个边缘分布函数再一次陈述
%1行100列的F_kx1\F_kx2\F_kx3\F_kx4，转置成100行1列Fkx1\Fkx2\Fkx3\Fkx4
Fkx1_weighted=F_kx1_weighted';Fkx2_weighted=F_kx2_weighted';Fkx3_weighted=F_kx3_weighted';Fkx4_weighted=F_kx4_weighted';
%等间隔选取的的100个点构成的向量值x1c的权重核分布估计函数：Fkx1_weighted、Fkx2_weighted、Fkx3_weighted、 Fkx4_weighted(其为n行1列)联立成4维随机向量，即100*4矩阵
%联立成100*4矩阵
U_weighted=[Fkx1_weighted(:) Fkx2_weighted(:) Fkx3_weighted(:) Fkx4_weighted(:)];

%各个测点残差的核分布函数值F_kx1q_weighted、F_kx2q_weighted、F_kx3q_weighted、 F_kx4q_weighted(其为n行1列)联立成4维随机向量，即n*4矩阵
%联立成n*4矩阵
T_weighted=[F_kx1q_weighted ,F_kx2q_weighted, F_kx3q_weighted, F_kx4q_weighted];








%%


% 1、四维Gumbel Copula函数
%四维Gumbel Copula函数形状参数求解
alpha_gumbel4D_weighted = estimate_gumbel_4D_alpha(U_weighted);
%定义极大化似然函数的目标函数
disp(['形状参数的极大似然估计值为：', num2str(alpha_gumbel4D_weighted)])

%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel的AIC 、BIC
[aic_k_gumbel4D_weighted, bic_k_gumbel4D_weighted] = compute_gumbel4D_aic_bic(U_weighted, alpha_gumbel4D_weighted);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel概率分布函数值与概率密度函数值
[cdf_k_gumbel4D_weighted, pdf_k_gumbel4D_weighted] = gumbel_copula_4d_non_recursive( Fkx1_weighted,Fkx2_weighted, Fkx3_weighted, Fkx4_weighted, alpha_gumbel4D_weighted);

%求解核分布估计函数值F_kx1q\F_kx2q\F_kx3q\F_kx4q的Gumbel的AIC 、BIC
[aic_gumbel4D_weighted, bic_gumbel4D_weighted] = compute_gumbel4D_aic_bic(T_weighted, alpha_gumbel4D_weighted);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的Gumbel函数概率分布函数值与概率密度函数值
[cdf_gumbel4D_weighted, pdf_gumbel4D_weighted] = gumbel_copula_4d_non_recursive( F_kx1q_weighted,F_kx2q_weighted, F_kx3q_weighted, F_kx4q_weighted, alpha_gumbel4D_weighted);


%%

% 2、四维Frank Copula函数
%四维Frank Copula函数形状参数求解
alpha_frank4D_weighted = estimate_frank4D_alpha(U_weighted);
%定义极大化似然函数的目标函数
disp(['形状参数的极大似然估计值为：', num2str(alpha_frank4D_weighted)])

%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Frank的AIC 、BIC
[aic_k_frank4D_weighted, bic_k_frank4D_weighted] = frank_copula_4d_aic_bic(U_weighted, alpha_frank4D_weighted);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的frank概率分布函数值与概率密度函数值
[cdf_k_frank4D_weighted, pdf_k_frank4D_weighted] = frank_copula_4d( Fkx1_weighted,Fkx2_weighted, Fkx3_weighted, Fkx4_weighted, alpha_frank4D_weighted);

%求解核分布估计函数值F_kx1q\F_kx2q\F_kx3q\F_kx4q的Frank的AIC 、BIC
[aic_frank4D_weighted, bic_frank4D_weighted] = frank_copula_4d_aic_bic(T_weighted, alpha_frank4D_weighted);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的frank函数概率分布函数值与概率密度函数值
[cdf_frank4D_weighted, pdf_frank4D_weighted] =frank_copula_4d( F_kx1q_weighted,F_kx2q_weighted, F_kx3q_weighted, F_kx4q_weighted, alpha_frank4D_weighted);


%%
% 3、四维Clayton Copula函数
%四维Clayton Copula函数形状参数求解
alpha_clayton4D_weighted = estimate_clayton_copula_4d_theta(U_weighted);
%定义极大化似然函数的目标函数
disp(['形状参数的极大似然估计值为：', num2str(alpha_clayton4D_weighted)])

%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Clayton的AIC 、BIC
[aic_k_clayton4D_weighted, bic_k_clayton4D_weighted] = clayton_copula_4d_aic_bic(U_weighted, alpha_clayton4D_weighted);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的clayton概率分布函数值与概率密度函数值
[cdf_k_clayton4D_weighted, pdf_k_clayton4D_weighted] = clayton_copula_4d( Fkx1_weighted,Fkx2_weighted, Fkx3_weighted, Fkx4_weighted, alpha_clayton4D_weighted);

%求解核分布估计函数值F_kx1q\F_kx2q\F_kx3q\F_kx4q的Clyaton的AIC 、BIC
[aic_clayton4D_weighted, bic_clayton4D_weighted] = clayton_copula_4d_aic_bic(T_weighted, alpha_clayton4D_weighted);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的clayton函数概率分布函数值与概率密度函数值
[cdf_clayton4D_weighted, pdf_clayton4D_weighted] =clayton_copula_4d( F_kx1q_weighted,F_kx2q_weighted, F_kx3q_weighted, F_kx4q_weighted, alpha_clayton4D_weighted);


