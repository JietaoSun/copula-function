
%残差数据Xi的取值范围内等间隔选取的的100个点构成的向量形成的核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4,联立成100*4矩阵
U=[F_kx1(:),F_kx2(:), F_kx3(:),F_kx4(:)];
%1行n列转置成n行1列
Fkx1=F_kx1';Fkx2=F_kx2';Fkx3=F_kx3';Fkx4=F_kx4';
%各个测点残差的核分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q联立成4维随机向量，即  n*4矩阵
T=[F_kx1q ,F_kx2q, F_kx3q, F_kx4q];
%%


% 1、四维Gumbel Copula函数
%四维Gumbel Copula函数形状参数求解
alpha_gumbel4D = estimate_gumbel_4D_alpha(U);
%定义极大化似然函数的目标函数
disp(['形状参数的极大似然估计值为：', num2str(alpha_gumbel4D)])

%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel的AIC 、BIC
[aic_k_gumbel4D, bic_k_gumbel4D] = compute_gumbel4D_aic_bic(U, alpha_gumbel4D);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel概率分布函数值与概率密度函数值
[cdf_k_gumbel4D, pdf_k_gumbel4D] = gumbel_copula_4d_non_recursive( Fkx1,Fkx2, Fkx3, Fkx4, alpha_gumbel4D);

%求解核分布估计函数值F_kx1q\F_kx2q\F_kx3q\F_kx4q的Gumbel的AIC 、BIC
[aic_gumbel4D, bic_gumbel4D] = compute_gumbel4D_aic_bic(T, alpha_gumbel4D);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的Gumbel函数概率分布函数值与概率密度函数值
[cdf_gumbel4D, pdf_gumbel4D] = gumbel_copula_4d_non_recursive( F_kx1q,F_kx2q, F_kx3q, F_kx4q, alpha_gumbel4D);


%%

% 2、四维Frank Copula函数
%四维Frank Copula函数形状参数求解
alpha_frank4D = estimate_frank4D_alpha(U);
%定义极大化似然函数的目标函数
disp(['形状参数的极大似然估计值为：', num2str(alpha_frank4D)])

%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Frank的AIC 、BIC
[aic_k_frank4D, bic_k_frank4D] = frank_copula_4d_aic_bic(U, alpha_frank4D);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的frank概率分布函数值与概率密度函数值
[cdf_k_frank4D, pdf_k_frank4D] = frank_copula_4d( Fkx1,Fkx2, Fkx3, Fkx4, alpha_frank4D);

%求解核分布估计函数值F_kx1q\F_kx2q\F_kx3q\F_kx4q的Frank的AIC 、BIC
[aic_frank4D, bic_frank4D] = frank_copula_4d_aic_bic(T, alpha_frank4D);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的frank函数概率分布函数值与概率密度函数值
[cdf_frank4D, pdf_frank4D] =frank_copula_4d( F_kx1q,F_kx2q, F_kx3q, F_kx4q, alpha_frank4D);


%%
% 3、四维Clayton Copula函数
%四维Clayton Copula函数形状参数求解
alpha_clayton4D = estimate_clayton_copula_4d_theta(U);
%定义极大化似然函数的目标函数
disp(['形状参数的极大似然估计值为：', num2str(alpha_clayton4D)])

%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Clayton的AIC 、BIC
[aic_k_clayton4D, bic_k_clayton4D] = clayton_copula_4d_aic_bic(U, alpha_clayton4D);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的clayton概率分布函数值与概率密度函数值
[cdf_k_clayton4D, pdf_k_clayton4D] = clayton_copula_4d( Fkx1,Fkx2, Fkx3, Fkx4, alpha_clayton4D);

%求解核分布估计函数值F_kx1q\F_kx2q\F_kx3q\F_kx4q的Clyaton的AIC 、BIC
[aic_clayton4D, bic_clayton4D] = clayton_copula_4d_aic_bic(T, alpha_clayton4D);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的clayton函数概率分布函数值与概率密度函数值
[cdf_clayton4D, pdf_clayton4D] =clayton_copula_4d( F_kx1q,F_kx2q, F_kx3q, F_kx4q, alpha_clayton4D);


