

% t1产生的各个边缘分布函数再一次陈述
%残差数据Xi的取值范围内等间隔选取的的100个点构成的向量形成的核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4,联立成100*4矩阵
%1行n列的F_kx1\F_kx2\F_kx3\F_kx4，转置成n行1列Fkx1\Fkx2\Fkx3\Fkx4
Fkx1=F_kx1';Fkx2=F_kx2';Fkx3=F_kx3';Fkx4=F_kx4';
U=[F_kx1(:),F_kx2(:), F_kx3(:),F_kx4(:)];
%各个测点残差的核分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q联立成4维随机向量
T=[F_kx1q ,F_kx2q, F_kx3q, F_kx4q];
 

%Gaussian t这两种Copula函数的未知参数本质上是Pearson线性相关系数，调用函数:copulafit）
% 1、计算多维Gaussian Copula

%多维Gaussian Copula的未知参数矩阵（Pearson线性相关系数矩阵）求解
Gaussian_Pearson=copulafit('Gaussian',U); %Gaussian_Pearson为（Pearson线性相关系数）未知参数矩阵
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel的AIC 、BIC
[aic_k_gaussian4D, bic_k_gaussian4D] = gaussian_copula4D_aic_bic(U, Gaussian_Pearson);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel概率分布函数值与概率密度函数值
F_k_Gaussian=copulacdf('Gaussian',U,Gaussian_Pearson); 
f_k_Gaussian=copulapdf('Gaussian',U,Gaussian_Pearson);

%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的AIC 、BIC
[aic_gaussian4D, bic_gaussian4D] = gaussian_copula4D_aic_bic(T, Gaussian_Pearson);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的Gumbel函数概率分布函数值与概率密度函数值
F_Gaussian=copulacdf('Gaussian',T,Gaussian_Pearson);
f_Gaussian=copulapdf('Gaussian',T,Gaussian_Pearson);
%%

%2 计算多维t Copula

%多维的未知参数矩阵（Pearson线性相关系数矩阵）求解
[t_Pearson,nuhat]=copulafit('t',U) ;%nuhat为tCopula函数的自由度。注意：数据量太小或数据间依赖性过强会出不来结果,t_Pearson为未知参数矩阵
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel的AIC 、BIC
[aic_k_t4D, bic_k_t4D] = t_copula4D_aic_bic(U, t_Pearson,nuhat);
%求解核分布估计函数值F_kx1\F_kx2\F_kx3\F_kx4的Gumbel概率分布函数值与概率密度函数值
F_k_t=copulacdf('t',U,t_Pearson,nuhat); 
f_k_t=copulapdf('t',U,t_Pearson,nuhat);
%%
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的AIC 、BIC
[aic_t4D, bic_t4D] = t_copula4D_aic_bic(T, t_Pearson,nuhat);
%求解观测值对应的边缘分布函数值F_kx1q、F_kx2q、F_kx3q、 F_kx4q的Gumbel函数概率分布函数值与概率密度函数值
F_t=copulacdf('t',T,t_Pearson,nuhat); 
f_t=copulapdf('t',T,t_Pearson,nuhat);





