function [AIC, BIC] = gaussian_copula4D_aic_bic(u, R)
% 参数：
%   - u：边缘分布函数值构成的n行4列矩阵
%   - R：4x4协方差矩阵
% 返回值：
%   - AIC：四维Gaussian Copula函数的AIC值
%   - BIC：四维Gaussian Copula函数的BIC值

n = size(u, 1); % 数据点数量

% % 通过边缘分布函数值计算标准正态分布分位数
% z = norminv(u);

% 计算对数似然值
log_likelihood = sum(log(copulapdf('Gaussian', u, R)),'all');

% 计算参数数量（协方差矩阵的上三角元素数量）
num_params = 4 * (4 - 1) / 2;

% 计算AIC和BIC值
AIC = -2 * log_likelihood + 2 * num_params;
BIC = -2 * log_likelihood + num_params * log(n);
end
