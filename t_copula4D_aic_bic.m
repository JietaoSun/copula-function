function [AIC, BIC] = t_copula4D_aic_bic(u, R, nu)
% 参数：
%   - u：边缘分布函数值构成的n行4列矩阵
%   - R：4x4相关矩阵
%   - nu：t Copula函数的自由度
% 返回值：
%   - AIC：四维t Copula函数的AIC值
%   - BIC：四维t Copula函数的BIC值

n = size(u, 1); % 数据点数量

% % 通过边缘分布函数值计算t分布分位数
% z = tinv(u, nu);

% 计算对数似然值
log_likelihood = sum(log(copulapdf('t', u, R,nu)),'all');

% 计算参数数量（相关矩阵的上三角元素数量 + 自由度参数）
num_params = 4 * (4 - 1) / 2 + 1;

% 计算AIC和BIC值
AIC = -2 * log_likelihood + 2 * num_params;
BIC = -2 * log_likelihood + num_params * log(n);
end
