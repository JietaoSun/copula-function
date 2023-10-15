function [AIC, BIC] = compute_gumbel4D_aic_bic(u, alpha_hat)
% 参数：
%   - u：边缘分布函数值构成的n行4列矩阵
%   - alpha_hat：估计的四维Gumbel Copula参数α
% 返回值：
%   - AIC：四维Gumbel Copula函数的AIC值
%   - BIC：四维Gumbel Copula函数的BIC值

% 获取数据的行数（样本量）和列数（维度）
[n, d] = size(u);

% 判断输入数据是否为四维
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% 计算四维Gumbel Copula函数的概率密度函数值
[~, c_4D] = gumbel_copula_4d_non_recursive(u(:,1), u(:,2), u(:,3), u(:,4), alpha_hat);

% 计算对数似然值
log_likelihood = sum(log(c_4D));

% 计算AIC和BIC
k = 1; % 参数个数
AIC = -2 * log_likelihood + 2 * k;
BIC = -2 * log_likelihood + log(n) * k;
end
