function theta_hat = estimate_clayton_copula_4d_theta(u)
% 参数：
%   - u：边缘分布函数值构成的n行4列矩阵
% 返回值：
%   - theta_hat：估计的四维Clayton Copula参数θ

% 获取数据的行数（样本量）和列数（维度）
[~, d] = size(u);

% 判断输入数据是否为四维
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% 定义对数似然函数
log_likelihood = @(theta) -sum(log(clayton_copula_4d_pdf(u(:, 1), u(:, 2), u(:, 3), u(:, 4), theta)));

% 设置初始值
theta0 = 0.1;

% 使用优化函数（如fminsearch）求解最大对数似然值
options = optimset('TolX', 1e-6, 'Display', 'off');
theta_hat = fminsearch(log_likelihood, theta0, options);
end
