function alpha_hat = estimate_frank4D_alpha(u)
% 参数：
%   - u：边缘分布函数值构成的n行4列矩阵
% 返回值：
%   - alpha_hat：估计的四维Frank Copula参数α

% 获取数据的行数和列数
[~, d] = size(u);

% 判断输入数据是否为四维
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% 定义对数似然函数
log_likelihood = @(alpha) -sum(log(frank_copula_4d_pdf(u(:,1), u(:,2), u(:,3), u(:,4), alpha)));

% 优化对数似然函数以估计参数α
options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000);
alpha0 = 0.1; % 初始参数值
alpha_lb = 0; % 参数下限
alpha_ub = Inf; % 参数上限
alpha_hat = fmincon(log_likelihood, alpha0, [], [], [], [], alpha_lb, alpha_ub, [], options);
end

