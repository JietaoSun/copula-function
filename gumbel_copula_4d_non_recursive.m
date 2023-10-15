function [C_4D, c_4D] = gumbel_copula_4d_non_recursive(u1, u2, u3, u4, alpha)
% 参数：
%   - u1, u2, u3, u4：边缘分布函数值，长度为n的向量
%   - alpha：已知的Gumbel Copula参数
% 返回值：
%   - C_4D：四维Gumbel Copula函数的CDF值
%   - c_4D：四维Gumbel Copula函数的PDF值

% 二维Gumbel Copula函数CDF的表达式
C_2D = @(u, v, alpha) exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1/alpha));

% 四维Gumbel Copula函数CDF的表达式
A = C_2D(u1, u2, alpha);
B = C_2D(u3, u4, alpha);
C_4D = C_2D(A, B, alpha);

% 二维Gumbel Copula函数PDF的表达式
c_2D = @(u, v, alpha) (alpha ./ ((-log(u)) .* (-log(v)))).*(u.*v).^(-alpha) .* (((-log(u)).^alpha + (-log(v)).^alpha).^(1/alpha - 1)) .* exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1/alpha));
c_2D_u1_u2 = c_2D(u1, u2, alpha);
c_2D_u3_u4 = c_2D(u3, u4, alpha);

% 四维Gumbel Copula函数PDF的表达式
c_4D = c_2D(A, B, alpha) .* c_2D_u1_u2 .* c_2D_u3_u4;
end
