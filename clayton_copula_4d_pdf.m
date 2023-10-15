function c_4D = clayton_copula_4d_pdf(u1, u2, u3, u4, alpha)
% 参数：
%   - u1, u2, u3, u4：边缘分布函数值，长度为n的向量
%   - alpha：已知的Clayton Copula参数
% 返回值：
%   - c_4D：四维Clayton Copula函数的PDF值

% 二维Clayton Copula函数CDF的表达式
C_2D = @(u, v, alpha) (max(u.^(-alpha) + v.^(-alpha) - 1, 0)).^(-1/alpha);

% 四维Clayton Copula函数CDF的表达式
A = C_2D(u1, u2, alpha);
B = C_2D(u3, u4, alpha);

% 二维Clayton Copula函数PDF的表达式
c_2D = @(u, v, alpha) ((alpha + 1) * (u .* v).^(-alpha - 1) .* (u.^(-alpha) + v.^(-alpha) - 1).^(-2 - 1/alpha)) .* (1 + (u.^(-alpha) + v.^(-alpha) - 1).^(-1/alpha));

c_2D_u1_u2 = c_2D(u1, u2, alpha);
c_2D_u3_u4 = c_2D(u3, u4, alpha);

% 四维Clayton Copula函数PDF的表达式
c_4D = c_2D(A, B, alpha) .* c_2D_u1_u2 .* c_2D_u3_u4;
end
