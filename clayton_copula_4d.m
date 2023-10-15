function [C_4D, c_4D] = clayton_copula_4d(u1, u2, u3, u4, theta)
% 参数：
%   - u1, u2, u3, u4：边缘分布函数值，长度为n的向量
%   - theta：已知的Clayton Copula参数
% 返回值：
%   - C_4D：四维Clayton Copula函数的CDF值
%   - c_4D：四维Clayton Copula函数的PDF值

% 二维Clayton Copula函数CDF的表达式
C_2D = @(u, v, theta) max((u.^(-theta) + v.^(-theta) - 1), 0).^(-1/theta);

% 四维Clayton Copula函数CDF的表达式
A = C_2D(u1, u2, theta);
B = C_2D(u3, u4, theta);
C_4D = C_2D(A, B, theta);

% 二维Clayton Copula函数PDF的表达式
c_2D = @(u, v, theta) ((theta + 1) .* (u .* v).^(-theta - 1) .* (u.^(-theta) + v.^(-theta) - 1).^(-theta - 2)) ./ (u .* v);
c_2D_u1_u2 = c_2D(u1, u2, theta);
c_2D_u3_u4 = c_2D(u3, u4, theta);

% 四维Clayton Copula函数PDF的表达式
c_4D = c_2D(A, B, theta) .* c_2D_u1_u2 .* c_2D_u3_u4;
end
