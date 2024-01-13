function [C_4D, c_4D] = clayton_copula_4d(u1, u2, u3, u4, theta)
% Parameters:
%   - u1, u2, u3, u4: Values of the marginal distribution functions, vectors of length n
%   - theta: Known parameter of the Clayton Copula
% Return values:
%   - C_4D: CDF values of the four-dimensional Clayton Copula function
%   - c_4D: PDF values of the four-dimensional Clayton Copula function

% Expression for the CDF of the two-dimensional Clayton Copula function
C_2D = @(u, v, theta) max((u.^(-theta) + v.^(-theta) - 1), 0).^(-1/theta);

% Expression for the CDF of the four-dimensional Clayton Copula function
A = C_2D(u1, u2, theta);
B = C_2D(u3, u4, theta);
C_4D = C_2D(A, B, theta);

% Expression for the PDF of the two-dimensional Clayton Copula function
c_2D = @(u, v, theta) ((theta + 1) .* (u .* v).^(-theta - 1) .* (u.^(-theta) + v.^(-theta) - 1).^(-theta - 2)) ./ (u .* v);
c_2D_u1_u2 = c_2D(u1, u2, theta);
c_2D_u3_u4 = c_2D(u3, u4, theta);

% Expression for the PDF of the four-dimensional Clayton Copula function
c_4D = c_2D(A, B, theta) .* c_2D_u1_u2 .* c_2D_u3_u4;
end
