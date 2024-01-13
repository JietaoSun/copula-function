function [C_4D, c_4D] = gumbel_copula_4d_non_recursive(u1, u2, u3, u4, alpha)
% Parameters:
%   - u1, u2, u3, u4: Values of the marginal distribution functions, vectors of length n
%   - alpha: Known parameter of the Gumbel Copula
% Return values:
%   - C_4D: CDF values of the four-dimensional Gumbel Copula function
%   - c_4D: PDF values of the four-dimensional Gumbel Copula function

% Expression for the CDF of the two-dimensional Gumbel Copula function
C_2D = @(u, v, alpha) exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1/alpha));

% Expression for the CDF of the four-dimensional Gumbel Copula function
A = C_2D(u1, u2, alpha);
B = C_2D(u3, u4, alpha);
C_4D = C_2D(A, B, alpha);

% Expression for the PDF of the two-dimensional Gumbel Copula function
c_2D = @(u, v, alpha) (alpha ./ ((-log(u)) .* (-log(v)))).*(u.*v).^(-alpha) .* (((-log(u)).^alpha + (-log(v)).^alpha).^(1/alpha - 1)) .* exp(-((-log(u)).^alpha + (-log(v)).^alpha).^(1/alpha));
c_2D_u1_u2 = c_2D(u1, u2, alpha);
c_2D_u3_u4 = c_2D(u3, u4, alpha);

% Expression for the PDF of the four-dimensional Gumbel Copula function
c_4D = c_2D(A, B, alpha) .* c_2D_u1_u2 .* c_2D_u3_u4;
end
