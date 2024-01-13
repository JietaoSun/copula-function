function [C_4D, c_4D] = frank_copula_4d(u1, u2, u3, u4, alpha)
% Parameters:
%   - u1, u2, u3, u4: Values of the marginal distribution functions, vectors of length n
%   - alpha: Known parameter of the Frank Copula
% Return values:
%   - C_4D: CDF values of the four-dimensional Frank Copula function
%   - c_4D: PDF values of the four-dimensional Frank Copula function

% Expression for the CDF of the two-dimensional Frank Copula function
C_2D = @(u, v, alpha) -(1/alpha)*log(1 + ((exp(-alpha*u)-1).*(exp(-alpha*v)-1))./(exp(-alpha)-1));

% Expression for the CDF of the four-dimensional Frank Copula function
A = C_2D(u1, u2, alpha);
B = C_2D(u3, u4, alpha);
C_4D = C_2D(A, B, alpha);

% Expression for the PDF of the two-dimensional Frank Copula function
c_2D = @(u, v, alpha) (alpha*exp(alpha*(u+v-2*u.*v)))./(((exp(alpha)-1)*(exp(alpha*u)+exp(alpha*v

