function c_4D = frank_copula_4d_pdf(u1, u2, u3, u4, alpha)
% Parameters:
%   - u1, u2, u3, u4: Values of the marginal distribution functions, vectors of length n
%   - alpha: Known parameter of the Frank Copula
% Return values:
%   - c_4D: PDF values of the four-dimensional Frank Copula function

% Expression for the CDF of the two-dimensional Frank Copula function
C_2D = @(u, v, alpha) -(1/alpha)*log(1 + ((exp(-alpha*u)-1).*(exp(-alpha*v)-1))./(exp(-alpha)-1));

% Expression for the CDF of the four-dimensional Frank Copula function
A = C_2D(u1, u2, alpha);
B = C_2D(u3, u4, alpha);

% Expression for the PDF of the two-dimensional Frank Copula function
c_2D = @(u, v, alpha) (alpha*exp(alpha*(u+v-2*u.*v)))./(((exp(alpha)-1)*(exp(alpha*u)+exp(alpha*v)-2*exp(alpha*u.*v)-1)).^2);
c_2D_u1_u2 = c_2D(u1, u2, alpha);
c_2D_u3_u4 = c_2D(u3, u4, alpha);

% Expression for the PDF of the four-dimensional Frank Copula function
c_4D = c_2D(A, B, alpha) .* c_2D_u1_u2 .* c_2D_u3_u4;
end

