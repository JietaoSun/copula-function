% Forming a 100*4 matrix of kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4 at 100 equally spaced points within the range of values for the residual data Xi
U = [F_kx1(:), F_kx2(:), F_kx3(:), F_kx4(:)];

% Transposing 1-row-n-column matrices F_kx1\F_kx2\F_kx3\F_kx4 to n-row-1-column Fkx1\Fkx2\Fkx3\Fkx4
Fkx1 = F_kx1'; Fkx2 = F_kx2'; Fkx3 = F_kx3'; Fkx4 = F_kx4';

% Combining the kernel distribution function values F_kx1q, F_kx2q, F_kx3q, and F_kx4q for various measurement points into a 4-dimensional random vector, forming an n*4 matrix
T = [F_kx1q, F_kx2q, F_kx3q, F_kx4q];
%%

% 1. Four-Dimensional Gumbel Copula Function
% Estimation of the shape parameter for the four-dimensional Gumbel Copula function
alpha_gumbel4D = estimate_gumbel_4D_alpha(U);
% Defining the objective function for maximizing the likelihood function
disp(['Maximum likelihood estimate of the shape parameter is: ', num2str(alpha_gumbel4D)])

% Calculating AIC and BIC for the Gumbel distribution of kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4
[aic_k_gumbel4D, bic_k_gumbel4D] = compute_gumbel4D_aic_bic(U, alpha_gumbel4D);
% Calculating the cumulative distribution and probability density functions of the Gumbel distribution for kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4
[cdf_k_gumbel4D, pdf_k_gumbel4D] = gumbel_copula_4d_non_recursive(Fkx1, Fkx2, Fkx3, Fkx4, alpha_gumbel4D);

% Calculating AIC and BIC for the observed marginal distribution values F_kx1q, F_kx2q, F_kx3q, and F_kx4q
[aic_gumbel4D, bic_gumbel4D] = compute_gumbel4D_aic_bic(T, alpha_gumbel4D);
% Calculating the cumulative distribution and probability density functions of the Gumbel distribution for observed marginal distribution values F_kx1q, F_kx2q, F_kx3q, and F_kx4q
[cdf_gumbel4D, pdf_gumbel4D] = gumbel_copula_4d_non_recursive(F_kx1q, F_kx2q, F_kx3q, F_kx4q, alpha_gumbel4D);
%%

% 2. Four-Dimensional Frank Copula Function
% Estimation of the shape parameter for the four-dimensional Frank Copula function
alpha_frank4D = estimate_frank4D_alpha(U);
% Defining the objective function for maximizing the likelihood function
disp(['Maximum likelihood estimate of the shape parameter is: ', num2str(alpha_frank4D)])

% Calculating AIC and BIC for the Frank distribution of kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4
[aic_k_frank4D, bic_k_frank4D] = frank_copula_4d_aic_bic(U, alpha_frank4D);
% Calculating the cumulative distribution and probability density functions of the Frank distribution for kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4
[cdf_k_frank4D, pdf_k_frank4D] = frank_copula_4d(Fkx1, Fkx2, Fkx3, Fkx4, alpha_frank4D);

% Calculating AIC and BIC for the observed marginal distribution values F_kx1q, F_kx2q, F_kx3q, and F_kx4q
[aic_frank4D, bic_frank4D] = frank_copula_4d_aic_bic(T, alpha_frank4D);
% Calculating the cumulative distribution and probability density functions of the Frank distribution for observed marginal distribution values F_kx1q, F_kx2q, F_kx3q, and F_kx4q
[cdf_frank4D, pdf_frank4D] = frank_copula_4d(F_kx1q, F_kx2q, F_kx3q, F_kx4q, alpha_frank4D);
%%

% 3. Four-Dimensional Clayton Copula Function
% Estimation of the shape parameter for the four-dimensional Clayton Copula function
alpha_clayton4D = estimate_clayton_copula_4d_theta(U);
% Defining the objective function for maximizing the likelihood function
disp(['Maximum likelihood estimate of the shape parameter is: ', num2str(alpha_clayton4D)])

% Calculating AIC and BIC for the Clayton distribution of kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4
[aic_k_clayton4D, bic_k_clayton4D] = clayton_copula_4d_aic_bic(U, alpha_clayton4D);
% Calculating the cumulative distribution and probability density functions of the Clayton distribution for kernel distribution estimate values F_kx1\F_kx2\F_kx3\F_kx4
[cdf_k_clayton4D, pdf_k_clayton4D] = clayton_copula_4d(Fkx1, Fkx2, Fkx3, Fkx4, alpha_clayton4D);

% Calculating AIC and BIC for the observed marginal distribution values F_kx1q, F_kx2q, F_kx3q, and F_kx4q
[aic_clayton4D, bic_clayton4D] = clayton_copula_4d_aic_bic(T, alpha_clayton4D);
% Calculating the cumulative distribution and probability density functions of the Clayton distribution for observed marginal distribution values F_kx1q, F_kx2q, F_kx3q, and F_kx4q
[cdf_clayton4D, pdf_clayton4D] = clayton_copula_4d(F_kx1q, F_kx2q, F_kx3q, F_kx4q, alpha_clayton4D);
