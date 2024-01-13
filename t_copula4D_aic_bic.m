function [AIC, BIC] = t_copula4D_aic_bic(u, R, nu)
% Parameters:
%   - u: n rows by 4 columns matrix of marginal distribution function values
%   - R: 4x4 correlation matrix
%   - nu: Degrees of freedom for the t Copula function
% Returns:
%   - AIC: AIC value for the 4-dimensional t Copula function
%   - BIC: BIC value for the 4-dimensional t Copula function

n = size(u, 1); % Number of data points

% Calculate log-likelihood value
log_likelihood = sum(log(copulapdf('t', u, R, nu)), 'all');

% Calculate the number of parameters (number of upper triangular elements in the correlation matrix + degrees of freedom parameter)
num_params = 4 * (4 - 1) / 2 + 1;

% Calculate AIC and BIC values
AIC = -2 * log_likelihood + 2 * num_params;
BIC = -2 * log_likelihood + num_params * log(n);
end
