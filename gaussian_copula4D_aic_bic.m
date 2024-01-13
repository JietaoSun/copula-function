function [AIC, BIC] = gaussian_copula4D_aic_bic(u, R)
% Parameters:
%   - u: An n-by-4 matrix of values from the marginal distribution functions
%   - R: 4x4 covariance matrix
% Return values:
%   - AIC: Akaike Information Criterion value for the four-dimensional Gaussian Copula function
%   - BIC: Bayesian Information Criterion value for the four-dimensional Gaussian Copula function

n = size(u, 1); % Number of data points

% Calculate the standard normal distribution quantiles from the marginal distribution function values
% z = norminv(u);

% Calculate the log-likelihood
log_likelihood = sum(log(copulapdf('Gaussian', u, R)),'all');

% Calculate the number of parameters (number of elements in the upper triangle of the covariance matrix)
num_params = 4 * (4 - 1) / 2;

% Calculate AIC and BIC values
AIC = -2 * log_likelihood + 2 * num_params;
BIC = -2 * log_likelihood + num_params * log(n);
end
