function [AIC, BIC] = clayton_copula_4d_aic_bic(u, theta_hat)
% Parameters:
%   - u: An n-by-4 matrix of values from the marginal distribution functions
%   - theta_hat: Estimated parameter θ for the four-dimensional Clayton Copula
% Return values:
%   - AIC: Akaike Information Criterion value for the four-dimensional Clayton Copula function
%   - BIC: Bayesian Information Criterion value for the four-dimensional Clayton Copula function

% Get the number of rows (sample size) and columns (dimensions) of the data
[n, d] = size(u);

% Check if the input data is four-dimensional
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% Calculate the probability density function values for the four-dimensional Clayton Copula function
[~, c_4D] = clayton_copula_4d(u(:, 1), u(:, 2), u(:, 3), u(:, 4), theta_hat);

% Calculate the log-likelihood
log_likelihood = sum(log(c_4D));

% Calculate AIC and BIC
k = 1; % Number of parameters
AIC = -2 * log_likelihood + 2 * k;
BIC = -2 * log_likelihood + log(n) * k;
end

