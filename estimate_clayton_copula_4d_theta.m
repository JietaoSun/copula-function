function theta_hat = estimate_clayton_copula_4d_theta(u)
% Parameters:
%   - u: An n-by-4 matrix of values from the marginal distribution functions
% Return values:
%   - theta_hat: Estimated parameter θ for the four-dimensional Clayton Copula

% Get the number of rows (sample size) and columns (dimensions) of the data
[~, d] = size(u);

% Check if the input data is four-dimensional
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% Define the log-likelihood function
log_likelihood = @(theta) -sum(log(clayton_copula_4d_pdf(u(:, 1), u(:, 2), u(:, 3), u(:, 4), theta)));

% Set the initial value
theta0 = 0.1;

% Use an optimization function (like fminsearch) to solve for the maximum log-likelihood
options = optimset('TolX', 1e-6, 'Display', 'off');
theta_hat = fminsearch(log_likelihood, theta0, options);
end
