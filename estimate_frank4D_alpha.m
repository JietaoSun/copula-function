function alpha_hat = estimate_frank4D_alpha(u)
% Parameters:
%   - u: An n-by-4 matrix of values from the marginal distribution functions
% Return values:
%   - alpha_hat: Estimated parameter α for the four-dimensional Frank Copula

% Get the number of rows and columns of the data
[~, d] = size(u);

% Check if the input data is four-dimensional
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% Define the log-likelihood function
log_likelihood = @(alpha) -sum(log(frank_copula_4d_pdf(u(:,1), u(:,2), u(:,3), u(:,4), alpha)));

% Optimize the log-likelihood function to estimate parameter α
options = optimset('Display', 'off', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000);
alpha0 = 0.1; % Initial parameter value
alpha_lb = 0; % Parameter lower bound
alpha_ub = Inf; % Parameter upper bound
alpha_hat = fmincon(log_likelihood, alpha0, [], [], [], [], alpha_lb, alpha_ub, [], options);
end
