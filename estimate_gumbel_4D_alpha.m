function alpha_hat = estimate_gumbel_4D_alpha(u)
% Parameters:
%   - u: An n-by-4 matrix of values from the marginal distribution functions
% Return values:
%   - alpha_hat: Estimated parameter α for the four-dimensional Gumbel Copula

% Get the number of rows and columns of the data
[~, d] = size(u);

% Check if the input data is four-dimensional
if d ~= 4
    error('The input matrix should have 4 columns.');
end

% Define the log-likelihood function
log_likelihood = @(alpha) -sum(log(gumbel_copula_4d_pdf(u(:,1), u(:,2), u(:,3), u(:,4), alpha)));

% Set optimization options
options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point', 'TolX', 1e-6, 'TolFun', 1e-6, 'MaxIter', 1000);

% Set initial parameter value and parameter range
alpha0 = 2; % Initial parameter value
alpha_lb = 1; % Parameter lower bound

% Optimize the log-likelihood function to estimate parameter α
alpha_hat = fmincon(log_likelihood, alpha0, [], [], [], [], alpha_lb, [], [], options);
end

