clear;
clc;
% 1. Extract the residual data of each monitoring point
Monitor=readmatrix('Monitoring point residual data_1.xlsx'); % Note: Read the residual data as a numerical matrix
X1=Monitor(:,1); % Extract the residual data from the 1st column
X2=Monitor(:,2);
X3=Monitor(:,3);
X4=Monitor(:,4); % Extract the residual data from the 4th column

% 2. Preliminary judgment on whether it follows a normal distribution, including skewness (using the skewness function, negative for left skew, i.e., the density function has a long left tail and the peak is to the right; 0 for symmetry; positive for right skew) and kurtosis (using the kurtosis function, normal distribution kurtosis is 3, if greater than 3, the density function curve is steeper than normal distribution near the peak; otherwise flatter than normal distribution) (reflects the steepness of the peak, whether the peak is sharp or flat and the tail is light or heavy)

x1s=skewness(X1); % Calculate the skewness of X1
x2s=skewness(X2); % Calculate the skewness of X2
x3s=skewness(X3); % Calculate the skewness of X3
x4s=skewness(X4); % Calculate the skewness of X4

x1k=kurtosis(X1); % Calculate the kurtosis of X1
x2k=skewness(X2); % Calculate the kurtosis of X2
x3k=skewness(X3); % Calculate the kurtosis of X3
x4k=skewness(X4); % Calculate the kurtosis of X4

% 3. Data normality test to further judge whether it follows a normal distribution (including jbtest (Jarque-Bera test), ktest (Kolmogorov-Smirnov test), lilletest (Lilliefors test),
% If all tests are passed, then h=0, otherwise h=1 (default p=0.01). If passed, use the normal estimate in parameter estimation (normfit function or Matlab statistical toolbox, etc., for unknown parameter estimation). If not passed, use non-parametric estimates such as kernel density estimation, kernel distribution estimation.)

% Perform jbtest (Jarque-Bera test) on Xi, hi=0 means normality test passed, hi=1 means test failed, and gives the corresponding significance level pi
[h1_jbtest,p1_jbtest]=jbtest(X1); 
[h2_jbtest,p2_jbtest]=jbtest(X2); 
[h3_jbtest,p3_jbtest]=jbtest(X3); 
[h4_jbtest,p4_jbtest]=jbtest(X4); 

% Perform ktest (Kolmogorov-Smirnov test) on Xi, hi=0 means normality test passed, hi=1 means test failed, and gives the corresponding significance level pi
% [h1_ktest,p1_ktest]=kstest(X1,[X1,normcdf(X1,mean(X1),std(X1))]);
% [h2_ktest,p2_ktest]=kstest(X2,[X2,normcdf(X2,mean(X2),std(X2))]);
% [h3_ktest,p3_ktest]=kstest(X3,[X3,normcdf(X3,mean(X3),std(X3))]);
% [h4_ktest,p4_ktest]=kstest(X4,[X4,normcdf(X4,mean(X4),std(X4))]);

% Perform lilletest (Lilliefors test) on Xi, hi=0 means normality test passed, hi=1 means test failed, and gives the corresponding significance level pi
[h1_lilletest,p1_lilletest]=lillietest(X1);
[h2_lilletest,p2_lilletest]=lillietest(X2);
[h3_lilletest,p3_lilletest]=lillietest(X3);
[h4_lilletest,p4_lilletest]=lillietest(X4);

% 4. Kernel Density Estimation Function and Kernel Distribution Estimation Function
% Calculation of Kernel Density Estimation Function
[f_kx1, ~, u1]=ksdensity(X1); % x1c is a vector of 100 equally spaced points selected within the range of values in the first column of residual data X1. f_kx1 is the corresponding vector of kernel density estimate values for x1c, and u1 is the optimal bandwidth calculated using the Gaussian kernel function based on the sample size and standard deviation.
[f_kx2, ~, u2]=ksdensity(X2);
[f_kx3, ~, u3]=ksdensity(X3);
[f_kx4, ~, u4]=ksdensity(X4);
u1; % View the optimal bandwidth value for X1
u2;
u3;
u4;
% Kernel Distribution Estimation Function
[F_kx1, x1c]=ksdensity(X1, 'function', 'cdf'); % x1c is a vector of 100 equally spaced points selected within the range of values in the first column of residual data X1. F_kx1 is the corresponding vector of kernel distribution estimate values for x1c.
[F_kx2, x2c]=ksdensity(X2, 'function', 'cdf'); % Kernel distribution estimate values vector F_kx2
[F_kx3, x3c]=ksdensity(X3, 'function', 'cdf'); % Kernel distribution estimate values vector F_kx3
[F_kx4, x4c]=ksdensity(X4, 'function', 'cdf'); % Kernel distribution estimate values vector F_kx4




% 5. Kernel Density Estimation and Fit with Frequency Histograms of Residual Data
% Fit of kernel density estimation function to the frequency distribution histogram of the first column of residual data
figure;
histogram(X1, 30, 'Normalization', 'pdf', 'FaceColor', [142/255 160/255 199/255]);  % Plot the histogram
hold on;
plot(x1c, f_kx1, '--', 'Color', [1/255 1/255 1/255], 'LineWidth', 6); % Plot the kernel density estimation function
xlabel('Residual (mm)', 'FontName', 'Times New Roman', 'FontSize', 56);
ylabel('Frequency/group spacing', 'FontName', 'Times New Roman', 'FontSize', 56);
legend('Frequency distribution histogram', 'Kernel density function {\itf}_1({\itx}_1)', 'Location', 'northeast', 'FontSize', 38, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50, 'LineWidth', 2, 'FontName', 'Times New Roman', 'LooseInset', [0, 0, 0, 0]);


% Fit of kernel density estimation function to the frequency distribution histogram of the second column of residual data
figure;
histogram(X2, 30, 'Normalization', 'pdf', 'FaceColor', [142/255 160/255 199/255]);  % Plot the histogram
hold on;
plot(x2c, f_kx2, '--', 'Color', [0 0 0], 'LineWidth', 6); % Plot the kernel density estimation function
xlabel('Residual (mm)', 'FontName', 'Times New Roman', 'FontSize', 56);
ylabel('Frequency/group spacing', 'FontName', 'Times New Roman', 'FontSize', 56);
legend('Frequency distribution histogram', 'Kernel density function {\itf}_2({\itx}_2)', 'Location', 'northeast', 'FontSize', 40, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50, 'LineWidth', 2, 'FontName', 'Times New Roman', 'LooseInset', [0, 0, 0, 0]);

% Fit of kernel density estimation function to the frequency distribution histogram of the third column of residual data
figure;
histogram(X3, 30, 'Normalization', 'pdf', 'FaceColor', [142/255 160/255 199/255]);  % Plot the histogram
hold on;
plot(x3c, f_kx3, '--', 'Color', [0 0 0], 'LineWidth', 6); % Plot the kernel density estimation function
xlabel('Residual (mm)', 'FontName', 'Times New Roman', 'FontSize', 56);
ylabel('Frequency/group spacing', 'FontName', 'Times New Roman', 'FontSize', 56);
legend('Frequency distribution histogram', 'Kernel density function {\itf}_3({\itx}_3)', 'Location', 'northeast', 'FontSize', 40, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50, 'LineWidth', 2, 'FontName', 'Times New Roman', 'LooseInset', [0, 0, 0, 0]);

% Fit of kernel density estimation function to the frequency distribution histogram of the fourth column of residual data
figure;
histogram(X4, 30, 'Normalization', 'pdf', 'FaceColor', [142/255 160/255 199/255]);  % Plot the histogram
hold on;
plot(x4c, f_kx4, '--', 'Color', [0 0 0], 'LineWidth', 6); % Plot the kernel density estimation function
xlabel('Residual (mm)', 'FontName', 'Times New Roman', 'FontSize', 56);
ylabel('Frequency/group spacing', 'FontName', 'Times New Roman', 'FontSize', 56);
legend('Frequency distribution histogram', 'Kernel density function {\itf}_4({\itx}_4)', 'Location', 'northeast', 'FontSize', 40, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50, 'LineWidth', 2, 'FontName', 'Times New Roman', 'LooseInset', [0, 0, 0, 0]);


% 6. Computing Kernel Density and Distribution Values for Residual Data X1 at Selected Points

% Calculate the kernel density function value at any point xiq for the residual data X1
f_kx1q=interp1(x1c, f_kx1, X1); % The interp1 function is used to compute the kernel density function value f_kx1q at any point for the variable X1.

% Calculate the kernel density function value at any point xiq for the residual data X2
f_kx2q=interp1(x2c, f_kx2, X2);

% Calculate the kernel density function value at any point xiq for the residual data X3
f_kx3q=interp1(x3c, f_kx3, X3);

% Calculate the kernel density function value at any point xiq for the residual data X4
f_kx4q=interp1(x4c, f_kx4, X4);

% Calculate the kernel distribution function value at any point xiq for the residual data X1
F_kx1q=interp1(x1c, F_kx1, X1, 'previous'); % The interp1 function is used to compute the kernel distribution function value F_kx1q at any point for the variable X1, with the 'previous' parameter indicating selecting the nearest value to xiq.

% Calculate the kernel distribution function value at any point xiq for the residual data X2
F_kx2q=interp1(x2c, F_kx2, X2, 'previous');

% Calculate the kernel distribution function value at any point xiq for the residual data X3
F_kx3q=interp1(x3c, F_kx3, X3, 'previous');

% Calculate the kernel distribution function value at any point xiq for the residual data X4
F_kx4q=interp1(x4c, F_kx4, X4, 'previous');

% Combine the kernel distribution function values at various measurement points into a random vector for subsequent feeding into the multidimensional Copula function
T=[F_kx1q, F_kx2q, F_kx3q, F_kx4q]; % T is the random vector formed by combining the kernel distribution function values at various measurement points.

