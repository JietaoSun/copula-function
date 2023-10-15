# copula-function
A comprehensive early warning indicator based on copula function and ISSR-MDF model
-t1" is the kernel density function and kernel distribution function used to construct the marginal density function and marginal distribution function files for the residuals of each measurement point.
-t2 is the process of obtaining t copula and Gaussian copula density function and distribution function without using ISSR-MDF model, which also involves the calculation of AIC BIC.
-t3 is a file of density and distribution functions of Gumbel copula, Frank copula and Clayton copula without ISSR-MDF model, which also involves the calculation of AIC BIC.
-t4 is the process of using ISSR-MDF model to obtain the density and distribution functions of t copula and Gaussian copula, which also involves the calculation of AIC BIC.
-t5" is the ISSR-MDF model to obtain the density function and distribution function of Gumbel copula, Frank copula and Clayton copula, which also involves the calculation of AIC BIC.
-EXCEL "监测点残差数据-正式2" is the overall (composite) early warning indicator not constructed by ISSR-MDF model.
EXCEL "监测点残差数据-正式3" is the overall (composite) early warning indicators constructed using the ISSR-MDF model, which also includes the standard deviation of the residuals of each measurement point, the safety and stability rate of each measurement point and the weight of the hazard, and the resulting new random variables.
-The rest of the files are copula function files, which are directly called from the "t" file.
