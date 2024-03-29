# copula-function
A data-driven early warning indicator based on copula function and ISSR-MDF model

1. t1 is the kernel density function and kernel distribution function used to construct the marginal density function and marginal distribution function files for the residuals of each measurement point.
2. t2 is the process of obtaining t copula and Gaussian copula density function and distribution function without using ISSR-MDF model, which also involves the calculation of AIC BIC.
3. t3 is a file of density and distribution functions of Gumbel copula, Frank copula and Clayton copula without ISSR-MDF model, which also involves the calculation of AIC BIC.
4. t4 is the process of using ISSR-MDF model to obtain the density and distribution functions of t copula and Gaussian copula, which also involves the calculation of AIC BIC.
5. t5" is the ISSR-MDF model to obtain the density function and distribution function of Gumbel copula, Frank copula and Clayton copula, which also involves the calculation of AIC BIC.
6. EXCEL "Monitoring point residual data_1" is the processed monitoring data of measuring points and the comprehensive early warning indicators not constructed by ISSR-MDF model.
7. EXCEL "Monitoring point residual data_2" is the processed monitoring data and the comprehensive early warning indicators constructed using the ISSR-MDF model, which also includes the standard deviation of the residuals of the measurement points, the safety stability rate and hazard weights of the measurement points, as well as the resulting new random variables.
8. The rest of the files are copula function files, which are directly called from the "t" file.
