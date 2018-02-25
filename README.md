# Fitting-the-yield-curve-with-Kalman-Filter

This code replicates the Kalman filter algorithm for one factor CIR model as described in 
Bolder, David Jamieson, Affine Term-Structure Models: Theory and Implementation (October 2001). 
Available at SSRN: https://ssrn.com/abstract=1082826 or http://dx.doi.org/10.2139/ssrn.1082826

1. SimulateCIR.R : This function simulates the CIR term structure of interest rate with specified CIR parameters. The function returns
monthly zero coupon rates for a period of 10 years.

2. LogLikelihoodCIR.R : This function returns the loglikelihood values of CIR parameters given the zero coupon rates.

3. MinimizeLL : This is the main function which calls the above two functions. The optimization function is iterated 
200 times and the estimated CIR parameters is the mean of the 200 observations.

Result: Estimating the term structure can be difficult as the instantaneous short rates are unobservable. Kalman filter
is an optimal estimation algorithm used to estimate states of a system when indirect measurement are present. Simulated and
estimated CIR paramters obtained from the code are close, indicating Kalman filter to be an effective tool for estimating the
term structure of interest rate. 
