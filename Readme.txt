Overview
------------------
This is an extension of the Mixed Freq. Data Sampling model ([Ghysels 2004](https://econpapers.repec.org/paper/circirwor/2004s-20.htm)) to multivariate settings.
Yet, this is still working-in-progress since the proposed model cannot find visible solution to parse matricies.


Function name
Input/ Description 
Output
mMIDAS_master
    • This is the main function

m1MIDAS_ADL
    • This function generates the estimation result.

    • Inputs:
[DataY,DataYdate];
[DataX,DataXdate]
    • Explanatory variable and regressors;
    • Date vector
Xlag;
Ylag
    • No. of lags in high-and-low freq. variables
Horizon = 3;

    • MIDAS lead/lag specification
[EstStart, EstEnd];
    • The start and end date of the estimation
※
Estimation result
mMixFreqData.m
This function aligns the date of both high-and-low freq. variables given the specified lag numbers
    • 3-dimensional matrix of regressors 
    • Vector of explanatory variable,
MIDAS_Estimate.m 
    • This is the main estimation algorithm which depends on the pre-specified polynomial function

midas_X.m
    • This function maps high-frequency data into low-frequency vector for a given pair of hyperparameters .
    • Matrix of 
    • Matrix of weights 
mfrvobj_adl.m;
Jacob.m
    • These functions are the objective functions which govern the non-linear least square (lsqnonlin.m);
    • mfrvobj_adl.m: generates the estimation residuals for given pairs of .
    • Jacob.m: constructs the Jacobian matrix,  
    • Residuals 
    • Jacobian matrix 
    • Matrix of optimized weights 

Mathematical details of the model can be found in Report.pdf

For any technical issues related to this package, please contact: fokmchubert@yahoo.com.hk

Thanks.
Hubert FMC
