#Overview
------------------
This is an extension of the Mixed Freq. Data Sampling model ([Ghysels 2004](https://econpapers.repec.org/paper/circirwor/2004s-20.htm)) to multivariate settings.
Yet, this is still working-in-progress since the proposed model cannot find visible solution to parse matricies.

## Function description
--------------------------
mMIDAS_master.m  : This is the main function
mMixFreqData.m   : This function aligns the date of both high-and-low freq. variables given the specified lag numbers
MIDAS_Estimate.m :This is the main estimation algorithm which depends on the pre-specified polynomial function
midas_X.m        : This function maps high-frequency data into low-frequency vector for a given pair of hyperparameters (i.e. Theta 1 & 2).
mfrvobj_adl.m    :These functions are the objective functions which govern the non-linear least square (lsqnonlin.m);
Jacob.m             • mfrvobj_adl.m: generates the estimation residuals for given pairs of .
                    • Jacob.m: constructs the Jacobian matrix
                   
Mathematical details of the model can be found in Report.pdf

For any technical issues related to this package, please contact: fokmchubert@yahoo.com.hk

Thanks.
Hubert FMC
