# Parsimonious-MV-HMRMs

This repository contains the code for fitting 98 parsimonious MV-HMRMFCs and 9604 parsimonious MV-HMRMRCs. In the following, you find a description of the functions (and their arguments) contained in the Functions.R file.

## Eigen.CWM_HMM_init ##

### Description ###

Runs the initialization of the EMC algorithm used for fitting the 9604 parsimonious MV-HMRMRCs. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

Eigen.CWM_HMM_init (Y, X, k, mod.row.Y = "all", mod.col.Y = "all", mod.row.X = "all", mod.col.X = "all", nstartR = 50, nThreads = 1)

### Arguments ###

* Y: An array of dimensions P x R x I x T, where P and R are the rows and columns of the responses, respectively, I refers to the statistical units, and T refers to the time points.
* X: An array of dimensions Q x R x I x T, where Q and R are the rows and columns of the covariates, respectively, I refers to the statistical units, and T refers to the time points. 
* k: A vector containing the numbers of groups to be tried. It must start from 1. 
* mod.row.Y: A character vector indicating the parsimonious structure for the row covariance matrix of the responses. Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV" or "all". When "all" is used, all of the 14 parsimonious structures are considered. 
* mod.col.Y: A character vector indicating the parsimonious structure for the column covariance matrix of the responses. Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV" or "all". When "all" is used, all of the 7 parsimonious structures are considered. 
* mod.row.X: A character vector indicating the parsimonious structure for the row covariance matrix of the covariates. Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV" or "all". When "all" is used, all of the 14 parsimonious structures are considered. 
* mod.col.X: A character vector indicating the parsimonious structure for the column covariance matrix of the covariates. Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV" or "all". When "all" is used, all of the 7 parsimonious structures are considered. 
* nstartR: An integer specifying the number of random starts to be considered. Default value is 50. 
* nThreads: A positive integer indicating the number of cores used for running in parallel.  

## Eigen.CWM_HMM_fit ##

### Description ###

Fits, by using the EMC algorithm, the 9604 parsimonious MV-HMRMRCs. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

Eigen.CWM_HMM_fit (Y, X, init.par = NULL, tol = 0.001, maxit = 500, nThreads = 1)

### Arguments ###

* Y: An array of dimensions P x R x I x T, where P and R are the rows and columns of the responses, respectively, I refers to the statistical units, and T refers to the time points.
* X: An array of dimensions Q x R x I x T, where Q and R are the rows and columns of the covariates, respectively, I refers to the statistical units, and T refers to the time points.
* init.par: The output of the Eigen.CWM_HMM_init() function for initializing the fitting algorithm.
* tol: Threshold for ECM algorithm convergence. Default value is 0.001.
* maxit: Maximum number of iterations for the ECM algorithm. Default value is 500.
* nThreads: A positive integer indicating the number of cores used for running in parallel.

## extract.bestM ##

### Description ###

This function extracts the top/best fitting models according to the Bayesian information criterion (BIC).

### Usage ###

extract.bestM (results, top = 1)

### Arguments ###

* results: The output of the Eigen.CWM_HMM_fit() function.
* top: A number indicating how many models to extract from the ranking provided by the BIC. 

## extract.shortEM ##

### Description ###

This function extracts the names and k of the best fitting models provided by the extract.bestM function.

### Usage ###

extract.shortEM (Extract.bestM)

### Arguments ###

* Extract.bestM: The output of the function extract.bestM().

## filt.init ##

### Description ###

This function extracts only the initializations to be used in the second step of our fitting strategy, according to the results provided by the extract.shortEM() function.

### Usage ###

filt.init (res.init, Extract.bestM)

### Arguments ###

* res.init: The output of the Eigen.CWM_HMM_init() function.
* Extract.bestM: The output of the extract.shortEM() function.

