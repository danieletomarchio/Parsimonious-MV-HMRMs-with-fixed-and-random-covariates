This repository contains the code for fitting the 9604 parsimonious MV-HMRMRCs and 98 parsimonious MV-HMRMFCs. In the following, you find a description of the functions (and their arguments) contained in the MV-HMRMRCs.R and MV-HMRMFCs files, respectively.

# Parsimonious MV-HMRMRCs

## Eigen.CWM_HMM_init ##

### Description ###

Runs the initialization of the ECM algorithm used for fitting the 9604 parsimonious MV-HMRMRCs. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

Eigen.CWM_HMM_init (Y, X, k, mod.row.Y = "all", mod.col.Y = "all", mod.row.X = "all", mod.col.X = "all", nstartR = 50, nThreads = 1)

### Arguments ###

* Y: An array of dimensions P x R x I x T, where P and R are the rows and columns of the responses, respectively, I refers to the statistical units, and T refers to the time points.
* X: An array of dimensions Q x R x I x T, where Q and R are the rows and columns of the covariates, respectively, I refers to the statistical units, and T refers to the time points. 
* k: A vector (or a number) containing the groups to be tried.
* mod.row.Y: A character vector indicating the parsimonious structure for the row covariance matrix of the responses. Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV" or "all". When "all" is used, all of the 14 parsimonious structures are considered. 
* mod.col.Y: A character vector indicating the parsimonious structure for the column covariance matrix of the responses. Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV" or "all". When "all" is used, all of the 7 parsimonious structures are considered. 
* mod.row.X: A character vector indicating the parsimonious structure for the row covariance matrix of the covariates. Possible values are: "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV" or "all". When "all" is used, all of the 14 parsimonious structures are considered. 
* mod.col.X: A character vector indicating the parsimonious structure for the column covariance matrix of the covariates. Possible values are: "II", "EI", "VI", "EE", "VE", "EV", "VV" or "all". When "all" is used, all of the 7 parsimonious structures are considered. 
* nstartR: An integer specifying the number of random starts to be considered. Default value is 50. 
* nThreads: A positive integer indicating the number of cores used for running in parallel.  

## Eigen.CWM_HMM_fit ##

### Description ###

Fits, by using the ECM algorithm, the 9604 parsimonious MV-HMRMRCs. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

Eigen.CWM_HMM_fit (Y, X, init.par = NULL, tol = 0.001, maxit = 500, nThreads = 1)

### Arguments ###

* Y: An array of dimensions P x R x I x T, where P and R are the rows and columns of the responses, respectively, I refers to the statistical units, and T refers to the time points.
* X: An array of dimensions Q x R x I x T, where Q and R are the rows and columns of the covariates, respectively, I refers to the statistical units, and T refers to the time points.
* init.par: The output of the Eigen.CWM_HMM_init() function.
* tol: Threshold for ECM algorithm convergence. Default value is 0.001.
* maxit: Maximum number of iterations for the ECM algorithm. Default value is 500.
* nThreads: A positive integer indicating the number of cores used for running in parallel.

## extract.bestM.CWM ##

### Description ###

This function extracts the top/best fitting models according to the Bayesian information criterion (BIC).

### Usage ###

extract.bestM.CWM (results, top = 1)

### Arguments ###

* results: The output of the Eigen.CWM_HMM_fit() function.
* top: A number indicating how many models to extract from the ranking provided by the BIC. 

## extract.shortEM ##

### Description ###

This function extracts the names and k of the best fitting models provided by the extract.bestM function.

### Usage ###

extract.shortEM (Extract.bestM)

### Arguments ###

* Extract.bestM: The output of the extract.bestM.CWM() function.

## filt.init ##

### Description ###

This function extracts only the initializations to be used in the second step of our fitting strategy, according to the results provided by the extract.shortEM() function.

### Usage ###

filt.init (res.init, Extract.bestM)

### Arguments ###

* res.init: The output of the Eigen.CWM_HMM_init() function.
* Extract.bestM: The output of the extract.shortEM() function.

## r.CWM_HMM ##

### Description ###

This function generates random observations from a MV-HMRMRC.

### Usage ###

r.CWM_HMM (num, t, PI, M, B, U.Y, V.Y, U.X, V.X, IP, seed = NULL)

### Arguments ###

* num: The number of statistical units.
* t: The number of time points.
* PI: The transition probability matrix of dimension k x k, where k is the number of hidden states.
* M: An array containing the mean matrices for the covariates having dimensions Q x R x k, where Q and R are the rows and columns of the covariates, respectively. 
* B: An array containing the regression coefficients having dimensions P x Q+1 x k, where P and Q are the rows of the responses and the covariates, respectively. The first column must contains the interecpts.
* U.Y: An array containing the row covariance matrices of the respones having dimensions P x P x k.
* V.Y: An array containing the column covariance matrices of the respones having dimensions R x R x k.
* U.X: An array containing the row covariance matrices of the covariates having dimensions Q x Q x k.
* V.X: An array containing the column covariance matrices of the covariates having dimensions R x R x k. 
* IP: A vector containing the initial probability weights.
* seed: The seed for random number generation.

# Parsimonious MV-HMRMFCs

## Eigen.FMR_HMM_init ##

### Description ###

Runs the initialization of the ECM algorithm used for fitting the 98 parsimonious MV-HMRMFCs. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

Eigen.FMR_HMM_init (Y, X, k, mod.row.Y = "all", mod.col.Y = "all", nstartR = 50, nThreads = 1)

### Arguments ###

The meaning of the arguments of this function is the same as that of the Eigen.CWM_HMM_init() function. See it for further details.

## Eigen.FMR_HMM_fit ##

### Description ###

Fits, by using the ECM algorithm, the 98 parsimonious MV-HMRMRFs. Parallel computing is implemented and highly recommended for a faster calculation.

### Usage ###

Eigen.FMR_HMM_fit (Y, X, init.par = NULL, tol = 0.001, maxit = 500, nThreads = 1)

### Arguments ###

The meaning of the arguments of this function is the same as that of the Eigen.CWM_HMM_fit() function. See it for further details.

## extract.bestM.FMR ##

### Description ###

This function extracts the top/best fitting models according to the Bayesian information criterion (BIC).

### Usage ###

extract.bestM.CWM (results, top = 1)

### Arguments ###

The meaning of the arguments of this function is the same as that of the extract.bestM.CWM() function. See it for further details.

## r.FMR_HMM ##

### Description ###

This function generates random observations from a MV-HMRMFC.

### Usage ###

r.FMR_HMM (num, t, PI, B, U.Y, V.Y, IP, seed = NULL)

### Arguments ###

The meaning of the arguments of this function is the same as that of the r.CWM_HMM() function. See it for further details.

# Example files

The two example files allow to simulate and fit the two families of MV-HMRMs.

# Real dataset

The real dataset has been extracted from http://dati.istat.it/#.
It contains the unemployment (X) and the labor force participation (Y) for I = 106 Italian provinces over T = 4 years spanning from 2018 to 2021.
For each province, the unemployment and the labor force participation are recorded in a two-factor format.
The first factor, gender, has two levels: males and females.
The second factor, age, has four levels driven by the age class: 15-24, 25-34, 35-49, and 50-74.
Therefore, both variables are presented in a four-way array format having dimensions 2 x 4 x 106 x 4.
