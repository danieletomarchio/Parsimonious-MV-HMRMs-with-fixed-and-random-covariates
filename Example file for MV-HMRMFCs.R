#### Setting the parameters ####

p<-3
r<-3
q<-3 
k<-3
t<-5
num<-100

initialw <- rep(1/k,k)

B <- array(NA, dim = c(p,q+1,k))
B[,,1]<- matrix(c(0,1,1,-1,
                  1,1,-1,1,
                  0,1,1,-1),nrow = p, ncol = q+1, byrow = TRUE) 
B[,,2]<- matrix(c(-7,1,1,-1,
                  -8,1,-1,1,
                  -7,1,1,-1),nrow = p, ncol = q+1, byrow = TRUE) 
B[,,3]<- matrix(c(7,1,1,-1,
                  8,1,-1,1,
                  7,1,1,-1),nrow = p, ncol = q+1, byrow = TRUE) 

PI     <- matrix(c(0.7, 0.2, 0.1,
                   0.1, 0.8, 0.1,
                   0.1, 0.3, 0.6), nrow=k, ncol=k,byrow = TRUE)

UY <- array(NA,dim = c(p,p,k))
UY[,,1] <- diag(c(0.64,2.12,0.25))
UY[,,2] <- diag(c(1.50,2.16,2.47))
UY[,,3] <- diag(c(7.16,2.51,3.56))

VY <- array(NA,dim = c(r,r,k))
VY[,,1] <- diag(c(0.58,1.92,0.90))
VY[,,2] <- diag(c(0.46,0.73,2.98))
VY[,,3] <- diag(c(0.49,0.37,5.50))

#### Simulate the data ####

data<- r.FMR_HMM(num = num, t=t, PI=PI, B = B, U.Y = UY, V.Y = VY, IP=initialw)

Y = data$Y
X = data$X

#### Fit the models ####

## Example for fitting one model ## 

res.init <- Eigen.FMR_HMM_init(Y=Y,X=X,k=k, mod.row.Y = "VVI", mod.col.Y = "VI", nThreads = 6)
res.fit <- Eigen.FMR_HMM_fit(Y=Y,X=X,init.par = res.init, nThreads = 1)

## Example for fitting all the models ## 

res.init <- Eigen.FMR_HMM_init(Y=Y,X=X,k=k, mod.row.Y = "all", mod.col.Y = "all", nThreads = 14)
res.fit <- Eigen.FMR_HMM_fit(Y=Y,X=X,init.par = res.init, nThreads = 14)
win <- extract.bestM.FMR(results = res.fit, top = 1)
