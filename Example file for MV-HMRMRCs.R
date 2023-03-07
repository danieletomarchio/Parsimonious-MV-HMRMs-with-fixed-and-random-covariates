#### Setting the parameters ####

p<-2
r<-5
q<-3 
k<-2
t<-5
num<-100

initialw <- rep(1/k,k)

B <- array(NA, dim = c(p,q+1,k))
B[,,1]<- matrix(c(2,1,1,-1,
                  3,1,-1,1),nrow = p, ncol = q+1, byrow = TRUE) 
B[,,2]<- matrix(c(-12,1,1,-1,
                  -13,1,-1,1),nrow = p, ncol = q+1, byrow = TRUE) 

PI     <- matrix(c(0.6, 0.4,
                   0.2, 0.8), nrow=k, ncol=k,byrow = TRUE)

UY <- array(NA,dim = c(p,p,k))
UY[,,1] <- diag(c(1.73,0.58))
UY[,,2] <- diag(c(1.63,9.80))

VY <- array(NA,dim = c(r,r,k))
VY[,,1] <- diag(c(0.74,1.48,0.74,1.11,1.11))
VY[,,2] <- diag(c(2.64,1.32,0.66,0.66,0.66))

M <- array(NA, dim = c(q,r,k))
M[,,1]<- matrix(c(-4,-3,-4,-3,-4,
                  -2,-3,-3,-2,-2,
                  -3,-3,-4,-4,-3),nrow = q, ncol = r, byrow = TRUE) 
M[,,2]<- M[,,1] + 9

UX <- array(NA,dim = c(q,q,k))
UX[,,1] <- diag(c(5.55,2.77,4.16))
UX[,,2] <- diag(c(1.00,2.00,0.50))

VX <- array(NA,dim = c(r,r,k))
VX[,,1] <- diag(c(1.06,1.77,1.06,0.71,0.71))
VX[,,2] <- diag(c(0.57,0.57,2.30,2.30,0.57))

#### Simulate the data ####

data<- r.CWM_HMM(num = num, t=t, PI=PI, M=M, B = B, U.Y = UY, V.Y = VY, U.X = UX, V.X = VX, IP=initialw)

Y = data$Y
X = data$X

#### Fit the models ####

## Example for fitting one model until convergence ## 

res.init <- Eigen.CWM_HMM_init(Y=Y,X=X,k=k, mod.row.Y = "VVI", mod.col.Y = "VI", mod.row.X = "VVI", mod.col.X = "VI", nThreads = 6)
res.fit <- Eigen.CWM_HMM_fit(Y=Y,X=X,init.par = res.init, nThreads = 1)

## Example for fitting all the models by using our fitting strategy ## 

n.top <- as.integer(98*98*0.01)

# First Part 

res.init <- Eigen.CWM_HMM_init(Y=Y,X=X,k=k, mod.row.Y = "all", mod.col.Y = "all", mod.row.X = "all", mod.col.X = "all", nThreads = 14)
res.fit.pt1 <- Eigen.CWM_HMM_fit(Y=Y,X=X,init.par = res.init, nThreads = 14, maxit = 5)
win.pt1 <- extract.bestM.CWM(results = res.fit.pt1, top = n.top)
win.pt2 <- extract.shortEM(Extract.bestM = win.pt1)
win.pt3 <- filt.init(res.init,win.pt2)

# Second Part 

res.fin<- Eigen.CWM_HMM_fit(Y=Y,X=X,init.par = win.pt3[[1]], nThreads = 14)
win.fin<- extract.bestM.CWM(results = res.fin,top = 1)





