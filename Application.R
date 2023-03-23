### Load the data here ###

## Fitting 

# Note that we run the analysis by parallelizing computation on 14 cores, i.e. nThreads = 14 in the code below.
# If this number of cores is not available, nThreads should be lowered.

res.init <- Eigen.CWM_HMM_init(Y=Y,X=X,k=1:7, mod.row.Y = "all", mod.col.Y = "all", mod.row.X = "all", mod.col.X = "all", nThreads = 14)
res.fit <- Eigen.CWM_HMM_fit(Y=Y,X=X,init.par = res.init, nThreads = 14, maxit=5)

n.top <- as.integer(98*98*0.01*7)

win.sht.1 <- extract.bestM.CWM(results = res.fit, top = n.top)
win.sht.2 <- extract.shortEM(Extract.bestM = win.sht.1)
win.sht.3 <- filt.init(res.init, win.sht.2)

res.fin <- win.fin <- win.tot <- list()

for (j in 1:length(win.sht.3)) {
  
  res.fin[[j]] <- Eigen.CWM_HMM_fit(Y=Y,X=X,init.par = win.sht.3[[j]], nThreads = 14)
  win.fin[[j]] <- extract.bestM.CWM(results = res.fin[[j]], top = 1)
  
}

win.temp <- numeric(length(win.sht.3))

for (j in 1:length(win.sht.3)) {
  
  win.temp[j] <- win.fin[[j]][[1]][["BIC"]]
  
}

win.tot <- win.fin[[which.min(win.temp)]][[1]]



