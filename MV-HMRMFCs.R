library(foreach)
library(doSNOW)
library(progress)
library(tensor)

## The following other packages must be installed in R for the running of our code:
# withr
# LaplacesDemon
# mclust
# tidyr
# data.table

### HMM-FMR Parsimonious Matrix - Variate ###

r.FMR_HMM <- function(num, t, PI, B, U.Y, V.Y, IP, seed = NULL) {
  run.mc.sim <- function(P, times, start = 1, IP) {
    
    # number of possible states
    num.states <- nrow(P)
    
    # stores the states X_t through time
    states <- numeric(times)
    
    # initialize variable for first state
    
    k <- length(IP)
    
    states[1] <- sample(1:k, size = 1, prob = IP)
    
    # states[1]    <- start
    
    for (t in 2:times) {
      
      # probability vector to simulate next state X_{t+1}
      p <- P[states[t - 1], ]
      
      ## draw from multinomial and determine state
      states[t] <- which(rmultinom(1, 1, p) == 1)
    }
    return(states)
  }
  
  row.Y <- dim(U.Y)[1]
  col.Y <- col.X <- dim(V.Y)[2]
  row.X <- dim(B)[2]-1
  k <- nrow(PI)
  
  Y <- array(NA, dim = c(row.Y, col.Y, num, t))
  X <- array(NA, dim = c(row.X, col.X, num, t))
  
  # each column stores the sequence of states for a single chains
  chain.states <- matrix(NA, ncol = num, nrow = t)
  
  # simulate chains
  for (c in seq_len(num)) {
    chain.states[, c] <- run.mc.sim(PI, t, start = start, IP = IP)
  }
  
  if (is.null(seed)) {
    rm(.Random.seed, envir = globalenv())
    seed <- round(runif(1, 1, .Machine$integer.max))
  } else {
    seed <- seed
  }
  scal=1
  withr::with_seed(seed, for (j in seq_len(t)) {
    for (i in 1:num) {
      X[, , i, j] <- LaplacesDemon::rmatrixnorm(M = matrix(c(runif(row.X*col.X, j*scal, (j*scal)+1)),row.X,col.X), U = diag(row.X), V = diag(col.X))
      X1 <- rbind(rep(1, r), X[, , i, j])
      Y[, , i, j] <- B[, , chain.states[j, i]] %*% X1 + LaplacesDemon::rmatrixnorm(M = matrix(0, p, r), U = U.Y[, , chain.states[j, i]], V = V.Y[, , chain.states[j, i]])
    }
  })
  
  chain.states <- t(chain.states)
  colnames(chain.states) <- paste("time", 1:t, sep = " ")
  
  return(list(Y = Y, X = X, obs.states = chain.states))
} # random generation

Eigen.FMR_HMM_init <- function(Y, X, k, mod.row.Y = "all", mod.col.Y = "all", nstartR = 50, nThreads = 1) {
  r_Pars_init <- function(Y, X, k, mod.row.Y, mod.col.Y, nstartR = 50) {
    dMVnorm <- function(X, M, U, V) {
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X
      
      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }
      
      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)
      
      return(pdf)
    }
    tr <- function(x) {
      return(sum(diag(x)))
    }
    transit <- function(Y) {
      
      # Y is a matrix N \times T
      
      N <- dim(Y)[1]
      T <- dim(Y)[2]
      K <- max(Y)
      
      if (K == 1) {
        PI <- matrix(1, 1, 1)
      }
      
      if (K > 1) {
        PI <- matrix(0, nrow = K, ncol = K)
        
        for (i in 1:N) {
          for (t in 2:T) {
            PI[Y[i, t - 1], Y[i, t]] <- PI[Y[i, t - 1], Y[i, t]] + 1
          }
        }
        PI <- diag(1 / rowSums(PI)) %*% PI
      }
      
      return(PI)
    }
    
    # Dimensions
    
    num0 <- dim(X)[3] 
    p <- nrow(Y) 
    r <- ncol(Y) 
    q <- nrow(X) + 1
    t <- dim(X)[4] 
    
    num <- num0 * t
    
    Xresh <- array(X, dim = c(q - 1, r, num)) 
    Yresh <- array(Y, dim = c(p, r, num)) 
    
    # Add intercept to X
    
    unos <- rep(1, r)
    X1resh <- array(0, dim = c(q, r, num))
    for (i in 1:num) {
      X1resh[, , i] <- rbind(unos, Xresh[, , i])
    }
    
    # Create some objects
    
    prior <- matrix(NA, nstartR, k)
    llk <- rep(NA, nstartR)
    
    coeff <- array(0, dim = c(p, q, k, nstartR))
    
    WRY <- array(0, dim = c(p, p, k, nstartR))
    WCY <- array(0, dim = c(r, r, k, nstartR))
    
    sigmaUY <- array(0, dim = c(p, p, k, nstartR))
    sigmaVY <- array(0, dim = c(r, r, k, nstartR))
    
    sel <- array(0, dim = c(p, r, k))
    dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
    post2 <- array(NA, c(k, k, num0, t - 1)) # transition probabilities
    f <- array(NA, c(num0, t, k, nstartR))
    A <- B <- array(NA, c(num0, t, k, nstartR)) 
    PI <- array(NA, dim = c(k, k, nstartR))
    
    ## Random initialization ##
    
    eu <- matrix(0, nrow = num, ncol = k)
    classy <- numeric(num)
    rand.start <- matrix(0, nstartR, k)
    
    withr::with_seed(p, for (i in 1:nstartR) {
      rand.start[i, ] <- sample(c(1:num), k)
    })
    
    for (h in 1:nstartR) {
      skip_to_next <- FALSE
      
      ### part 0 ###
      
      tryCatch(
        {
          sec <- rand.start[h, ]
          
          for (j in 1:k) {
            sel[, , j] <- Yresh[, , sec[j]]
          }
          
          for (j in 1:k) {
            for (i in 1:(num)) {
              eu[i, j] <- norm((Yresh[, , i] - sel[, , j]), type = "F")
            }
          }
          
          for (i in 1:(num)) {
            classy[i] <- which.min(eu[i, ])
          }
          
          z <- mclust::unmap(classy)
          
          ### part 1 ###
          
          # Regression Coefficients + Partial U & V #
          
          for (j in 1:k) {
            
            # Y #
            
            tempcoeff1 <- tensor::tensor(aperm(tensor::tensor(Yresh * z[, j][slice.index(Yresh, 3)], diag(r), 2, 1), c(1, 3, 2)), aperm(X1resh, c(2, 1, 3)), c(2, 3), c(1, 3))
            tempcoeff2 <- tensor::tensor(aperm(tensor::tensor(X1resh * z[, j][slice.index(X1resh, 3)], diag(r), 2, 1), c(1, 3, 2)), aperm(X1resh, c(2, 1, 3)), c(2, 3), c(1, 3))
            
            coeff[, , j, h] <- tempcoeff1 %*% solve(tempcoeff2)
            
            BX <- tensor::tensor(coeff[, , j, h], X1resh, 2, 1)
            
            WRY[, , j, h] <- tensor::tensor(aperm(tensor::tensor((Yresh - BX) * z[, j][slice.index(Yresh - BX, 3)], diag(r), 2, 1), c(1, 3, 2)), aperm((Yresh - BX), c(2, 1, 3)), c(2, 3), c(1, 3))
            WCY[, , j, h] <- tensor::tensor(aperm(tensor::tensor(aperm(Yresh - BX, c(2, 1, 3)) * z[, j][slice.index(aperm(Yresh - BX, c(2, 1, 3)), 3)], solve(WRY[, , j, h]), 2, 1), c(1, 3, 2)), (Yresh - BX), c(2, 3), c(1, 3))
            
          }
          
          # Row covariance matrix Y #
          
          if (mod.row.Y == "EII") {
            if (k == 1) {
              phiY <- tr(WRY[, , , h]) / ((num) * p * r)
            } else {
              phiY <- tr(rowSums(WRY[, , , h], dims = 2)) / ((num) * p * r)
            }
            
            for (j in 1:k) {
              sigmaUY[, , j, h] <- phiY * diag(1, p, p)
            }
          }
          
          if (mod.row.Y == "EEI") {
            if (k == 1) {
              deltaU <- diag(diag(WRY[, , , h]), p, p) / (det(diag(diag(WRY[, , , h]), p, p)))^(1 / p)
              
              phiY <- (det(diag(diag(WRY[, , , h]), p, p)))^(1 / p) / ((num) * r)
            } else {
              deltaU <- diag(diag(rowSums(WRY[, , , h], dims = 2)), p, p) / (det(diag(diag(rowSums(WRY[, , , h], dims = 2)), p, p)))^(1 / p)
              
              phiY <- (det(diag(diag(rowSums(WRY[, , , h], dims = 2)), p, p)))^(1 / p) / ((num) * r)
            }
            
            for (j in 1:k) {
              sigmaUY[, , j, h] <- phiY * deltaU
            }
          }
          
          if (mod.row.Y == "EEE") {
            if (k == 1) {
              sigmaUY[, , 1, h] <- WRY[, , , h] / (num * r)
            } else {
              for (j in 1:k) {
                sigmaUY[, , j, h] <- rowSums(WRY[, , , h], dims = 2) / (num * r)
              }
            }
          }
          
          # Column covariance matrix Y #
          
          if (mod.col.Y == "II") {
            for (j in 1:k) {
              sigmaVY[, , j, h] <- diag(1, r, r)
            }
          }
          
          if (mod.col.Y == "EI") {
            if (k == 1) {
              deltaVY <- diag(diag(WCY[, , , h]), r, r) / (det(diag(diag(WCY[, , , h]), r, r)))^(1 / r)
            } else {
              deltaVY <- diag(diag(rowSums(WCY[, , , h], dims = 2)), r, r) / (det(diag(diag(rowSums(WCY[, , , h], dims = 2)), r, r)))^(1 / r)
            }
            
            for (j in 1:k) {
              sigmaVY[, , j, h] <- deltaVY
            }
          }
          
          if (mod.col.Y == "EE") {
            if (k == 1) {
              sigmaVY[, , 1, h] <- WCY[, , , h] / (det(WCY[, , , h]))^(1 / r)
            } else {
              for (j in 1:k) {
                sigmaVY[, , j, h] <- rowSums(WCY[, , , h], dims = 2) / ((det(rowSums(WCY[, , , h], dims = 2)))^(1 / r))
              }
            }
          }
          
          # Density computation #
          
          for (j in 1:k) {
            for (i in 1:num) {
              dens[i, j] <- dMVnorm(X = Yresh[, , i], M = coeff[, , j, h] %*% X1resh[, , i], U = sigmaUY[, , j, h], V = sigmaVY[, , j, h])
            }
          }
          
          if (k == 1) {
            prior[h, ] <- 1
          } else {
            prior[h, ] <- colMeans(z)
          }
          
          ### part 2 ###
          
          f[, , , h] <- array(dens, c(num0, t, k))
          
          PI[, , h] <- transit(matrix(classy, num0, t))
          
          A[, 1, , h] <- matrix(rep(prior[h, ], each = num0), ncol = k) * f[, 1, , h] 
          B[, t, , h] <- 1
          if(k==1){
            
            for (T in 2:t) {
              A[, T, , h] <- A[, T - 1, , h] %*% as.matrix(PI[, , h],1,1) * f[, T, , h]
              tp <- t - T + 1
              B[, tp, , h] <- t(as.matrix(PI[, , h],1,1) %*% t(B[, tp + 1, , h] * f[, tp + 1, , h]))
            }
            
          }else{
            
            for (T in 2:t) {
              A[, T, , h] <- A[, T - 1, , h] %*% PI[, , h] * f[, T, , h]
              tp <- t - T + 1
              B[, tp, , h] <- t(PI[, , h] %*% t(B[, tp + 1, , h] * f[, tp + 1, , h]))
            }
            
          }
          
          part <- rowSums(matrix(A[, t, , h], nrow = num0, ncol = k))
          part[part == 0] <- .Machine$double.xmin
          llk[h] <- sum(log(part)) # log-likelihood
          
          if(any(prior[h, ])<=0.05)
            llk[h] <- NA
        },
        error = function(e) {
          skip_to_next <<- TRUE
        }
      )
      
      if (skip_to_next) {
        next
      }
    }
    
    df <- data.frame(llk = llk, pos = c(1:nstartR))
    df <- tidyr::drop_na(df)
    df <- df[!is.infinite(rowSums(df)), ]
    bestR <- head(data.table::setorderv(df, cols = "llk", order = -1), n = 1)$pos
    
    CF <- array(coeff[, , , bestR], dim = c(p, q, k))
    WR <- array(sigmaUY[, , , bestR], dim = c(p, p, k))
    WC <- array(sigmaVY[, , , bestR], dim = c(r, r, k))
    PI2 <- matrix(PI[, , bestR], k, k)
    ef <- array(f[, , , bestR], c(num0, t, k))
    a <- array(A[, , , bestR], c(num0, t, k))
    b <- array(B[, , , bestR], c(num0, t, k))
    lik <- llk[bestR]
    
    return(list(
      model = c(mod.row.Y, mod.col.Y), prior = prior[bestR, ],
      coeff = CF, sigmaUY = WR, sigmaVY = WC,
      PI = PI2, f = ef, A = a, B = b, llk = lik
    ))
  }
  
  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }
  
  if (any(mod.row.Y == "all")) {
    mod.row.Y.i <- c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "VEE", "EVE", "VVE", "EEV", "VEV", "EVV", "VVV")
  } else {
    mod.row.Y.i <- mod.row.Y
  }
  if (any(mod.col.Y == "all")) {
    mod.col.Y.i <- c("II", "EI", "VI", "EE", "VE", "EV", "VV")
  } else {
    mod.col.Y.i <- mod.col.Y
  }
  
  req.model <- expand.grid(mod.row.Y.i, mod.col.Y.i)
  names(req.model) <- paste(rep("V", each =  ncol(req.model)), 1:ncol(req.model), sep = "")
  
  nest.EII.Y <- c("EII", "VII")
  nest.EEI.Y <- c("EEI", "VEI", "EVI", "VVI")
  nest.EEE.Y <- c("EEE", "VEE", "EVE", "EEV", "VVE", "VEV", "EVV", "VVV")
  nest.II.Y <- c("II")
  nest.EI.Y <- c("EI", "VI")
  nest.EE.Y <- c("EE", "VE", "EV", "VV")
  
  list.comb <- data.frame(matrix(NA, nrow = nrow(req.model), ncol = 2))
  for (i in 1:nrow(req.model)) {
    if (req.model[i, 1] %in% nest.EII.Y) {
      list.comb[i, 1] <- "EII"
    }
    if (req.model[i, 1] %in% nest.EEI.Y) {
      list.comb[i, 1] <- "EEI"
    }
    if (req.model[i, 1] %in% nest.EEE.Y) {
      list.comb[i, 1] <- "EEE"
    }
    
    if (req.model[i, 2] %in% nest.II.Y) {
      list.comb[i, 2] <- "II"
    }
    if (req.model[i, 2] %in% nest.EI.Y) {
      list.comb[i, 2] <- "EI"
    }
    if (req.model[i, 2] %in% nest.EE.Y) {
      list.comb[i, 2] <- "EE"
    }
    
  }
  
  list.comb2 <- unique(list.comb)
  
  oper <- vector(mode = "list", length = length(k))
  
  for (g in 1:length(k)) {
    print(paste("Initializing Parsimonious MV-HMRMFCs with k =", k[g]))
    
    cluster <- makeCluster(nThreads, type = "SOCK")
    registerDoSNOW(cluster)
    
    pb <- progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = nrow(list.comb2),
      complete = "=", # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">", # Current bar character
      width = 100
    )
    progress <- function(n) {
      pb$tick()
    }
    opts <- list(progress = progress)
    
    oper[[g]] <- foreach(l = 1:nrow(list.comb2), .combine = "comb", .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- r_Pars_init(
        Y = Y, X = X, k = k[g], mod.row.Y = list.comb2[l, 1], mod.col.Y = list.comb2[l, 2],nstartR = nstartR)
      
      list(res)
    }
    
    if (g == 1) {
      oper2 <- foreach(i = 1:nrow(req.model), .combine = "comb", .multicombine = TRUE, .init = list(list())) %dopar% {
        for (j in 1:nrow(list.comb2)) {
          if (all(list.comb[i, ] == oper[[1]][[1]][[j]][["model"]])) {
            res <- j
          }
        }
        
        list(res)
      }
    }
    
    stopCluster(cluster)
    registerDoSEQ()
  }
  
  return(list(
    results = oper,
    k = k,
    req.model = req.model,
    init.used = list.comb,
    index = unlist(oper2[[1]])
  ))
} # initializing function

Eigen.FMR_HMM_fit <- function(Y, X, init.par = NULL, tol = 0.001, maxit = 500, nThreads = 1) {
  k <- init.par[[2]]
  list.comb <- init.par[[3]]
  pt.mod <- nrow(list.comb)
  
  tol2 <- 0.001
  maxit2 <- 100
  
  comb <- function(x, ...) {
    lapply(
      seq_along(x),
      function(i) c(x[[i]], lapply(list(...), function(y) y[[i]]))
    )
  }
  oper <- vector(mode = "list", length = length(k))
  time <- numeric(length(k))
  
  Pars_MVN_FMR <- function(Y, X, k, init.par = NULL, mod.row.Y = NULL, mod.col.Y = NULL, tol = 0.001, tol2 = 0.001, maxit = 500, maxit2 = 100) {
    ptm <- proc.time()
    
    # Functions
    
    dMVnorm <- function(X, M, U, V) {
      tr <- function(x) {
        return(sum(diag(x)))
      }
      
      num <- dim(X)[3] # sample size
      p <- nrow(X) # rows of X
      r <- ncol(X) # columns of X
      
      if (is.na(num)) {
        X <- as.matrix(X)
        delta <- tr(solve(U) %*% (X - M) %*% solve(V) %*% t(X - M))
      } else {
        delta <- sapply(1:num, function(i) tr(solve(U) %*% (X[, , i] - M) %*% solve(V) %*% t(X[, , i] - M)))
      }
      
      pdf <- (2 * pi)^(-(p * r) / 2) * det(U)^(-r / 2) * det(V)^(-p / 2) * exp(-1 / 2 * delta)
      
      return(pdf)
    }
    tr <- function(x) {
      return(sum(diag(x)))
    }
    
    # Dimensions
    
    num0 <- dim(X)[3] 
    p <- nrow(Y) 
    r <- ncol(Y) 
    q <- nrow(X) + 1 
    t <- dim(X)[4] 
    
    num <- num0 * t
    
    Xresh <- array(X, dim = c(q - 1, r, num)) 
    Yresh <- array(Y, dim = c(p, r, num)) 
    
    # Add intercept to X
    
    unos <- rep(1, r)
    X1resh <- array(0, dim = c(q, r, num))
    for (i in 1:num) {
      X1resh[, , i] <- rbind(unos, Xresh[, , i])
    }
    
    ## Objects related to the two covariance matrices
    
    WRY <- array(0, dim = c(p, p, k))
    WCY <- array(0, dim = c(r, r, k))
    
    phiY.k <- phiX.k <- numeric(k)
    deltaUY.k <- gammaUY.k <- array(0, dim = c(p, p, k))
    deltaVY.k <- gammaVY.k <- array(0, dim = c(r, r, k))
    
    temp.numdeltaR <- ftemp.r <- tempomegaR <- V_EVV.UY.K <- array(0, dim = c(p, p, k)) # for VEI, VEE, VEV, MM object
    temp.phi2 <- numeric(k) # for EVI, EVE, EVV
    tempW_EEV <- tempW_EV <- vector("list", k) # for EEV, VEV, EV
    ftemp.c <- tempomegaC <- array(0, dim = c(r, r, k)) # MM object
    
    ## Other objects
    
    post <- dens <- array(0, c(num, k), dimnames = list(1:(num), paste("comp.", 1:k, sep = "")))
    post2 <- array(NA, c(k, k, num0, t - 1)) # transition probabilities
    
    # Preliminary definition of convergence criteria for EM/MM algorithms
    
    check <- 0
    check2 <- 0
    loglik.old <- -Inf
    loglik.new <- NULL
    ll <- NULL
    mark <- 1  
    MM.r.old <- -Inf
    MM.c.old <- -Inf
    m.iter <- 0
    m.iter2 <- 0
    
    ### Algorithm ###
    
    coeff <- init.par$coeff
    sigmaUY <- init.par$sigmaUY
    sigmaVY <- init.par$sigmaVY
    prior <- init.par$prior
    f <- init.par$f
    A <- init.par$A
    B <- init.par$B
    PI <- init.par$PI
    
    ei.UY <- ei.VY <- vector(mode = "list", length = k)
    
    for (j in 1:k) {
      ei.UY[[j]] <- eigen(sigmaUY[, , j])
      ei.VY[[j]] <- eigen(sigmaVY[, , j])
      
      phiY.k[j] <- prod(ei.UY[[j]][["values"]])^(1 / p)
      deltaUY.k[, , j] <- diag(ei.UY[[j]][["values"]] / phiY.k[j])
      gammaUY.k[, , j] <- ei.UY[[j]][["vectors"]]
      
      deltaVY.k[, , j] <- diag(ei.VY[[j]][["values"]])
      gammaVY.k[, , j] <- ei.VY[[j]][["vectors"]]
      
    }
    
    if (mod.row.Y %in% c("EVE", "VVE")) {
      gammaU <- rowSums(gammaUY.k, dims = 2) / k
    }
    if (mod.col.Y == "VE") {
      gammaV <- rowSums(gammaVY.k, dims = 2) / k
    }
    
    ### Estimation ###
    
    while (check < 1) {
      m.iter <- m.iter + 1
      
      ### E - STEP ###
      
      for(j in 1:k){
        
        numer <- (A[,,j]*B[,,j])
        numer[numer<5e-324] <- 5e-324
        denom <- apply(A*B,c(1,2),sum)
        denom[denom<5e-324] <- 5e-324
        
        post[,j] <- numer/denom # posterior probabilities (z)
        
        for(n in 1:num0){
          for(T in 1:(t-1)){
            
            post2[,,n,T] <- diag(A[n,T,],nrow=k,ncol=k) %*% PI %*% diag(B[n,T+1,]*f[n,T+1,],nrow=k,ncol=k)
            post2[,,n,T] <- post2[,,n,T]/sum(post2[,,n,T]) # posterior probabilities (zz)
            post2[,,n,T][is.nan(post2[,,n,T])] <- 0
            
          }
        }
        
      }  
      
      ### M - STEP ###
      
      for (j in 1:k) {
        
        tempcoeff1 <- tensor::tensor(aperm(tensor::tensor(Yresh * post[, j][slice.index(Yresh, 3)], solve(sigmaVY[, , j]), 2, 1), c(1, 3, 2)), aperm(X1resh, c(2, 1, 3)), c(2, 3), c(1, 3))
        tempcoeff2 <- tensor::tensor(aperm(tensor::tensor(X1resh * post[, j][slice.index(X1resh, 3)], solve(sigmaVY[, , j]), 2, 1), c(1, 3, 2)), aperm(X1resh, c(2, 1, 3)), c(2, 3), c(1, 3))
        
        coeff[, , j] <- tempcoeff1 %*% solve(tempcoeff2)
        
      }
      
      # ROWS COVARIANCE MATRIX Y #
      
      for (j in 1:k) {
        
        BX <- tensor::tensor(coeff[, , j], X1resh, 2, 1)
        WRY[, , j] <- tensor::tensor(aperm(tensor::tensor((Yresh - BX) * post[, j][slice.index(Yresh - BX, 3)], solve(sigmaVY[, , j]), 2, 1), c(1, 3, 2)), aperm((Yresh - BX), c(2, 1, 3)), c(2, 3), c(1, 3))
        
      }
      
      if (mod.row.Y == "EII") {
        phiY <- tr(rowSums(WRY, dims = 2)) / ((num) * p * r)
        
        for (j in 1:k) {
          sigmaUY[, , j] <- phiY * diag(1, p, p)
        }
      }
      
      if (mod.row.Y == "VII") {
        for (j in 1:k) {
          phiY.k[j] <- tr(WRY[, , j]) / (p * r * sum(post[, j]))
          sigmaUY[, , j] <- phiY.k[j] * diag(1, p, p)
        }
      }
      
      if (mod.row.Y == "EEI") {
        deltaU <- diag(diag(rowSums(WRY, dims = 2)), p, p) / (det(diag(diag(rowSums(WRY, dims = 2)), p, p)))^(1 / p)
        
        phiY <- (det(diag(diag(rowSums(WRY, dims = 2)), p, p)))^(1 / p) / ((num) * r)
        
        for (j in 1:k) {
          sigmaUY[, , j] <- phiY * deltaU
        }
      }
      
      if (mod.row.Y == "VEI") {
        for (j in 1:k) {
          temp.numdeltaR[, , j] <- (1 / phiY.k[j]) * WRY[, , j]
        }
        
        deltaU <- diag(diag(rowSums(temp.numdeltaR, dims = 2)), p, p) / (det(diag(diag(rowSums(temp.numdeltaR, dims = 2)), p, p)))^(1 / p)
        
        for (j in 1:k) {
          phiY.k[j] <- (tr(solve(deltaU) %*% WRY[, , j])) / (p * r * sum(post[, j]))
          
          sigmaUY[, , j] <- phiY.k[j] * deltaU
        }
      }
      
      if (mod.row.Y == "EVI") {
        for (j in 1:k) {
          deltaUY.k[, , j] <- diag(diag(WRY[, , j]), p, p) / (det(diag(diag(WRY[, , j]), p, p)))^(1 / p)
          
          temp.phi2[j] <- det(diag(diag(WRY[, , j]), p, p))^(1 / p)
        }
        
        phiY <- sum(temp.phi2) / ((num) * r)
        
        for (j in 1:k) {
          sigmaUY[, , j] <- phiY * deltaUY.k[, , j]
        }
      }
      
      if (mod.row.Y == "VVI") {
        for (j in 1:k) {
          deltaUY.k[, , j] <- diag(diag(WRY[, , j]), p, p) / (det(diag(diag(WRY[, , j]), p, p)))^(1 / p)
          
          phiY.k[j] <- det(diag(diag(WRY[, , j]), p, p))^(1 / p) / (r * sum(post[, j]))
          
          sigmaUY[, , j] <- phiY.k[j] * deltaUY.k[, , j]
        }
      }
      
      if (mod.row.Y == "EEE") {
        for (j in 1:k) {
          sigmaUY[, , j] <- rowSums(WRY, dims = 2) / (num * r)
        }
      }
      
      if (mod.row.Y == "VEE") {
        for (j in 1:k) {
          temp.numdeltaR[, , j] <- (1 / phiY.k[j]) * WRY[, , j]
        }
        
        deltaU <- rowSums(temp.numdeltaR, dims = 2) / ((det(rowSums(temp.numdeltaR, dims = 2)))^(1 / p))
        
        for (j in 1:k) {
          phiY.k[j] <- tr(solve(deltaU) %*% WRY[, , j]) / (p * r * sum(post[, j]))
          
          sigmaUY[, , j] <- phiY.k[j] * deltaU
        }
      }
      
      if (mod.row.Y == "EVE") {
        while (check2 < 1) {
          m.iter2 <- m.iter2 + 1
          
          for (j in 1:k) {
            ftemp.r[, , j] <- tcrossprod(solve(deltaUY.k[, , j]), gammaU) %*% WRY[, , j] - max(eigen(WRY[, , j])$values) * tcrossprod(solve(deltaUY.k[, , j]), gammaU)
          }
          
          f <- rowSums(ftemp.r, dims = 2)
          
          MM.r.new <- tr(f %*% gammaU)
          
          if ((abs(MM.r.new - MM.r.old)) < tol2 | m.iter2 == maxit2) {
            check2 <- 1
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          } else {
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          }
          
          MM.r.old <- MM.r.new
        }
        
        m.iter2 <- 0
        check2 <- 0
        MM.r.old <- -Inf
        
        for (j in 1:k) {
          deltaUY.k[, , j] <- diag(diag(crossprod(gammaU, WRY[, , j]) %*% gammaU), p, p) / (det(diag(diag(crossprod(gammaU, WRY[, , j]) %*% gammaU), p, p)))^(1 / p)
          
          temp.phi2[j] <- tr(gammaU %*% tcrossprod(solve(deltaUY.k[, , j]), gammaU) %*% WRY[, , j])
        }
        
        phiY <- sum(temp.phi2) / ((num) * p * r)
        
        for (j in 1:k) {
          sigmaUY[, , j] <- phiY * gammaU %*% tcrossprod(deltaUY.k[, , j], gammaU)
        }
      }
      
      if (mod.row.Y == "VVE") {
        while (check2 < 1) {
          m.iter2 <- m.iter2 + 1
          
          for (j in 1:k) {
            ftemp.r[, , j] <- tcrossprod(solve(deltaUY.k[, , j]), gammaU) %*% WRY[, , j] - max(eigen(WRY[, , j])$values) * tcrossprod(solve(deltaUY.k[, , j]), gammaU)
          }
          
          f <- rowSums(ftemp.r, dims = 2)
          
          MM.r.new <- tr(f %*% gammaU)
          
          if ((abs(MM.r.new - MM.r.old)) < tol2 | m.iter2 == maxit2) {
            check2 <- 1
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          } else {
            res.svd <- svd(f)
            gammaU <- tcrossprod(res.svd$v, res.svd$u)
          }
          
          MM.r.old <- MM.r.new
        }
        
        m.iter2 <- 0
        check2 <- 0
        MM.r.old <- -Inf
        
        for (j in 1:k) {
          deltaUY.k[, , j] <- diag(diag(crossprod(gammaU, WRY[, , j]) %*% gammaU), p, p) / (det(diag(diag(crossprod(gammaU, WRY[, , j]) %*% gammaU), p, p)))^(1 / p)
          phiY.k[j] <- (det(diag(diag(crossprod(gammaU, WRY[, , j]) %*% gammaU), p, p))^(1 / p)) / (r * sum(post[, j]))
          sigmaUY[, , j] <- phiY.k[j] * gammaU %*% tcrossprod(deltaUY.k[, , j], gammaU)
        }
      }
      
      if (mod.row.Y == "EEV") {
        for (j in 1:k) {
          tempW_EEV[[j]] <- eigen(WRY[, , j])
          
          gammaUY.k[, , j] <- tempW_EEV[[j]][["vectors"]]
          
          tempomegaR[, , j] <- diag(tempW_EEV[[j]][["values"]], p, p)
        }
        
        deltaU <- rowSums(tempomegaR, dims = 2) / ((det(rowSums(tempomegaR, dims = 2)))^(1 / p))
        
        phiY <- ((det(rowSums(tempomegaR, dims = 2)))^(1 / p)) / ((num) * r)
        
        for (j in 1:k) {
          sigmaUY[, , j] <- phiY * gammaUY.k[, , j] %*% tcrossprod(deltaU, gammaUY.k[, , j])
        }
      }
      
      if (mod.row.Y == "VEV") {
        for (j in 1:k) {
          tempW_EEV[[j]] <- eigen(WRY[, , j])
          
          gammaUY.k[, , j] <- tempW_EEV[[j]][["vectors"]]
          
          tempomegaR[, , j] <- diag(tempW_EEV[[j]][["values"]], p, p)
          
          temp.numdeltaR[, , j] <- (1 / phiY.k[j]) * tempomegaR[, , j]
        }
        
        deltaU <- rowSums(temp.numdeltaR, dims = 2) / ((det(rowSums(temp.numdeltaR, dims = 2)))^(1 / p))
        
        for (j in 1:k) {
          phiY.k[j] <- tr(tempomegaR[, , j] %*% solve(deltaU)) / (p * r * sum(post[, j]))
          
          sigmaUY[, , j] <- phiY.k[j] * gammaUY.k[, , j] %*% tcrossprod(deltaU, gammaUY.k[, , j])
        }
      }
      
      if (mod.row.Y == "EVV") {
        for (j in 1:k) {
          V_EVV.UY.K[, , j] <- WRY[, , j] / ((det(WRY[, , j]))^(1 / p))
          
          temp.phi2[j] <- det(WRY[, , j])^(1 / p)
        }
        
        phiY <- sum(temp.phi2) / ((num) * r)
        
        for (j in 1:k) {
          sigmaUY[, , j] <- phiY * V_EVV.UY.K[, , j]
        }
      }
      
      if (mod.row.Y == "VVV") {
        for (j in 1:k) {
          sigmaUY[, , j] <- WRY[, , j] / (r * sum(post[, j]))
        }
      }
      
      # COLUMNS COVARIANCE MATRIX Y #
      
      for (j in 1:k) {
        
        BX <- tensor::tensor(coeff[, , j], X1resh, 2, 1)
        WCY[, , j] <- tensor::tensor(aperm(tensor::tensor(aperm(Yresh - BX, c(2, 1, 3)) * post[, j][slice.index(aperm(Yresh - BX, c(2, 1, 3)), 3)], solve(sigmaUY[, , j]), 2, 1), c(1, 3, 2)), (Yresh - BX), c(2, 3), c(1, 3))
        
      }
      
      if (mod.col.Y == "II") {
        for (j in 1:k) {
          sigmaVY[, , j] <- diag(1, r, r)
        }
      }
      
      if (mod.col.Y == "EI") {
        deltaVY <- diag(diag(rowSums(WCY, dims = 2)), r, r) / (det(diag(diag(rowSums(WCY, dims = 2)), r, r)))^(1 / r)
        
        for (j in 1:k) {
          sigmaVY[, , j] <- deltaVY
        }
      }
      
      if (mod.col.Y == "VI") {
        for (j in 1:k) {
          sigmaVY[, , j] <- diag(diag(WCY[, , j]), r, r) / (det(diag(diag(WCY[, , j]), r, r)))^(1 / r)
        }
      }
      
      if (mod.col.Y == "EE") {
        for (j in 1:k) {
          sigmaVY[, , j] <- rowSums(WCY, dims = 2) / ((det(rowSums(WCY, dims = 2)))^(1 / r))
        }
      }
      
      if (mod.col.Y == "VE") {
        while (check2 < 1) {
          m.iter2 <- m.iter2 + 1
          
          for (j in 1:k) {
            ftemp.c[, , j] <- tcrossprod(solve(deltaVY.k[, , j]), gammaV) %*% WCY[, , j] - max(eigen(WCY[, , j])$values) * tcrossprod(solve(deltaVY.k[, , j]), gammaV)
          }
          
          f.C <- rowSums(ftemp.c, dims = 2)
          
          MM.c.new <- tr(f.C %*% gammaV)
          
          if ((abs(MM.c.new - MM.c.old)) < tol2 | m.iter2 == maxit2) {
            check2 <- 1
            res.svd.C <- svd(f.C)
            gammaV <- tcrossprod(res.svd.C$v, res.svd.C$u)
          } else {
            res.svd.C <- svd(f.C)
            gammaV <- tcrossprod(res.svd.C$v, res.svd.C$u)
          }
          
          MM.c.old <- MM.c.new
        }
        
        m.iter2 <- 0
        check2 <- 0
        MM.c.old <- -Inf
        
        for (j in 1:k) {
          deltaVY.k[, , j] <- diag(diag(crossprod(gammaV, WCY[, , j]) %*% gammaV), r, r) / (det(diag(diag(crossprod(gammaV, WCY[, , j]) %*% gammaV), r, r)))^(1 / r)
        }
        
        for (j in 1:k) {
          sigmaVY[, , j] <- gammaV %*% tcrossprod(deltaVY.k[, , j], gammaV)
        }
      }
      
      if (mod.col.Y == "EV") {
        for (j in 1:k) {
          tempW_EV[[j]] <- eigen(WCY[, , j])
          
          gammaVY.k[, , j] <- tempW_EV[[j]][["vectors"]]
          
          tempomegaC[, , j] <- diag(tempW_EV[[j]][["values"]], r, r)
        }
        
        deltaVY <- rowSums(tempomegaC, dims = 2) / ((det(rowSums(tempomegaC, dims = 2)))^(1 / r))
        
        for (j in 1:k) {
          sigmaVY[, , j] <- gammaVY.k[, , j] %*% tcrossprod(deltaVY, gammaVY.k[, , j])
        }
      }
      
      if (mod.col.Y == "VV") {
        for (j in 1:k) {
          sigmaVY[, , j] <- WCY[, , j] / ((det(WCY[, , j]))^(1 / r))
        }
      }
      
      prior  <- colMeans(matrix(post[1:num0,],nrow=num0,ncol=k))     # initial prob.
      post2a <- as.matrix(apply(post2,c(1,2),sum))
      PI     <- diag(1/rowSums(post2a),nrow=k,ncol=k) %*% post2a    # Transition probability matrix
      
      for (j in 1:k) {
        for (i in 1:num) {
          dens[i, j] <- dMVnorm(X = Yresh[, , i], M = coeff[, , j] %*% X1resh[, , i], U = sigmaUY[, , j], V = sigmaVY[, , j]) 
        }
      }
      
      f      <- array(dens,c(num0,t,k))
      A[,1,] <- matrix(rep(prior,each=num0),ncol=k)*f[,1,] 
      B[,t,] <- 1
      for(T in 2:t){
        
        A[,T,]  <- A[,T-1,] %*% PI * f[,T,]
        tp      <- t-T+1
        B[,tp,] <- t(PI%*%t(B[,tp+1,]*f[,tp+1,]))
        
      }
      
      part <- rowSums(matrix(A[,t,], nrow=num0, ncol=k))
      part[part<5e-324] <- 5e-324
      
      loglik.new  <- sum(log(part)) # log-likelihood
      ll <- c(ll,loglik.new) 
      
      # stopping rule
      
      if ((loglik.new - loglik.old) < tol) {
        check <- 1
      }
      
      if (m.iter > 1) {
        if (loglik.new < loglik.old) {
          mark <- 1  
        } else {
          mark <- 0
        }  
      }
      
      if (m.iter == maxit) {
        check <- 1
      }
      
      loglik.old <- loglik.new
    }
    
    #### Output ####
    
    # --------------------- #
    # Classification Matrix #
    # --------------------- #
    
    group  <- apply(post,1,which.max)
    groupT <- array(group,c(num0,t),dimnames=list(1:num0,paste("time",1:t,sep=" ")))
    
    # -------------------- #
    # Information criteria #
    # -------------------- #
    
    # Number of parameters
    
    coeff.par <- (p * q) * k
    
    if (mod.row.Y == "EII") {
      rowparY <- 1
    }
    if (mod.row.Y == "VII") {
      rowparY <- k
    }
    if (mod.row.Y == "EEI") {
      rowparY <- p
    }
    if (mod.row.Y == "VEI") {
      rowparY <- k + (p - 1)
    }
    if (mod.row.Y == "EVI") {
      rowparY <- 1 + k * (p - 1)
    }
    if (mod.row.Y == "VVI") {
      rowparY <- k * p
    }
    if (mod.row.Y == "EEE") {
      rowparY <- p * (p + 1) / 2
    }
    if (mod.row.Y == "VEE") {
      rowparY <- k - 1 + p * (p + 1) / 2
    }
    if (mod.row.Y == "EVE") {
      rowparY <- 1 + k * (p - 1) + p * (p - 1) / 2
    }
    if (mod.row.Y == "VVE") {
      rowparY <- k * p + p * (p - 1) / 2
    }
    if (mod.row.Y == "EEV") {
      rowparY <- p + k * p * (p - 1) / 2
    }
    if (mod.row.Y == "VEV") {
      rowparY <- k + (p - 1) + (k * p * (p - 1) / 2)
    }
    if (mod.row.Y == "EVV") {
      rowparY <- 1 + k * (p * ((p + 1) / 2) - 1)
    }
    if (mod.row.Y == "VVV") {
      rowparY <- k * p * (p + 1) / 2
    }
    
    if (mod.col.Y == "II") {
      colparY <- 0
    }
    if (mod.col.Y == "EI") {
      colparY <- r - 1
    }
    if (mod.col.Y == "VI") {
      colparY <- k * (r - 1)
    }
    if (mod.col.Y == "EE") {
      colparY <- r * ((r + 1) / 2) - 1
    }
    if (mod.col.Y == "VE") {
      colparY <- k * (r - 1) + r * (r - 1) / 2
    }
    if (mod.col.Y == "EV") {
      colparY <- (r - 1) + k * r * (r - 1) / 2
    }
    if (mod.col.Y == "VV") {
      colparY <- k * (r * ((r + 1) / 2) - 1)
    }
    
    weights <- k - 1
    pipar <- k*(k-1)
    
    npar <- coeff.par + rowparY + colparY + weights + pipar
    
    name <- c(mod.row.Y, mod.col.Y)
    
    # to be minimized
    
    BIC <- -2 * loglik.new + npar * log(num)

    ptm2 <- proc.time() - ptm
    time <- ptm2[3]
    
    return(list(
      name = name, prior = prior, coeff = coeff, sigmaUY = sigmaUY, sigmaVY = sigmaVY,
      PI = round(PI,digits = 3), loglik = loglik.new, mark = mark, check = check, npar = npar, iter = m.iter, time = time,
      BIC = BIC, class = groupT
    ))
  }
  
  for (g in 1:length(k)) {
    ptm <- proc.time()
    
    print(paste("Fitting Parsimonious MV-HMRMFCs with k =", k[g]))
    
    cluster <- makeCluster(nThreads, type = "SOCK")
    registerDoSNOW(cluster)
    
    pb <- progress_bar$new(
      format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
      total = pt.mod,
      complete = "=", # Completion bar character
      incomplete = "-", # Incomplete bar character
      current = ">", # Current bar character
      width = 100
    )
    progress <- function(n) {
      pb$tick()
    }
    opts <- list(progress = progress)
    
    oper[[g]] <- foreach(l = 1:pt.mod, .combine = "comb", .multicombine = TRUE, .init = list(list()), .options.snow = opts) %dopar% {
      res <- tryCatch(Pars_MVN_FMR(Y = Y, X = X, k = k[g], init.par = init.par[[1]][[g]][[1]][[init.par[[5]][l]]],
                                   mod.row.Y = as.character(list.comb[l, 1]), mod.col.Y = as.character(list.comb[l, 2]),
                                   tol = tol, tol2 = tol2, maxit = maxit, maxit2 = maxit2), error = function(e) {NA})
      
      list(res)
    }
    
    stopCluster(cluster)
    registerDoSEQ()
    
    ptm2 <- proc.time() - ptm
    time[g] <- ptm2[3]
  }
  
  return(list(results = oper, c.time = time, models = list.comb))
} # fitting function 

extract.bestM.FMR <- function(results,top=1) {
  
  k <- length(results[["results"]])
  num.mod <- length(results[["results"]][[1]][[1]])
  list.mod <- results[["models"]]
  list.mod2 <- do.call("rbind", replicate(k, list.mod, simplify = FALSE))
  count.k <- sort(rep(1:k, num.mod))
  count.mod <- rep(1:num.mod, k)
  list.mod3 <- data.frame(list.mod2, count.k, count.mod)
  
  allBIC <- numeric(k * num.mod)
  
  cont <- 0
  
  for (j in 1:k) {
    for (i in 1:num.mod) {
      if (!all(is.na(results[["results"]][[j]][[1]][[i]]))) {
        cont <- cont + 1
        allBIC[cont] <- -results[["results"]][[j]][[1]][[i]][["BIC"]]
      } else {
        cont <- cont + 1
        allBIC[cont] <- NA
      }
    }
  }
  
  topBIC <- which(allBIC>=sort(allBIC, decreasing = T)[top])
  topBIC.order <- order(allBIC[topBIC], decreasing = T)
  tempBIC <- list.mod3[topBIC[topBIC.order], ]
  bestBIC <- vector(mode = "list", length = top)
  for (i in 1:top) {
    bestBIC[[i]] <- results[["results"]][[as.numeric(tempBIC[i,3])]][[1]][[as.numeric(tempBIC[i,4])]]
  }
  
  return(bestBIC = bestBIC)
  
} # extract best fitting model 


