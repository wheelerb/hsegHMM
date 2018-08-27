################################################################################################
##                                         hsegHMM:                                           ##
##                               Hidden Markov Model-Based                                    ##
##                   Allele-Specific Copy Number Alteration Analysis                          ##
##                           accounting for hypersegmentation                                 ##
##                     (R-code supported by Hyoyoung Choo-Wosoba)                             ##
##                            This r-code shows hsegHMM-N model                               ##
##                            of allele-specific SCNA analysis                                ##
################################################################################################

# History: Jun 29 2017 Initial coding

hsegHMM_N <- function(logR, logOR, purity=0.8, ploidy=1.5, logR.var=0.5, logOR.var=0.5, 
                    genoStates=c("", "A", "AA", "AB", "AAB", "AAA", "AAAB", 
                               "AABB", "AAAA", "AAAAB", "AAABB", "AAAAA"),
                    prob0=NULL, transProb=NULL, maxiter=100, stopTol=0.01, minLogOR2=1e-6,
                    optim.control=list(trace=0)) {

  if (length(logR) != length(logOR)) stop("ERROR: logR and logOR must have the same length")
  zk           <- unique(removeWhiteSpace(genoStates))
  J            <- length(zk)
  if (!J) stop("ERROR with genoStates")
  names(zk)    <- NULL
  lr           <- logR
  logor2       <- logOR^2
  pr0          <- purity
  lpr0         <- log((1-pr0)/pr0)
  sigma20_lr   <- logR.var
  lsigma20_lr  <- log(sigma20_lr) 
  lpsi0        <- log(ploidy, base=2)
  lsigma20_lor <- log(logOR.var) 
  r0           <- prob0
  if (is.null(r0)) r0 <- rep(1/J, J)
  if (length(r0) != J) stop("ERROR with prob0")
  hatp0        <- transProb
  if (is.null(hatp0)) hatp0 <- matrix(c(rep(c(1-(J-1)/5000,rep(1/5000,J)),(J-1)),1-(J-1)/5000),J,J) 
  if (nrow(hatp0) != J) stop("ERROR with transProb")
  if (ncol(hatp0) != J) stop("ERROR with transProb")
  M            <- as.integer(maxiter)
  if (M < 0) stop("ERROR: with maxiter")

  # Remove missing. logOR can have missing values
  temp   <- is.finite(lr)
  temp[is.na(temp)] <- FALSE
  lr     <- lr[temp]
  logor2 <- logor2[temp] 

  # Replace small values of logor2 with minLogOR2
  temp   <- logor2 < minLogOR2
  temp[is.na(temp)] <- FALSE
  if (any(temp)) logor2[temp] <- minLogOR2
  n      <- length(lr)

  # Compute mz, pz, ctz0 and check for errors with the genotype states
  temp1       <- callGenoStates(zk) 
  ctz0        <- temp1$ctz0
  mz          <- temp1$mz
  pz          <- temp1$pz
  unq_alleles <- temp1$unq_alleles
  
  # Initialize
  converged  <- FALSE
  notMissing <- is.finite(logor2)
  noMissing  <- all(notMissing)

  # conditional expectation of log-likelihood function given the posterior probability 
  logL_ <- function(parm){

    pr0 <- 1/(1+exp(parm[2]))
   
    # Modify for calculation of the hessian
    if (length(parm) > 3) {
      sigma20_lr <- parm[3]
      p3         <- parm[4] # We are including sigma20_lr, so shift
    } else {
      # sigma20_lr is defined outside this function
      p3         <- parm[3] # Optim call has 3 parms, sigma20_lr has closed-form solution
    }

    t1      <- -0.5*log(2*pi*sigma20_lr)
    t2      <- 0.5/sigma20_lr
    t3      <- exp(-p3)
    p1      <- parm[1]
    t4VecN  <- logor2*t3
    ncpVecJ <- (log((mz*pr0+(1-pr0)))-log(pz*pr0+(1-pr0)))^2*t3
    t5VecJ  <- log((2+pr0*(ctz0-2)), base=2)
   
    dmat <- dchisq(rep(t4VecN, times=J), df=1, ncp=rep(ncpVecJ, each=n), log=T)
    dim(dmat) <- c(n, J)

    if (noMissing) {
      MAT <- pzks1*(t1 - t2*(matrix(lr, byrow=FALSE, nrow=n, ncol=J) - matrix(t5VecJ, byrow=TRUE, nrow=n, ncol=J) + p1)^2 - p3 + dmat)
    } else {
      MAT <- pzks1*(t1 - t2*(matrix(lr, byrow=FALSE, nrow=n, ncol=J) - matrix(t5VecJ, byrow=TRUE, nrow=n, ncol=J) + p1)^2)
      MAT[notMissing, ] <- MAT[notMissing, , drop=FALSE] + pzks1[notMissing, , drop=FALSE]*(-p3 + dmat[notMissing, , drop=FALSE])
    } 
    logf <- sum(MAT)

    return(-logf)   

  } # END: logL_ 

  # Local function for computing ca_all and akh_all
  locFunca <- function() {
    c1             <- 1/sum(a1d)
    a1h            <- c1*a1d
    akh_all        <- matrix(NA, n,J)
    akh_all[1,]    <- a1h
    ca_all         <- numeric(n)
    ca_all[1]      <- c1

    t1             <- sqrt(sigma20_lr)
    t2             <- exp(-lsigma20_lor)
    tmat           <- dnorm(rep(lr, times=J), rep(mu0_lr, each=n), t1)
    dim(tmat)      <- c(n, J)
    ncpVecJ        <- mu0_lor^2*t2
    t2dchiMat      <- t2*dchisq(rep(logor2*t2, times=J), df=1, ncp=rep(ncpVecJ, each=n))
    dim(t2dchiMat) <- c(n, J)
    tmat[notMissing, ] <- tmat[notMissing, , drop=FALSE]*t2dchiMat[notMissing, , drop=FALSE]
    
    # Note: a1h gets updated after each k
    for (k in 2:n){
      a1hatMat    <- matrix(a1h, byrow=FALSE, nrow=J, ncol=J)*hatp0
      akd         <- colSums(a1hatMat)*tmat[k, ]
      cka         <- 1/sum(akd)
      ca_all[k]   <- cka
      akh         <- cka*akd                       
      akh_all[k,] <- akh
      a1h         <- akh
    }

    list(ca_all=ca_all, akh_all=akh_all)

  } # END: locFunca

  # Local function for computing cb_all and bkh_all
  locFuncb <- function() {

    bnd         <- rep(1,J)
    cn          <- 1/sum(bnd)
    bnh         <- cn*bnd
    cb_all      <- numeric(n)
    cb_all[n]   <- cn
    bkh_all     <- matrix(NA,n,J)
    bkh_all[n,] <- bnh

    t1             <- sqrt(sigma20_lr)
    t2             <- exp(-lsigma20_lor)
    tmat           <- dnorm(rep(lr, times=J), rep(mu0_lr, each=n), t1)
    dim(tmat)      <- c(n, J)
    ncpVecJ        <- mu0_lor^2*t2
    t2dchiMat      <- t2*dchisq(rep(logor2*t2, times=J), df=1, ncp=rep(ncpVecJ, each=n))
    dim(t2dchiMat) <- c(n, J)
    tmat[notMissing, ] <- tmat[notMissing, , drop=FALSE]*t2dchiMat[notMissing, , drop=FALSE]

    for (k in (n-1):1) {
      bnhatMat    <- matrix(bnh, byrow=TRUE, nrow=J, ncol=J)*hatp0
      tempmat     <- matrix(tmat[k+1, ], byrow=TRUE, nrow=J, ncol=J)
      bkd         <- rowSums(bnhatMat*tempmat)
      ckb         <- 1/sum(bkd)
      cb_all[k]   <- ckb
      bnh         <- ckb*bkd
      bkh_all[k,] <- bnh
    } 

    list(cb_all=cb_all, bkh_all=bkh_all)

  } # END: locFuncb

  # Local function for updating matrix of transition probs
  locFunc_hatp0 <- function() {
    
    t1 <- sqrt(sigma20_lr)    
    t2 <- exp(-lsigma20_lor)
    t3 <- sum(akh_all[n,])
    tveca     <- reverseCumSum(log(ca_all))
    tvecb     <- reverseCumSum(log(cb_all))
    tvec      <- exp(tveca[3:n]-tvecb[2:(n-1)])*ca_all[2:(n-1)]
    ncpVecJ   <- mu0_lor^2*t2
    tmat      <- dnorm(rep(lr, times=J), rep(mu0_lr, each=n), t1)
    dim(tmat) <- c(n, J)
    MAT       <- tmat*bkh_all/t3
    tmat      <- dchisq(rep(logor2*t2, times=J), df=1, ncp=rep(ncpVecJ, each=n))
    dim(tmat) <- c(n, J)
    MAT[notMissing, ] <- MAT[notMissing, , drop=FALSE]*t2*tmat[notMissing, , drop=FALSE]

#pzzks[,,k-1] <- if (is.na(logor2[k])) {
#exp(sum(log(ca_all[(k+1):n]))-sum(log(cb_all[k:n])))*ca_all[k]*
#(hatp0*(akh_all[k-1,]%*%matrix(bkh_all[k,]*dnorm(lr[k],mu0_lr, t1), nrow=1)))/t3
#} else {
#exp(sum(log(ca_all[(k+1):n]))-sum(log(cb_all[k:n])))*ca_all[k]*
#(hatp0*(akh_all[k-1,]%*%matrix(bkh_all[k,]*(dnorm(lr[k],mu0_lr, t1)*t2*
#dchisq(logor2[k]*t2,df=1, ncp=ncpVecJ)), nrow=1)))/t3
#}
#pzzks[,,(n-1)] <- if (is.na(logor2[n])) {
#sapply(1:J, function(i) sapply(1:J, function(j) 
# exp(sum(log(ca_all[n:n])))*hatp0[j,i]*akh_all[n-1,j]*dnorm(lr[n],mu0_lr[i], sqrt(sigma20_lr))/sum(akh_all[n,])))
#} else {
#sapply(1:J, function(i) sapply(1:J, function(j) 
# exp(sum(log(ca_all[n:n])))*hatp0[j,i]*akh_all[n-1,j]*dnorm(lr[n],mu0_lr[i], sqrt(sigma20_lr))*
# exp(-lsigma20_lor)*dchisq(logor2[n]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor))/sum(akh_all[n,])))
#}
  
    rm(tmat, tveca, tvecb)
    gc()

    # adjusted joint posterior probability of states (genotypes)
    pzzks <- array(NA, c(J,J,n-1))
    for( k in 2:(n-1)){
      pzzks[,,k-1] <- tvec[k-1]*(hatp0*(akh_all[k-1,]%*%matrix(MAT[k, ], nrow=1)))
    }
    
    vec <- ca_all[n]*dnorm(lr[n], mu0_lr, t1)/sum(akh_all[n,])
    if (notMissing[n]) vec <- vec*t2*dchisq(logor2[n]*t2, df=1, ncp=ncpVecJ)
    pzzks[, , (n-1)] <- hatp0*matrix(akh_all[n-1,], nrow=J, ncol=J, byrow=FALSE)*matrix(vec, nrow=J, ncol=J, byrow=TRUE)

    for (g in 1:J) {
      sumg <- sum(pzzks[g,,])
      for (l in 1:J) {
        hatp0[g, l] <- sum(pzzks[g,l,])/sumg
      }
    }

    hatp0

  } # END: locFunc_hatp0

  # Calculate mean of logR and logOR
  #mu0_lr <- sapply(1:J, function(j)  log((2+pr0*(ctz0[j]-2)), base=2)-lpsi0)
  mu0_lr  <- log((2+pr0*(ctz0-2)), base=2) - lpsi0
  mu0_lor <- log((mz*pr0+(1-pr0)))-log(pz*pr0+(1-pr0))

  ################# The proposed hidden Markov Model (HMM) procedure ####################################################################################
  
  logL_Y <- numeric(M+1)
  # Forward algorithm based on scaled HMM as a part of E-step in E-M algorithm
  #a1d <- if(is.na(logor2[1])) {sapply(1:J,function(i) dnorm(lr[1],mu0_lr[i], sqrt(sigma20_lr))*r0[i])} else {
  #  sapply(1:J,function(i) dnorm(lr[1],mu0_lr[i], sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor), df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor))*r0[i])
  #}

  if (is.na(logor2[1])) {
    a1d <- dnorm(lr[1], mu0_lr, sqrt(sigma20_lr))*r0
  } else {
    a1d <- dnorm(lr[1], mu0_lr, sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor), df=1, ncp=mu0_lor^2*exp(-lsigma20_lor))*r0         
  }

  # Old code
  #c1 <- 1/sum(a1d)
  #a1h <- c1*a1d
  #akh_all <- matrix(NA, n,J)
  #akh_all[1,] <- a1h
  #ca_all <- numeric(n)
  #ca_all[1] <- c1
  #for (k in 2:n){
  #  akd <- if(is.na(logor2[k])) {sapply(1:J, function(i)  sum(sapply(1:J, function(j) a1h[j]*hatp0[j,i]*dnorm(lr[k],mu0_lr[i], sqrt(sigma20_lr)))))} else
  #  {sapply(1:J, function(i) sum(sapply(1:J, function(j) a1h[j]*hatp0[j,i]*dnorm(lr[k],mu0_lr[i], sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[k]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor)))))}
  #  cka <- 1/sum(akd)
  #  ca_all[k] <- cka
  #  akh <- cka*akd                       
  #  akh_all[k,] <- akh
  #  a1h <- akh
  #}
  #old1 <- ca_all
  #old2 <- akh_all

  # New code
  temp1   <- locFunca()
  ca_all  <- temp1$ca_all
  akh_all <- temp1$akh_all

  logL_Y[1] <- sum(sum(-log(ca_all[1:n]))+log(sum(akh_all[n,]))) # calculate logLikelihood value, logL(Y|all the first set of estimates)
  cat(paste("Iteration 0: loglike = ", logL_Y[1],"\n", sep=""))

  if (!M) {
    temp <- list(loglike=logL_Y[1], 
       allele1=unq_alleles[1], allele2=unq_alleles[2],  
       alleleFreq1=mz, alleleFreq2=pz, copyNumber=ctz0,
       genoStates=zk, prob0=r0, transProb=hatp0)
    return(temp)
  }  


  ############## Iteration procedure including E-M algorithm ##############################################################################################
  for (mm in 1:M){
    # Backward algorithm based on scaled HMM as a part of E-step in E-M algorithm
  
    # Old code
    #bnd <- rep(1,J)
    #cn <- 1/sum(bnd)
    #bnh <- cn*bnd
    #cb_all <- numeric(n)
    #cb_all[n] <- cn
    #bkh_all <- matrix(NA,n,J)
    #bkh_all[n,] <- bnh
    #for (k in (n-1):1) {
    #  bkd <- if(is.na(logor2[k+1])) {sapply(1:J, function(l) sum(sapply(1:J, function(i) bnh[i]*hatp0[l,i]*dnorm(lr[k+1],mu0_lr[i], sqrt(sigma20_lr)))))} else
    #{                                sapply(1:J, function(l) sum(sapply(1:J, function(i) bnh[i]*hatp0[l,i]*dnorm(lr[k+1],mu0_lr[i], sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[k+1]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor)))))}
    #
    #  ckb <- 1/sum(bkd)
    #  cb_all[k] <- ckb
    #  bkh <- ckb*bkd
    # bkh_all[k,] <- bkh
    #  bnh <- bkh
    #} 
    #old1 <- cb_all
    #old2 <- bkh_all

    # New code
    temp1   <- locFuncb() 
    cb_all  <- temp1$cb_all
    bkh_all <- temp1$bkh_all

    # scaled posterior probability of states (genotypes)
    #pzks <- t(sapply(1:n, function(i) {
    #  pzki <- akh_all[i,]*bkh_all[i,]/sum(akh_all[n,])
    #  return(pzki)} ))
    pzks <- akh_all*bkh_all/sum(akh_all[n,])

    # adjusted posterior probability 
    # Old code
    #pzks1 <- t(cbind(sapply(1:(n-1), function(i) exp(sum(log(ca_all[(i+1):n]))-sum(log(cb_all[i:n])))*pzks[i,]), pzks[n,]/cb_all[n]))

    # New code
    pzks1            <- matrix(data=NA, nrow=n, ncol=J)
    colnames(pzks1)  <- zk
    pzks1[n ,]       <- pzks[n,]/cb_all[n]
    temp1            <- reverseCumSum(log(ca_all))
    temp2            <- reverseCumSum(log(cb_all))
    pzks1[1:(n-1), ] <- exp(temp1[2:n] - temp2[1:(n-1)])*pzks[1:(n-1), ]

    # adjusted joint posterior probability of states (genotypes)
    # Old code
    #pzzks <- array(NA, c(J,J,n-1))
    #for( k in 2:(n-1)){
    #  pzzks[,,k-1] <- if (is.na(logor2[k])) {exp(sum(log(ca_all[(k+1):n]))-sum(log(cb_all[k:n])))*ca_all[k]*(hatp0*
    #              (akh_all[k-1,]%*%matrix(bkh_all[k,]*dnorm(lr[k],mu0_lr, sqrt(sigma20_lr)), nrow=1)))/sum(akh_all[n,])} else
    #                {exp(sum(log(ca_all[(k+1):n]))-sum(log(cb_all[k:n])))*ca_all[k]*(
    #                  hatp0*(akh_all[k-1,]%*%matrix(bkh_all[k,]*(dnorm(lr[k],mu0_lr, sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[k]*exp(-lsigma20_lor),df=1, 
    #                    ncp=mu0_lor^2*exp(-lsigma20_lor))), nrow=1)))/sum(akh_all[n,])}
    #}
    #pzzks[,,(n-1)] <- if (is.na(logor2[n])) {sapply(1:J, function(i)  sapply(1:J, function(j) exp(sum(log(ca_all[n:n])))*hatp0[j,i]*akh_all[n-1,j]*dnorm(lr[n],mu0_lr[i], sqrt(sigma20_lr))/sum(akh_all[n,])))} else
    #  {sapply(1:J, function(i) sapply(1:J, function(j) exp(sum(log(ca_all[n:n])))*hatp0[j,i]*akh_all[n-1,j]*dnorm(lr[n],mu0_lr[i], sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[n]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor))/sum(akh_all[n,])))}
    #hatp00 <- sapply(1:J, function(g) sapply(1:J, function(l) sum(pzzks[l,g,])/(sum(pzzks[l,,]))))
    
    # New code
    hatp0 <- locFunc_hatp0()

    # estimate the variance of logR
    #sigma20_lr <- sum(sapply(1:J, function(i) sum(pzks1[,i]*(lr-mu0_lr[i])^2)))/sum(pzks1) 
    sigma20_lr <- sum(pzks1*(matrix(lr, nrow=n, ncol=J, byrow=FALSE)-matrix(mu0_lr, nrow=n, ncol=J, byrow=TRUE))^2)/sum(pzks1) 
    r0         <- pzks1[1,]/sum(pzks1[1,]) # estimate the initial probability of states

    # obtain all the other estimates to maximize the conditional expectation of logLikelihood with optim function 
    m_est <- optim(c(lpsi0,lpr0,lsigma20_lor), logL_, method="L-BFGS-B", control=optim.control)

    lpsi0        <- m_est$par[1]     # log2(estimated ploidy)
    lpr0         <- m_est$par[2]     # logit(estimated tumor purity)
    lsigma20_lor <- m_est$par[3]     # log(estimated variance of logOR)
    pr0          <- 1/(1+exp(lpr0))  # estimated tumor purity

    # update mean functions of logR and logOR given the updated estimates
    #mu0_lr <- sapply(1:J, function(j)  log((2+pr0*(ctz0[j]-2)), base=2)-lpsi0)
    mu0_lr  <- log((2+pr0*(ctz0-2)), base=2)-lpsi0
    mu0_lor <- log((mz*pr0+(1-pr0)))-log(pz*pr0+(1-pr0))
  
    # Forward algorithm based on scaled HMM with the updated estimates
    #a1d <- if(is.na(logor2[1])) {sapply(1:J,function(i) dnorm(lr[1],mu0_lr[i], sqrt(sigma20_lr))*r0[i])} else {
    #  sapply(1:J,function(i) dnorm(lr[1],mu0_lr[i], sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor), df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor))*r0[i])
    #}

    if (notMissing[1]) {
      a1d <- dnorm(lr[1], mu0_lr, sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor), df=1, ncp=mu0_lor^2*exp(-lsigma20_lor))*r0           
    } else {
      a1d <- dnorm(lr[1], mu0_lr, sqrt(sigma20_lr))*r0               
    }

    # old code
    #c1 <- 1/sum(a1d)
    #a1h <- c1*a1d
    #akh_all <- matrix(NA, n,J)
    #akh_all[1,] <- a1h
    #ca_all <- numeric(n)
    #ca_all[1] <- c1
    #for (k in 2:n){
    #  akd <- if(is.na(logor2[k])) {sapply(1:J, function(i)  sum(sapply(1:J, function(j) a1h[j]*hatp0[j,i]*dnorm(lr[k],mu0_lr[i], sqrt(sigma20_lr)))))} else
    #  {sapply(1:J, function(i) sum(sapply(1:J, function(j) a1h[j]*hatp0[j,i]*dnorm(lr[k],mu0_lr[i], sqrt(sigma20_lr))*exp(-lsigma20_lor)*dchisq(logor2[k]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor)))))}
    #  cka <- 1/sum(akd)
    #  ca_all[k] <- cka
    #  akh <- cka*akd                       
    #  akh_all[k,] <- akh
    #  a1h <- akh
    #}

    # New code
    temp1   <- locFunca()
    ca_all  <- temp1$ca_all
    akh_all <- temp1$akh_all

    # update the logLikelihood function, logL(Y|all the updated estimates)
    ll2 <- sum(sum(-log(ca_all[1:n]))+log(sum(akh_all[n,]))) 
    if (!is.finite(ll2)) stop("ERROR: log-likelihood is not finite")

    logL_Y[mm+1] <- ll2
    cat(paste("Iteration ", mm, ": loglike = ", ll2,"\n", sep=""))

    # convergence criteria by the difference between the updated logL and previous logL #
    if (abs(logL_Y[mm]-ll2) < stopTol) {
      converged <- TRUE
      break
    }
    
  } # END: for (mm in 1:M)

  # Find genotype status which gives the maximum posterior probability at each location 
  #idx_hgtype <- sapply(1:n, function(k) which.max(pzks1[k,]))
  idx_hgtype <- max.col(pzks1)
  #hat_ctz <- sapply(1:n, function(k) ctz0[idx_hgtype[k]])
  hat_ctz <- ctz0[idx_hgtype]

  # estimate expectation of logR and logOR based on estimates acquired from hsegHMM
  #hat_logr  <- sapply(1:n, function(k) log((2+pr0*(hat_ctz[k]-2)), base=2)-lpsi0)
  hat_logr  <- log((2+pr0*(hat_ctz-2)), base=2)-lpsi0
  hat_logor <- log(mz[idx_hgtype]*pr0+(1-pr0))-log(pz[idx_hgtype]*pr0+(1-pr0))

  # Compute AIC and BIC (4 parms)
  parms <- c(lpsi0, lpr0, sigma20_lr, lsigma20_lor)
  AIC   <- get_AIC(ll2, length(parms), J)
  BIC   <- get_BIC(ll2, length(parms), J, logR)

  # Get the hessian and covariance matrix (unscaled)
  hess   <- optimHess(parms, logL_)
  covmat <- get_AsymCovMat("N", parms, hess)

  list(converged=converged, loglike=ll2, 
       allele1=unq_alleles[1], allele2=unq_alleles[2],  
       alleleFreq1=mz, alleleFreq2=pz, copyNumber=ctz0,
       post.prob=pzks1,
       which.max.post.prob=idx_hgtype,
       logR_hat=hat_logr, logOR_hat=hat_logor,
       purity_hat=pr0, ploidy_hat=2^lpsi0, logR.var_hat=sigma20_lr, logOR.var_hat=exp(lsigma20_lor), 
       genoStates=zk, prob0=r0, transProb=hatp0, AIC=AIC, BIC=BIC, covariance=covmat)

} # END: hsegHMM_N
