
# History: May 20 2017 Make code more efficient
#          Jun 02 2017 Add options for lower/upper bounds of df
#          Jun 09 2017 Allow missing values for logOR

hsegHMM_T <- function(logR, logOR, purity=0.8, ploidy=1.5, logR.var=0.5, 
                logOR.var=0.5, df=3, 
                genoStates=c("", "A", "AA", "AB", "AAB", "AAA", "AAAB", 
                    "AABB", "AAAA", "AAAAB", "AAABB", "AAAAA"),
                prob0=NULL, transProb=NULL, maxiter=100, stopTol=0.01, 
                minLogOR2=1e-6, df.min=0.0001, df.max=100, 
                optim.control=list(trace=0)) {

  if (length(logR) != length(logOR)) {
    stop("ERROR: logR and logOR must have the same length")
  }
  zk           <- unique(removeWhiteSpace(genoStates))
  J            <- length(zk)
  if (!J) stop("ERROR with genoStates")
  names(zk)    <- NULL
  lr           <- logR
  logor2       <- logOR^2
  pr0          <- purity
  lpr0         <- log((1-pr0)/pr0)
  v0           <- df
  sigma20_lr   <- logR.var
  lsigma20_lr  <- log(sigma20_lr) 
  lpsi0        <- log(ploidy, base=2)
  lsigma20_lor <- log(logOR.var) 
  r0           <- prob0
  if (is.null(r0)) r0 <- rep(1/J, J)
  if (length(r0) != J) stop("ERROR with prob0")
  hatp0        <- transProb
  if (is.null(hatp0)) {
    hatp0 <- matrix(c(rep(c(1-(J-1)/5000,rep(1/5000,J)),(J-1)),1-(J-1)/5000),J,J) 
  }
  if (nrow(hatp0) != J) stop("ERROR with transProb")
  if (ncol(hatp0) != J) stop("ERROR with transProb")
  M            <- as.integer(maxiter)
  if (M < 0) stop("ERROR: with maxiter")
  if (df.min < 0) stop("ERROR: df.min must be non-negative")
  if (df.max < df.min) stop("ERROR: df.max < df.min")

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

    pr0        <- 1/(1+exp(parm[2]))
    sigma20_lr <- exp(parm[3])
    #mu0_lr    <- sapply(1:J, function(j)  log((2+pr0*(ctz0[j]-2)), base=2)-parm[1])
    mu0_lr     <- log((2+pr0*(ctz0-2)), base=2)-parm[1]
    mu0_lor    <- log((mz*pr0+(1-pr0)))-log(pz*pr0+(1-pr0))
    
    parm4 <- parm[4]
    parm5 <- parm[5]
    temp1 <- -0.5*log(2*pi*sigma20_lr)+lgamma((parm5+1)/2)-lgamma(parm5/2)+0.5*parm5*log(parm5/2)
    temp2 <- 0.5*(parm5+1)
    emp4  <- exp(-parm4)
    ncp   <- mu0_lor^2*emp4      
    val   <- temp1 - temp2*log(0.5*(parm5 +(lr - matrix(mu0_lr, nrow=n, ncol=J, byrow=TRUE))^2/sigma20_lr))
    val   <- val*pzks1
    tvec  <- logor2*emp4
    dmat <- dchisq(rep(tvec, times=J), df=1, ncp=rep(ncp, each=n), log=T)
    dim(dmat) <- c(n, J)
    if (noMissing) {
      val <- val - pzks1*(parm4 - dmat)
    } else {
      val[notMissing, ] <- val[notMissing, , drop=FALSE] - pzks1[notMissing, , drop=FALSE]*(parm4 - dmat[notMissing, , drop=FALSE])
    }

    logf <- sum(val)

    return(-logf)   

  } # END: logL_ 

  # Local function for computing ca_all and akh_all
  locFunca <- function() {
    c1          <- 1/sum(a1d)
    a1h         <- c1*a1d
    akh_all     <- matrix(NA, n,J)
    akh_all[1,] <- a1h
    ca_all      <- numeric(n)
    ca_all[1]   <- c1
    temp1       <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tempg1      <- gamma(v0/2)
    temp2       <- 0.5*v0
    temp3       <- 2*sigma20_lr
    temp4       <- (v0+1)/2
    temp5       <- exp(-lsigma20_lor)
    tempncp     <- mu0_lor^2*temp5  # A vector

    # Get all chi-sq values in a matrix
    dmat <- dchisq(rep(logor2*temp5, times=J),df=1, ncp=rep(tempncp, each=n))
    dim(dmat) <- c(n, J)
    tmat <- temp1/(tempg1*(temp2+((lr-matrix(mu0_lr, nrow=n, ncol=J, byrow=TRUE))^2/(temp3)))^(temp4))
    if (noMissing) {
      tmat <- tmat*temp5*dmat
    } else {
      tmat[notMissing, ] <- tmat[notMissing, , drop=FALSE]*temp5*dmat[notMissing, , drop=FALSE]
    }

    # Note: a1h gets updated after each k
    for (k in 2:n) {
      tempmat     <- matrix(tmat[k,], nrow=J, ncol=J, byrow=TRUE)
      akd         <- colSums(a1h*hatp0*tempmat)
      cka         <- 1/sum(akd)
      ca_all[k]   <- cka
      a1h         <- cka*akd                       
      akh_all[k,] <- a1h
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
    temp1       <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tempg1      <- gamma(v0/2)
    temp2       <- 0.5*v0
    temp3       <- 2*sigma20_lr
    temp4       <- (v0+1)/2
    temp5       <- exp(-lsigma20_lor)
    tempncp     <- mu0_lor^2*temp5  # A vector

    # Get all chi-sq values in a matrix
    dmat <- dchisq(rep(logor2*temp5, times=J),df=1, ncp=rep(tempncp, each=n))
    dim(dmat) <- c(n, J)
    tmat <- temp1/(tempg1*(temp2+((lr-matrix(mu0_lr, nrow=n, ncol=J, byrow=TRUE))^2/(temp3)))^(temp4))
    if (noMissing) {
      tmat <- tmat*temp5*dmat
    } else {
      tmat[notMissing, ] <- tmat[notMissing, , drop=FALSE]*temp5*dmat[notMissing, , drop=FALSE]
    }

    for (k in (n-1):1) {
      tempmat     <- matrix(tmat[k+1, ], nrow=J, ncol=J, byrow=TRUE)
      tempmat2    <- matrix(bnh, nrow=J, ncol=J, byrow=TRUE)
      bkd         <- rowSums(tempmat2*hatp0*tempmat)
      ckb         <- 1/sum(bkd)
      cb_all[k]   <- ckb
      bnh         <- ckb*bkd
      bkh_all[k,] <- bnh
    } 

    list(cb_all=cb_all, bkh_all=bkh_all)

  } # END: locFuncb

  # Local function for updating matrix of transition probs
  locFunc_hatp0 <- function() {
    
    temp1   <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5
    tempg1  <- gamma(v0/2)
    temp2   <- 0.5*v0
    temp3   <- 2*sigma20_lr
    temp4   <- (v0+1)/2
    temp5   <- exp(-lsigma20_lor)
    tempncp <- mu0_lor^2*temp5  # A vector
    temp6   <- sum(akh_all[n,])
    tveca   <- reverseCumSum(log(ca_all))
    tvecb   <- reverseCumSum(log(cb_all))
    tvec    <- exp(tveca[3:n]-tvecb[2:(n-1)])*ca_all[2:(n-1)]/temp6

    # Get all dchisq values
    dchimat <- matrix(data=NA, nrow=n, ncol=J)
    for (j in 1:J) dchimat[, j] <- dchisq(logor2*temp5,df=1, ncp=tempncp[j])
    basemat <- temp2+((lr - matrix(mu0_lr, nrow=n, ncol=J, byrow=TRUE))^2/temp3)
    if (noMissing) {
      tempmat <- bkh_all*((temp1/(tempg1*basemat^temp4))*temp5*dchimat)
    } else {
      tempmat <- bkh_all*((temp1/(tempg1*basemat^temp4)))
      tempmat[notMissing, ] <- bkh_all[notMissing, ]*((temp1/(tempg1*basemat[notMissing, ]^temp4))*temp5*dchimat[notMissing, ])
    }
    rm(dchimat, basemat, tveca, tvecb)
    gc()

    # adjusted joint posterior probability of states (genotypes)
    pzzks <- array(NA, c(J,J,n-1))
    for( k in 2:(n-1)){
      pzzks[,,k-1] <- tvec[k-1]*(hatp0*(akh_all[k-1,]%*%matrix(tempmat[k, ], nrow=1)))
    }
    
    vec <- temp1/(tempg1*(temp2+((lr[n]-mu0_lr)^2/(temp3)))^(temp4))/temp6
    if (notMissing[n]) vec <- vec*temp5*dchisq(logor2[n]*temp5,df=1, ncp=mu0_lor^2*temp5)
    pzzks[, , (n-1)] <- ca_all[n]*hatp0*matrix(akh_all[n-1,], nrow=J, ncol=J, byrow=FALSE)*matrix(vec, nrow=J, ncol=J, byrow=TRUE)

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

  a1d <- if (is.na(logor2[1])) {
           (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lr)^2/(2*sigma20_lr)))^((v0+1)/2))*r0
         } else {
           (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lr)^2/(2*sigma20_lr)))^((v0+1)/2))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor),df=1, ncp=mu0_lor^2*exp(-lsigma20_lor))*r0                 
         }

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


  ############## Iteration procedure including E-M algorithm ######################################
  for (mm in 1:M){
    # Backward algorithm based on scaled HMM as a part of E-step in E-M algorithm
    temp1   <- locFuncb() 
    cb_all  <- temp1$cb_all
    bkh_all <- temp1$bkh_all

    # scaled posterior probability of states (genotypes)
    pzks <- akh_all*bkh_all/sum(akh_all[n,])

    # adjusted posterior probability 
    pzks1            <- matrix(data=NA, nrow=n, ncol=J)
    colnames(pzks1)  <- zk
    pzks1[n ,]       <- pzks[n,]/cb_all[n]
    temp1            <- reverseCumSum(log(ca_all))
    temp2            <- reverseCumSum(log(cb_all))
    pzks1[1:(n-1), ] <- exp(temp1[2:n] - temp2[1:(n-1)])*pzks[1:(n-1), ]

    hatp0 <- locFunc_hatp0()

    r0 <- pzks1[1,]/sum(pzks1[1,]) # estimate the initial probability of states
    # obtain all the other estimates to maximize the conditional expectation of logLikelihood with optim function 
    lower <- c(rep(-Inf,4), df.min)
    upper <- c(rep(Inf,4),  df.max)
    m_est <- optim(c(lpsi0,lpr0,lsigma20_lr, lsigma20_lor, v0), logL_, method="L-BFGS-B", control=optim.control,lower=lower, upper=upper)

    lpsi0        <- m_est$par[1]     # log2(estimated ploidy)
    lpr0         <- m_est$par[2]     # logit(estimated tumor purity)
    lsigma20_lr  <- m_est$par[3]     # log(estimated variance of logR)
    lsigma20_lor <- m_est$par[4]     # log(estimated variance of logoR)
    v0           <- m_est$par[5]     # degree of freedom for a t-distribution of logR
    sigma20_lr   <- exp(lsigma20_lr) # estimated variance of logR
    pr0          <- 1/(1+exp(lpr0))  # estimated tumor purity

    # update mean functions of logR and logOR given the updated estimates
    #mu0_lr <- sapply(1:J, function(j)  log((2+pr0*(ctz0[j]-2)), base=2)-lpsi0)
    mu0_lr  <- log((2+pr0*(ctz0-2)), base=2)-lpsi0
    mu0_lor <- log((mz*pr0+(1-pr0)))-log(pz*pr0+(1-pr0))
  
    # Forward algorithm based on scaled HMM with the updated estimates
    #a1d <- if(is.na(logor2[1])) {sapply(1:J,function(i) (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lr[i])^2/(2*sigma20_lr)))^((v0+1)/2))*r0[i]) } else {
    #                             sapply(1:J,function(i) (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lr[i])^2/(2*sigma20_lr)))^((v0+1)/2))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor))*r0[i])                   
    #}

    if (notMissing[1]) {
      a1d <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lr)^2/(2*sigma20_lr)))^((v0+1)/2))*exp(-lsigma20_lor)*dchisq(logor2[1]*exp(-lsigma20_lor),df=1, ncp=mu0_lor^2*exp(-lsigma20_lor))*r0                 
    } else {
      a1d <- (0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[1]-mu0_lr)^2/(2*sigma20_lr)))^((v0+1)/2))*r0                  
    }

    # Old code
    #c1 <- 1/sum(a1d)
    #a1h <- c1*a1d
    #akh_all <- matrix(NA, n,J)
    #akh_all[1,] <- a1h
    #ca_all <- numeric(n)
    #ca_all[1] <- c1
    #for (k in 2:n){
    #  akd <- if(is.na(logor2[k])) {sapply(1:J, function(i) sum(sapply(1:J, function(j) a1h[j]*hatp0[j,i]*(0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[k]-mu0_lr[i])^2/(2*sigma20_lr)))^((v0+1)/2)))))} else
    #                              {sapply(1:J, function(i) sum(sapply(1:J, function(j) a1h[j]*hatp0[j,i]*(0.5*v0)^(0.5*v0)*gamma((v0+1)/2)*(2*pi*sigma20_lr)^-0.5/(gamma(v0/2)*(0.5*v0+((lr[k]-mu0_lr[i])^2/(2*sigma20_lr)))^((v0+1)/2))*exp(-lsigma20_lor)*dchisq(logor2[k]*exp(-lsigma20_lor),df=1, ncp=mu0_lor[i]^2*exp(-lsigma20_lor)))))}
    #
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

  # Compute AIC and BIC 
  parms <- c(lpsi0, lpr0, lsigma20_lr, lsigma20_lor, v0)
  AIC   <- get_AIC(ll2, length(parms), J)
  BIC   <- get_BIC(ll2, length(parms), J, logR)

  # Get the hessian and covariance matrix (unscaled)
  hess   <- optimHess(parms, logL_)
  covmat <- get_AsymCovMat("T", parms, hess)

  list(converged=converged, loglike=ll2, 
       allele1=unq_alleles[1], allele2=unq_alleles[2],  
       alleleFreq1=mz, alleleFreq2=pz, copyNumber=ctz0,
       post.prob=pzks1,
       which.max.post.prob=idx_hgtype,
       logR_hat=hat_logr, logOR_hat=hat_logor,
       purity_hat=pr0, ploidy_hat=2^lpsi0, logR.var_hat=sigma20_lr, logOR.var_hat=exp(lsigma20_lor), df_hat=v0,
       genoStates=zk, prob0=r0, transProb=hatp0, AIC=AIC, BIC=BIC, covariance=covmat)

} # END: hsegHMM_T
