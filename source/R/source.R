
# Function to remove leading/trailing white space
removeWhiteSpace <- function(str, leading=1, trailing=1) {

  if ((leading) && (trailing)) {
    ret <- gsub("^\\s+|\\s+$", "", str, perl=TRUE)
  } else if (leading) {
    ret <- gsub("^\\s+", "", str, perl=TRUE)
  } else if (trailing) {
    ret <- gsub("\\s+$", "", str, perl=TRUE)
  } else {
    ret <- str
  }

  ret

} # END: removeWhiteSpace

# Function to compute a cummaulative sum in "reverse order"
reverseCumSum <- function(vec) {

  n   <- length(vec)
  ret <- vec[n:1]
  ret <- cumsum(ret)
  ret <- ret[n:1]

  ret
}

# Function to check the genotype states
checkGenoStates <- function(gs, unq) {

  temp <- nchar(gs) > 0
  gs   <- gs[temp]
  gs2  <- gsub(unq[1], "@", gs, fixed=TRUE)
  gs2  <- gsub(unq[2], unq[1], gs2, fixed=TRUE)
  gs2  <- gsub("@", unq[2], gs2, fixed=TRUE)
    
  temp <- gs2 %in% gs
  if (any(temp)) {
    err <- gs2[temp]
    str <- paste(err, collapse=",", sep="")
    str <- paste("ERROR: with genoStates ", str, sep="")
    stop(str)
  }

  NULL

} # END: checkGenoStates

# Function to compute copy number, etc
callGenoStates <- function(zk) {

  J     <- length(zk)
  ctz0  <- nchar(zk)
  tlist <- strsplit(zk, "", fixed=TRUE)
  unq   <- unique(unlist(tlist))
  temp  <- nchar(unq) > 0
  unq   <- unq[temp]
  nunq  <- length(unq)
  if ((!nunq) || (nunq > 2)) stop("ERROR with genoStates")
  mz    <- rep(NA, J)
  pz    <- rep(0, J)
  # Compute mz, pz and also normalize the genotype states
  for (j in 1:J) {
    vec   <- tlist[[j]]
    mz[j] <- sum(vec %in% unq[1]) 
    if (nunq > 1) pz[j] <- sum(vec %in% unq[2])
    zk[j] <- paste(sort(vec), collapse="", sep="") 
  }  
  unq_alleles <- unq
  if (nunq == 1) unq_alleles <- c(unq_alleles, unq_alleles)
  if (nunq > 1) checkGenoStates(zk, unq_alleles) 

  list(ctz0=ctz0, mz=mz, pz=pz, unq_alleles=unq_alleles)

} # END: callGenoStates

# Function to compute the AIC
get_AIC <- function(loglike, nparm, nstate) {

  J   <- nstate
  aic <- 2*((J-1)*J+nparm+J-1)-2*loglike
  aic

} # END: get_AIC

# Function to compute the AIC
get_BIC <- function(loglike, nparm, nstate, logR) {

  J   <- nstate
  bic <- log(length(logR))*((J-1)*J+nparm+J-1)-2*loglike
  bic

} # END: get_BIC

# Function to return parameter names
get_VarNames <- function(WHICH) {

  if (WHICH == "T") {
    ret <- c("ploidy", "purity", "logR.var", "logOR.var", "df")
  } else {
    ret <- c("ploidy", "purity", "logR.var", "logOR.var")
  }
  ret

} # END: get_VarNames

# Function to return the score vector needed for asymptotic covariance matrix
get_scoreVec <- function(WHICH, parms) {

  # WHICH:   "T" or "N"

  # Order for case T:
  #m_est <- optim(c(lpsi0,lpr0,lsigma20_lr, lsigma20_lor, v0), logL_,
  # Order for case N:

  lpsi0        <- parms[1]
  lpr0         <- parms[2]
  lsigma20_lr  <- parms[3]
  lsigma20_lor <- parms[4]

  if (WHICH == "T") {
    ret <- c(2^lpsi0*log(2), -(1+exp(lpr0))^-2*exp(lpr0), 
    exp(lsigma20_lr), exp(lsigma20_lor), 1)
  } else {
    ret <- c(2^lpsi0*log(2), -(1+exp(lpr0))^-2*exp(lpr0), 1, 
    exp(lsigma20_lor))
  }
  
  ret

} # END: get_scoreVec

# Function to return the asymptotic covar matrix
get_AsymCovMat <- function(WHICH, parms, hess) {

  # WHICH:   "T" or "N"

  var_scale   <- try(solve(-hess))
  if (!("try-error" %in% class(var_scale))) {
    vec <- get_scoreVec(WHICH, parms) 
    mat <- diag(vec)
    ret <- t(mat) %*% var_scale %*% mat
  } else {
    ret <- matrix(data=NA, nrow=length(parms), ncol=length(parms)) 
  }
  cnames <- get_VarNames(WHICH)
  colnames(ret) <- cnames
  rownames(ret) <- cnames

  ret

} # END: get_AsymCovMat


