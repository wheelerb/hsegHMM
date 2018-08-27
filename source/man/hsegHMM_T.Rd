
\name{hsegHMM_T}
\alias{hsegHMM_T}
\title{The hsegHMM procedure for t-distributed log(Ratio) values}
\description{
   A hidden Markov model approach that accounts for hypersegmentation 
    for allele-specific SCNA analysis for t-distributed log(Ratio) values                             
}
\usage{ hsegHMM_T(logR, logOR, purity=0.8, ploidy=1.5, logR.var=0.5, logOR.var=0.5, df=3,
                    genoStates=c("", "A", "AA", "AB", "AAB", "AAA", "AAAB", 
                               "AABB", "AAAA", "AAAAB", "AAABB", "AAAAA"),
                    prob0=NULL, transProb=NULL, maxiter=100, stopTol=0.01, minLogOR2=1e-6,
                    df.min=0.0001, df.max=100, optim.control=list(trace=0))
}
\arguments{
  \item{logR}{ Vector of (non-missing) log(Ratio) values. No default }
  \item{logOR}{ Vector of log(Odds Ratio) values. Note that this vector can have
               missing values and will be squared inside the function. No default. }
  \item{purity}{ Initial value for the tumor purity. The default is 0.8.}
  \item{ploidy}{ Initial value for the ploidy. The default is 1.5.}
  \item{logR.var}{ Initial value for the variance component of logR. The default is 0.5.}
  \item{logOR.var}{Initial value for the variance component of logOR. The default is 0.5.}
  \item{df}{Initial value for the degrees of freedom. The default is 3.}
  \item{genoStates}{Character vector of hidden genotype states. The default is \cr
     c("", "A", "AA", "AB", "AAB", "AAA", "AAAB", \cr
       "AABB", "AAAA", "AAAAB", "AAABB", "AAAAA")}. 
  \item{prob0}{NULL or a vector of initial probabilities for \code{genoStates}. 
               If NULL, then it will be set to \code{rep(1/length(genoStates), length(genoStates))}. 
               The default is NULL.}
  \item{transProb}{NULL or a matrix of transition probabilities. 
               If NULL, then it will be set to 
               \code{matrix(c(rep(c(1-(J-1)/5000,rep(1/5000,J)),(J-1)),1-(J-1)/5000),J,J)},
               where \code{J = length(genoStates)}.
               The default is NULL.}
  \item{maxiter}{Maximum number of iterations for the algorithm. The default is 100.}
  \item{stopTol}{Stopping tolerance for the algorithm. The algorithm will stop
                when two succesive log-likelihood values differ by less than \code{stopTol}.
                The default is 0.01.}
  \item{minLogOR2}{Minimum value for \code{logOR^2} to prevent numerical difficulties in
                  the algorithm.
                The default is 1e-6.}
  \item{df.min}{Lower bound for the \code{df} parameter in the optimization. The default is 0.0001.}
  \item{df.max}{Upper bound for the \code{df} parameter in the optimization. The default is 100.}
  \item{optim.control}{List for the \code{control} option in the \code{\link[stats]{optim}} function. 
             The default is list(trace=0).}
}
\details{
 Missing values are allowed for \code{logOR} values as long as \code{logR} is observed. 
 The optimization is performed
  using the L-BFGS-B method in the \code{\link[stats]{optim}} function. 
  To prevent numerical errors in the \code{\link[base]{gamma}} function, \code{df.min} should be 
  a positive value and \code{df.max} a finite positive value. For the other parameters,
  the lower and upper bounds can be infinite.
 See the reference for details of the algorithm.
}
\value{
  A list with the following names and descriptions: 

  \itemize{
    \item \code{converged} Convergence status
    \item \code{loglike} Final value of the log-likelihood
    \item \code{allele1} Allele 1
    \item \code{allele2} Allele 2
    \item \code{alleleFreq1} Frequency of allele 1
    \item \code{alleleFreq2} Frequency of allele 2
    \item \code{copyNumber} Copy number
    \item \code{post.prob} Matrix of posterior probabilities for each genotype state
    \item \code{which.max.post.prob} Indices for the genotype status which
                  gives the maximum posterior probability.
    \item \code{logR_hat} The expected value of \code{logR}
    \item \code{logOR_hat} The expected value of \code{logOR}
    \item \code{purity_hat} The expected value of \code{purity}
    \item \code{ploidy_hat} The expected value of \code{ploidy}
    \item \code{logR.var_hat} The expected value of \code{logR.var}
    \item \code{logOR.var_hat} The expected value of \code{logOR.var}
    \item \code{df_hat} The expected value of \code{df}
    \item \code{genoStates} The genotype states
    \item \code{prob0} The initial probability of the genotype states
    \item \code{transProb} The matrix of transition probabilities
    \item \code{AIC} Akaike information criterion
    \item \code{BIC} Bayesian information criterion
    \item \code{covariance} Covariance matrix for all parameters
  } % END: itemize

}
\references{Choo-Wosoba, H., Albert, P.S., Zhu, B.
  hsegHMM: Hidden Markov Model-based Allele-specific Copy Number Alteration Analysis
  Accounting for Hypersegmentation
}
\author{Hyoyoung Choo-Wosoba, Paul S. Albert, and Bin Zhu <bin.zhu@nih.gov>}
\examples{
  data(hseg_data, package="hsegHMM")

  hsegHMM_T(lr, logor)

}
\keyword{ hsegHMM, SCNA } 
\seealso{\code{\link{hsegHMM_N}}}
