
\name{hsegHMM}
\alias{hsegHMM}
\docType{package}
\title{
 Hidden Markov Model-based Allele-specific Copy Number Alteration Analysis Accounting for Hypersegmentation
}
\description{
 An R package for the hsegHMM model of allele-specific SCNA analysis
}
\details{
Somatic copy number alternation (SCNA) is a common feature of the
cancer genome and is associated with cancer etiology and prognosis. The allele-specific
SCNA analysis of a tumor sample aims to identify the allele-specific copy numbers of
both alleles, adjusting for the ploidy and the tumor purity. Next generation sequencing
platforms produce abundant read counts at the base-pair resolution across the exome or
whole genome which is susceptible to hypersegmentation, a phenomenon where
numerous regions with very short length are falsely identified as SCNA.
This package employs a robust hidden Markov model approach that accounts for
hypersegmentation for allele-specific SCNA analysis, and an efficient E-M
algorithm procedure that uses a forward-backward algorithm for evaluating the E-step.
The main functions that perform this method are 
\code{\link{hsegHMM_N}} for normally distributed log(ratio) values, and 
\code{\link{hsegHMM_T}} for t-distributed log(ratio) values.
}
\references{Choo-Wosoba, H., Albert, P.S., Zhu, B.
  hsegHMM: Hidden Markov Model-based Allele-specific Copy Number Alteration Analysis
  Accounting for Hypersegmentation
}
\author{Hyoyoung Choo-Wosoba, Paul S. Albert, and Bin Zhu <bin.zhu@nih.gov>}
\keyword{package}