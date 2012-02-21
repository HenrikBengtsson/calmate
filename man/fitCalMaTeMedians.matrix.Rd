%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitCalMaTeMedians.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitCalMaTeMedians.matrix}
\alias{fitCalMaTeMedians.matrix}
\alias{fitCalMaTeMedians.matrix}

\alias{fitCalMaTeMedians}

\title{Calibrates SNP loci according to the CalMaTe method, not using "rlm" function.}

\description{
 Calibrates SNP loci according to the CalMaTe method, not using "rlm" function.
 \emph{Note: This is an internal function of the package, which is kept
  only kept to provide easy access to the internal fit functions.
  It it actually not elsewhere in the package, and should nor by others.}
}

\usage{\method{fitCalMaTeMedians}{matrix}(T, references, fB1=1/3, fB2=2/3, ...)}

\arguments{
 \item{T}{A 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of allele specific copy numbers (ASCNs),
    where 2 is the number alleles and I is the number of samples.}
 \item{references}{A \code{\link[base]{integer}} \code{\link[base]{vector}} with elements in [1,I] specifying
    which samples should be used as the reference set.}
 \item{fB1, fB2}{Thresholds for calling genotypes AA, AB, BB from the
    allele B fractions.}
 \item{...}{Additional arguments passed to the internal fit functions.}
}

\value{
  Returns a 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of calibrated ASCNs.
}

\keyword{methods}
