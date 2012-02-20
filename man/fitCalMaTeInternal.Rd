%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitCalMaTeInternal.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitCalMaTeInternal}
\alias{fitCalMaTeInternal}

\alias{fitCalMaTeV1}
\alias{fitCalMaTeV2}

\title{Algorithms to fit the CalMaTe model}

\description{
 Algorithms to fit the CalMaTe model.
 \emph{Note: These are internal functions of the package.
       They should not be used elsewhere.}
}

\usage{
  fitCalMaTeV1(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...)
  fitCalMaTeV2(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...)
}

\arguments{
 \item{dataT}{A 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of allele specific copy numbers (ASCNs),
    where 2 is the number alleles and I is the number of samples.}
 \item{references}{A \code{\link[base]{integer}} \code{\link[base]{vector}} with elements in [1,I] specifying
    which samples should be used as the reference set.}
 \item{fB1, fB2}{Thresholds for calling genotypes AA, AB, BB from the
    allele B fractions.}
 \item{maxIter}{The maximum number of iterations without converging
    before the algorithm quits.}
 \item{...}{Not used.}
}

\value{
  Returns a 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of calibrated ASCNs.
}

\section{Flavor "v1"}{
  \emph{To be documented.}
}

\section{Flavor "v2"}{
  \emph{To be documented.}
}

\seealso{
  These functions are called by \code{\link{calmateByThetaAB}}().
}


\keyword{documentation}
\keyword{internal}
