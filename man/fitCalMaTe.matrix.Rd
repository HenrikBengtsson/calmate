%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitCalMaTe.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{fitCalMaTe.matrix}
\alias{fitCalMaTe.matrix}
\alias{fitCalMaTe}

\title{Calibrates SNP loci according to the CalMaTe method}

\description{
 Calibrates SNP loci according to the CalMaTe method.
 \emph{Note: This is an internal function of the package, which is kept
  only kept to provide easy access to the internal fit functions.
  It it actually not elsewhere in the package, and should nor by others.}
}

\usage{
\method{fitCalMaTe}{matrix}(dataT, references, flavor=c("v2", "v1"), ...)
}

\arguments{
 \item{dataT}{A 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of allele specific copy numbers (ASCNs),
    where 2 is the number alleles and I is the number of samples.}
 \item{references}{A \code{\link[base]{integer}} \code{\link[base]{vector}} with elements in [1,I] specifying
    which samples should be used as the reference set.}
 \item{...}{Additional arguments passed to the internal fit functions.}
 \item{flavor}{A \code{\link[base]{character}} string specifying which internal fit function
    (flavor of the CalMaTe algorithm) to use for fitting the model.}
}

\value{
  Returns a 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of calibrated ASCNs.
}

\section{Flavors}{
  For backward compatibility, we try to keep all major versions of
  the CalMaTe algorithm available.  Older versions can be used by
  specifying argument \code{flavor}.
  For more information about the different flavors,
  see \code{\link{fitCalMaTeInternal}}.
}

\references{
 [1] M. Ortiz-Estevez, A. Aramburu, H. Bengtsson, P. Neuvial and A. Rubio, \emph{CalMaTe: A method and software to improve allele-specific copy number of SNP arrays for downstream segmentation}, Bioinformatics, 2012 (in print).
 \cr
}

\seealso{
  For further information on the internal fit functions,
  see \code{\link{fitCalMaTeInternal}}.
}


\keyword{methods}
\keyword{internal}
