%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.package.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{calmate-package}
\alias{calmate-package}
\alias{calmate}
\docType{package}

\title{Package calmate}


 \description{
   A multi-array post-processing method of allele-specific copy-number estimates (ASCNs)
 }

 \section{Requirements}{
   This package depends on a set of packages that are all available
   via CRAN.  It has been tested and verified to run on all common
   operating systems on which R runs, including Linux, Windows and OSX.
 }

 \section{Installation and updates}{
   To install this package, do \code{install.packages("calmate")}.
 }

 \section{To get started}{
  \enumerate{
   \item To process SNP and non-polymorphic signals, see \code{\link{calmateByTotalAndFracB}}().  If you are working solely with SNP signals, \code{\link{calmateByThetaAB}}() is also available, but we recommend the former.
   \item For processing data in the aroma framework, see \code{\link{CalMaTeCalibration}}.
  }
 }

 \section{How to cite}{
   Please cite [1] when using CalMaTe.
 }

 \author{
  Maria Ortiz [aut, ctb],
  Ander Aramburu [ctb],
  Henrik Bengtsson [aut, cre, cph],
  Pierre Neuvial [aut, ctb],
  Angel Rubio [aut, ctb].
 }

 \section{License}{
  LGPL (>= 2.1).
 }

 \references{
  [1] M. Ortiz-Estevez, A. Aramburu, H. Bengtsson, P. Neuvial and A. Rubio, \emph{CalMaTe: A method and software to improve allele-specific copy number of SNP arrays for downstream segmentation}, Bioinformatics, 2012 (in print).
 \cr
 }
\keyword{package}
