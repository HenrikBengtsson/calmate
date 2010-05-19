%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  weightedCalMaTeByASCN.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{weightedCalMaTeByASCN.array}
\alias{weightedCalMaTeByASCN.array}
\alias{weightedCalMaTeByASCN.array}

 \alias{weightedCalMaTeByASCN}

 \title{Internal CalMaTe fit function}

 \description{
  Internal CalMaTe fit function.
 }

 \usage{\method{weightedCalMaTeByASCN}{array}(data, ..., verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} array, where J is the number of SNPs,
          2 is the number of alleles, and I is the number of samples.}
  \item{...}{Additional arguments passed to internal
             \code{calMaTeWeighted().}}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns an Jx2xI \code{\link[base]{numeric}} array.
 }
\keyword{methods}
