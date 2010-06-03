%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  thetaAB2TotalAndFracB.R
% 
% on Thu Jun 03 12:10:04 2010.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{thetaAB2TotalAndFracB.array}
\alias{thetaAB2TotalAndFracB.array}
\alias{thetaAB2TotalAndFracB.array}

 \alias{thetaAB2TotalAndFracB}
 \alias{totalAndFracB2ThetaAB}
 \alias{totalAndFracB2ThetaAB.array}

 \title{Converts an Jx2xI array between (thetaA,thetaB) and (total,fracB) formats}

 \description{
  Converts an Jx2xI array between (thetaA,thetaB) and (total,fracB) formats.
 }

 \usage{\method{thetaAB2TotalAndFracB}{array}(data, ..., verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} array, where J is the number of SNPs,
          2 is the number of alleles, and I is the number of samples.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{
# Load example (thetaA,thetaB) signals
library("R.utils");
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

data <- thetaAB2TotalAndFracB(theta);
str(data);

theta2 <- totalAndFracB2ThetaAB(data);
str(theta2);

stopifnot(all.equal(theta2, theta));
}

 \value{
   Returns an Jx2xI \code{\link[base]{numeric}} array.
 }


\keyword{methods}
\keyword{internal}
