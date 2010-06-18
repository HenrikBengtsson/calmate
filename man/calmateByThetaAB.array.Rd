%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  calmateByThetaAB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{calmateByThetaAB.array}
\alias{calmateByThetaAB.array}
\alias{calmateByThetaAB.array}

 \alias{calmateByThetaAB}

 \title{Internal CalMaTe fit function}

 \description{
  Internal CalMaTe fit function.
 }

 \usage{\method{calmateByThetaAB}{array}(data, ..., verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} array, where J is the number of SNPs,
          2 is the number of alleles, and I is the number of samples.}
  \item{...}{Additional arguments passed to internal
             \code{calmate().}}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns an Jx2xI \code{\link[base]{numeric}} array.
 }                                 

 \examples{
library(calmate);
library("R.utils");

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Calculate (CA,CB)
thetaR <- rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE);
C <- 2*theta/thetaR;

# Calibrate (CA,CB)
CC <- calmateByThetaAB(theta);

# Plot to "random" arrays
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
for (ii in c(1,5)) {
  plot(C[,,ii], xlim=Clim, ylim=Clim);
  title(main=dimnames(C)[[3]][ii]);
  plot(CC[,,ii], xlim=Clim, ylim=Clim);
  title(main=sprintf("\%s\ncalibrated", dimnames(C)[[3]][ii]));
}
}

 \seealso{
  To calibrate (total,fracB) data,
  see \code{\link[calmate:calmateByTotalAndFracB.array]{*calmateByTotalAndFracB}()}.
 }
\keyword{methods}
