%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  weightedCalMaTeByThetaAB.R
% 
% on Thu Jun 03 11:18:03 2010.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{weightedCalMaTeByThetaAB.array}
\alias{weightedCalMaTeByThetaAB.array}
\alias{weightedCalMaTeByThetaAB.array}

 \alias{weightedCalMaTeByThetaAB}

 \title{Internal CalMaTe fit function}

 \description{
  Internal CalMaTe fit function.
 }

 \usage{\method{weightedCalMaTeByThetaAB}{array}(data, ..., verbose=FALSE)}

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

 \examples{
# Load example (thetaA,thetaB) signals
library("R.utils");
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Calculate (CA,CB)
thetaR <- rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE);
C <- 2*theta/thetaR;

# Calibrate (CA,CB)
CC <- weightedCalMaTeByThetaAB(theta);

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
  see \code{\link[CalMaTe2:weightedCalMaTeByTotalAndFracB.array]{*weightedCalMaTeByTotalAndFracB}()}.
 }
\keyword{methods}
