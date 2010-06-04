%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  weightedCalMaTeByTotalAndFracB.R
% 
% on Fri Jun 04 17:27:21 2010.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{weightedCalMaTeByTotalAndFracB.array}
\alias{weightedCalMaTeByTotalAndFracB.array}
\alias{weightedCalMaTeByTotalAndFracB.array}

 \alias{weightedCalMaTeByTotalAndFracB}

 \title{Internal CalMaTe fit function}

 \description{
  Internal CalMaTe fit function.
 }

 \usage{\method{weightedCalMaTeByTotalAndFracB}{array}(data, ..., verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} array, where J is the number of SNPs,
                      2 is total and fracB, and I is the number of samples.}
  \item{...}{Additional arguments passed to
         \code{\link[calmate:weightedCalMaTeByThetaAB.array]{*weightedCalMaTeByThetaAB}()}.}
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

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta);

# CalMaTe calibration
dataC <- weightedCalMaTeByTotalAndFracB(data);

theta <- data[,"total",];
thetaR <- rowMedians(theta, na.rm=TRUE);
data[,"total",] <- 2*theta/thetaR;

# Plot to "random" arrays
Clim <- c(0,4);
Blim <- c(0,1);
subplots(4, ncol=2, byrow=TRUE);
for (ii in c(1,5)) {
  plot(data[,,ii], xlim=Clim, ylim=Blim);
  title(main=dimnames(C)[[3]][ii]);
  plot(dataC[,,ii], xlim=Clim, ylim=Blim);
  title(main=sprintf("\%s\ncalibrated", dimnames(C)[[3]][ii]));
}
}

 \seealso{
  To calibrate (thetaA,thetaB) data,
  see \code{\link[calmate:weightedCalMaTeByThetaAB.array]{*weightedCalMaTeByThetaAB}()}.
 }
\keyword{methods}
