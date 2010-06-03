%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  weightedCalMaTe.R
% 
% on Thu Jun 03 16:07:14 2010.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{weightedCalMaTe.matrix}
\alias{weightedCalMaTe.matrix}
\alias{weightedCalMaTe.matrix}


\title{Internal Calmate function}

\description{
 Internal Calmate function
}

\usage{\method{weightedCalMaTe}{matrix}(dataA, dataB, refs=0, ..., verbose=FALSE)}

\arguments{
 \item{dataA}{Matrix with the number of copies of A allele of J SNPs in the I samples.}
 \item{dataB}{Matrix with the number of copies of B allele of J SNPs in the I samples.}
 \item{refs}{Vector indicating the reference samples.}
}

\value{
  Returns a list with J elements (number of SNPs) an each containing a 2xI \code{\link[base]{numeric}} array.
}

\examples{
# Load example (thetaA,thetaB) signals
library("R.utils");
library("aroma.affymetrix");
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Calculate (CA,CB)
thetaR <- rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE);
C <- 2*theta/thetaR;

# Calibrate (CA,CB)
CC <- weightedCalMaTe(theta[,"A",],theta[,"B",]);

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
 see \code{\link[CalMaTe2:weightedCalMaTeByTotalAndFracB.matrix]{*weightedCalMaTeByTotalAndFracB}()}.
 see \code{\link[CalMaTe2:refineCNrlmWeighted.matrix]{*refineCNrlmWeighted}()}.
}
\keyword{methods}
