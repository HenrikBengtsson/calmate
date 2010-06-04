%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  calmate.R
% 
% on Fri Jun 04 17:27:21 2010.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{calmate.matrix}
\alias{calmate.matrix}
\alias{calmate.matrix}


\title{Internal Calmate function}

\description{
 Internal Calmate function
}

\usage{\method{calmate}{matrix}(dataA, dataB, refs=0, maxIter=50, ..., verbose=FALSE)}

\arguments{
 \item{dataA}{Matrix with the number of copies of A allele of J SNPs in the I samples.}
 \item{dataB}{Matrix with the number of copies of B allele of J SNPs in the I samples.}
 \item{refs}{Vector indicating the reference samples.}
 \item{maxIter}{Maximum number of iterations for "rlm" function in refineCN. Initially set to 50.}
 \item{...}{Additional arguments passed to
        \code{\link[calmate:refineCN.matrix]{*refineCN}()}.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
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
CC <- calmate(theta[,"A",],theta[,"B",]);
matAB <- matrix( unlist(CC), nrow=(2*length(CC)),ncol=ncol(CC[[1]]))

# Plot to "random" arrays
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
for (ii in c(1,5)) {
  plot(C[,,ii], xlim=Clim, ylim=Clim);
  title(main=dimnames(C)[[3]][ii]);
  plot(matAB[seq(1,nrow(matAB),2),ii], matAB[seq(2,nrow(matAB),2),ii], xlim=Clim, ylim=Clim);
  title(main=sprintf("\%s\ncalibrated", dimnames(C)[[3]][ii]));
}



}

\seealso{
 To calibrate (total,fracB) data,
 see \code{\link[calmate:calmateByTotalAndFracB.matrix]{*calmateByTotalAndFracB}()}.
 see \code{\link[calmate:refineCN.matrix]{*refineCN}()}.
}
\keyword{methods}
