%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  calmate.R
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

\usage{\method{calmate}{matrix}(dataA, dataB, refs=0, ..., verbose=FALSE)}

\arguments{
 \item{dataA}{A \code{\link[base]{numeric}} JxI copy-number \code{\link[base]{matrix}},
       where J is the number of SNPs and I is the number samples.}
 \item{dataB}{A \code{\link[base]{numeric}} JxI copy-number \code{\link[base]{matrix}}.}
 \item{refs}{A \code{\link[base]{logical}} \code{\link[base]{vector}} indicating the reference samples.}
 \item{...}{Additional arguments passed to \code{\link{refineCN.list}}().}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[base]{list}} with J elements (SNPs) each containing
  a 2xI \code{\link[base]{numeric}} \code{\link[base]{matrix}}.
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
res <- calmate(theta[,"A",],theta[,"B",]);
CC <- unlist(res, use.names=FALSE);
dim(CC) <- dim(C)[c(2,3,1)];
CC <- aperm(CC, perm=c(3,1,2));
dimnames(CC) <- dimnames(C);

# Plot two "random" arrays
Clim <- c(0,4);
subplots(4, ncol=2, byrow=FALSE);
for (ii in c(1,5)) {
  sampleName <- dimnames(C)[[3]][ii];
  sampleLabel <- sprintf("Sample #\%d ('\%s')", ii, sampleName);
  plot(C[,,ii], xlim=Clim, ylim=Clim);
  title(main=sampleLabel);
  plot(CC[,,ii], xlim=Clim, ylim=Clim);
  title(main=sprintf("\%s\ncalibrated", sampleLabel));
}
}

\seealso{
 See \code{\link{refineCN.list}}().
}


\keyword{methods}
\keyword{internal}
