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

 \title{Normalize allele-specific copy numbers (CA,CB)}

 \description{
  Normalize allele-specific copy numbers (CA,CB).
 }

 \usage{\method{calmateByThetaAB}{array}(data, references=NULL, ..., truncate=FALSE, verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} array, where J is the number of SNPs,
          2 is the number of alleles, and I is the number of samples.}
  \item{references}{A \code{\link[base]{logical}} \code{\link[base]{vector}} of length I, or an index \code{\link[base]{vector}}
          with values in [0,I]. If \code{\link[base]{NULL}}, all samples are considered to
          be reference samples (==1:I).}
  \item{...}{Additional arguments passed to internal
             \code{calmate().}}
  \item{truncate}{If \code{\link[base:logical]{TRUE}}, final ASCNs are forced to be non-negative,
     which preserving the total CNs.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns an Jx2xI \code{\link[base]{numeric}} array.
 }                                 

 \examples{
library("calmate");
library("R.utils");

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Calculate (CA,CB)
thetaR <- rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE);
C <- 2*theta/thetaR;

# Calibrate (CA,CB) by CalMaTe
CC <- calmateByThetaAB(theta);

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
  To calibrate (total,fracB) data,
  see \code{\link[calmate:calmateByTotalAndFracB.array]{*calmateByTotalAndFracB}()}.
 }
\keyword{methods}
