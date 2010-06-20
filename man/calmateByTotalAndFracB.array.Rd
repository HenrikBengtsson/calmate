%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  calmateByTotalAndFracB.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{calmateByTotalAndFracB.array}
\alias{calmateByTotalAndFracB.array}
\alias{calmateByTotalAndFracB.array}

 \alias{calmateByTotalAndFracB}

 \title{Normalize allele-specific copy numbers (total,fracB)}

 \description{
  Normalize allele-specific copy numbers (total,fracB), where total is the total (non-polymorphic) signal and
  fracB is the allele B fraction.
  It is only loci with a non-missing (\code{\link[base]{NA}}) fracB value that are
  considered to be SNPs and normalized by CalMaTe.  The other loci
  are left untouched.
 }

 \usage{\method{calmateByTotalAndFracB}{array}(data, ..., verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} \code{\link[base]{array}}, where J is the number of loci,
                      2 is total and fracB, and I is the number of samples.}
  \item{...}{Additional arguments passed to
         \code{\link[calmate:calmateByThetaAB.array]{*calmateByThetaAB}()}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns an Jx2xI \code{\link[base]{numeric}} \code{\link[base]{array}}.
 }

 \examples{
library("calmate");
library("R.utils");

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate");
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta);

# Calibrate (total,fracB) by CalMaTe
dataC <- calmateByTotalAndFracB(data);

# Calculate copy-number ratios
theta <- data[,"total",];
thetaR <- rowMedians(theta, na.rm=TRUE);
data[,"total",] <- 2*theta/thetaR;

# Plot two "random" arrays
Clim <- c(0,4);
Blim <- c(0,1);
subplots(4, ncol=2, byrow=FALSE);
for (ii in c(1,5)) {
  sampleName <- dimnames(data)[[3]][ii];
  sampleLabel <- sprintf("Sample #\%d ('\%s')", ii, sampleName);
  plot(data[,,ii], xlim=Clim, ylim=Blim);
  title(main=sampleLabel);
  plot(dataC[,,ii], xlim=Clim, ylim=Blim);
  title(main=sprintf("\%s\ncalibrated", sampleLabel));
}
}

 \seealso{
  To calibrate (thetaA,thetaB) or (CA,CB) signals,
  see \code{\link[calmate:calmateByThetaAB.array]{*calmateByThetaAB}()}.
 }
\keyword{methods}
