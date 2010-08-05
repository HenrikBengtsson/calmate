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

 \usage{\method{calmateByThetaAB}{array}(data, references=NULL, ..., truncate=FALSE, refAvgFcn=NULL, verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} \code{\link[base]{array}}, where J is the number of SNPs,
          2 is the number of alleles, and I is the number of samples.}
  \item{references}{An index \code{\link[base]{vector}} in [1,I] or a \code{\link[base]{logical}} \code{\link[base]{vector}}
     of length I specifying which samples are used when calculating the
     reference signals.  If \code{\link[base]{NULL}}, all samples are used.}
  \item{...}{Additional arguments passed to internal \code{fitCalMaTe().}}
  \item{truncate}{If \code{\link[base:logical]{TRUE}}, final ASCNs are forced to be non-negative,
     which preserving the total CNs.}
  \item{refAvgFcn}{(optional) A \code{\link[base]{function}} that takes a JxI \code{\link[base]{numeric}} \code{\link[base]{matrix}}
     an argument \code{na.rm} and returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of length J.
     It should calculate some type of average for each of the J rows, e.g.
     \code{\link[matrixStats]{rowMedians}}.
     If specified, then the total copy numbers of the calibrated ASCNs
     are standardized toward (twice) the average of the total copy numbers
     of the calibrated reference ASCNs.}
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


# Assert that it also works with a single unit
dummy <- calmateByThetaAB(theta[1,,,drop=FALSE]);
stopifnot(length(dim(dummy)) == 3);
}

 \seealso{
  To calibrate (total,fracB) data,
  see \code{\link[calmate:calmateByTotalAndFracB.array]{*calmateByTotalAndFracB}()}.
 }
\keyword{methods}
