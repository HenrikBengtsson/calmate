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

 \usage{\method{calmateByTotalAndFracB}{array}(data, references=NULL, ..., refAvgFcn=NULL, verbose=FALSE)}

 \arguments{
  \item{data}{An Jx2xI \code{\link[base]{numeric}} \code{\link[base]{array}}, where J is the number of loci,
              2 is total and fracB (in that order, if unnamed), and
              I is the number of samples.}
  \item{references}{A \code{\link[base]{logical}} or \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying which
     samples should be used as the reference set.
     By default, all samples are considered. If not NULL at least 3 samples.}
  \item{...}{Additional arguments passed to \code{\link[calmate:calmateByThetaAB.array]{*calmateByThetaAB}()}.}
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
   Returns an Jx2xI \code{\link[base]{numeric}} \code{\link[base]{array}}
   with the same dimension names as argument \code{data}.
 }

 \examples{
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


# Assert that it also works with a single unit
dummy <- calmateByTotalAndFracB(data[1,,,drop=FALSE]);
stopifnot(length(dim(dummy)) == 3);
}

 \references{
  [1] M. Ortiz-Estevez, A. Aramburu, H. Bengtsson, P. Neuvial and A. Rubio, \emph{CalMaTe: A method and software to improve allele-specific copy number of SNP arrays for downstream segmentation}, Bioinformatics, 2012 (in print).
 \cr
 }

 \seealso{
  To calibrate (thetaA,thetaB) or (CA,CB) signals,
  see \code{\link[calmate:calmateByThetaAB.array]{*calmateByThetaAB}()}.
 }
\keyword{methods}
