%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  refineCN.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{refineCN.list}
\alias{refineCN.list}
\alias{refineCN.list}


\title{CalMaTe calibration function}

\description{
 CalMaTe calibration function
}

\usage{\method{refineCN}{list}(input, fB1=0.33, fB2=0.66, maxIter=50, ..., verbose=FALSE)}

\arguments{
 \item{input}{A list with two elements. First a 2xI \code{\link[base]{numeric}} array, where 2 is the number of alleles,
              and I is the number of samples. Second an array pointed out the reference samples.}
 \item{fB1}{Lower heterozygous threshold. Initially set to .33.}
 \item{fB2}{Higher heterozygous threshold. Initially set to .66.}
 \item{maxIter}{Maximum number of iterations for "rlm" function. Initially set to 50.}
 \item{...}{Not used.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns an 2xI \code{\link[base]{numeric}} array.
}

\examples{
# Load example (thetaA,thetaB) signals for one SNP in I samples and no given references
library("R.utils");
library("aroma.affymetrix");
path <- system.file("exData", package="calmate");
oneSNPdata <- loadObject("oneSNPdata,2x40.Rbin", path=path);

# Calibrating the initial data with our method
calibratedData <- refineCN(oneSNPdata);

# Normalizing the initial data to compare it to our method
inputData <- oneSNPdata[[1]]$inputData[[1]];          
nSamples <- length(inputData)/2;
dataA <- inputData[1:nSamples];
dataB <- inputData[(nSamples+1):(2*nSamples)];
fracB <- (dataB)/(dataA+dataB)
dataNorm <- 2*(dataA+dataB)/median(dataA+dataB);
dataA <- dataNorm*(1-fracB)
dataB <- dataNorm*fracB

# Plot both results to see the difference
plot(dataA, dataB, xlim = c(0,3),ylim =c(0,3))
points(calibratedData[1,], calibratedData[2,],col="red")


}

\seealso{
 To calibrate (total,fracB) data,
 see \code{\link[calmate:calmateByTotalAndFracB.list]{*calmateByTotalAndFracB}()}.
}
\keyword{methods}
