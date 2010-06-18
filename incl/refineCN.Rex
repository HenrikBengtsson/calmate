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
fracB <- dataB / (dataA+dataB);
dataNorm <- 2*(dataA+dataB)/median(dataA+dataB);
dataA <- dataNorm * (1-fracB);
dataB <- dataNorm * fracB;

# Plot both results to see the difference
plot(dataA, dataB, xlim=c(0,3), ylim=c(0,3));
points(calibratedData[1,], calibratedData[2,], col="red");


