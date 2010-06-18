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
  sampleLabel <- sprintf("Sample #%d ('%s')", ii, sampleName);
  plot(C[,,ii], xlim=Clim, ylim=Clim);
  title(main=sampleLabel);
  plot(CC[,,ii], xlim=Clim, ylim=Clim);
  title(main=sprintf("%s\ncalibrated", sampleLabel));
}
