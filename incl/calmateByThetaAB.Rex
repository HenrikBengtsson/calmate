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
  sampleLabel <- sprintf("Sample #%d ('%s')", ii, sampleName);
  plot(C[,,ii], xlim=Clim, ylim=Clim);
  title(main=sampleLabel);
  plot(CC[,,ii], xlim=Clim, ylim=Clim);
  title(main=sprintf("%s\ncalibrated", sampleLabel));
}


# Assert that it also works with a single unit
dummy <- calmateByThetaAB(theta[1,,,drop=FALSE]);
stopifnot(length(dim(dummy)) == 3);
