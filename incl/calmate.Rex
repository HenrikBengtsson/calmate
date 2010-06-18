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
matAB <- matrix(unlist(CC), nrow=(2*length(CC)), ncol=ncol(CC[[1]]));

# Plot to "random" arrays
Clim <- c(0,4);
subplots(4, ncol=2, byrow=TRUE);
for (ii in c(1,5)) {
  plot(C[,,ii], xlim=Clim, ylim=Clim);
  title(main=dimnames(C)[[3]][ii]);
  plot(matAB[seq(1,nrow(matAB),2),ii], matAB[seq(2,nrow(matAB),2),ii],
                                                 xlim=Clim, ylim=Clim);
  title(main=sprintf("%s\ncalibrated", dimnames(C)[[3]][ii]));
}
