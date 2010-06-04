# Load example (thetaA,thetaB) signals
library("R.utils");
path <- system.file("exData", package="calmate"); 
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta);

# CalMaTe calibration
dataC <- calmateByTotalAndFracB(data);

theta <- data[,"total",];
thetaR <- rowMedians(theta, na.rm=TRUE);
data[,"total",] <- 2*theta/thetaR;

# Plot to "random" arrays
Clim <- c(0,4);
Blim <- c(0,1);
subplots(4, ncol=2, byrow=TRUE);
for (ii in c(1,5)) {
  plot(data[,,ii], xlim=Clim, ylim=Blim);
  title(main=dimnames(C)[[3]][ii]);
  plot(dataC[,,ii], xlim=Clim, ylim=Blim);
  title(main=sprintf("%s\ncalibrated", dimnames(C)[[3]][ii]));
}
