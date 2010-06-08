#CalMaTe's Test Script
# - - - - - - - - - - - - - - - - - - - - - - -
# load the test data
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#
# In this case, the variable data has the J=100 SNPs in the I=40 samples used in the manuscript
# - - - - - - - - - - - - - - - - - - - - - - -

library("calmate");
library("R.utils");

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate"); 
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Scaling the initial data to CN scale
thetaR <- rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE);
C <- 2*theta/thetaR;

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta);
dimnames(data)[[2]][2] <- "fracB";

# It returns a list where the elements of the list are the SNPs and their total
# copy number and fracB      
dataC <- calmateByTotalAndFracB(data);

# Allele specific copy numbers
dataAC <- dataC[,1,]*(1-dataC[,2,]);
dataBC <- dataC[,1,]*dataC[,2,];

# Comparing allele specific copy number results before and after CalMaTe calibration for sample 3 
nSample = 3;                 
plot(C[,1,nSample],C[,2,nSample], xlim = c(0,3), ylim = c(0,3));
points(dataAC[,nSample],dataBC[,nSample], col="blue");

# Comparing fracB results before and after CalMaTe calibration for sample 3 
nSample = 3;                 
x11();plot(C[,2,nSample], ylim = c(0,1));
points(dataC[,2,nSample], col="blue");
