#CalMaTe's Test Script
# - - - - - - - - - - - - - - - - - - - - - - -
# load the test data
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#
# In this case, the variable data has the J=100 SNPs in the I=40 samples used in the manuscript
# - - - - - - - - - - - - - - - - - - - - - - -

library("calmate")

# Load example (thetaA,thetaB) signals
library("R.utils");
path <- system.file("exData", package="calmate"); 
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path);

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta);

# it returns a list where the elements of the list are the SNPs and their allele
# specific copy number values.

dataC <- weightedCalMaTeByTotalAndFracB(data);

#compare freqB results before and after CalMaTe calibration for sample 3 

nSp = 3;                 
plot(data[,1,nSp]/median(data[,1,nSp]),data[,2,nSp]/median(data[,2,nSp]), xlim = c(0,3), ylim = c(0,3))
points(dataC[,1,nSp],dataC[,2,nSp], col="blue")
