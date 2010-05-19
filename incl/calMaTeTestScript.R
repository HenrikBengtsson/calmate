#CalMaTe's Test Script

library(aroma.affymetrix)
library("calmate")
library(MASS)

source("weightedCalMaTe.R")
source("refineCN_rlmWeighted.R")
source("weightedCalMaTeByASCN.R")

# - - - - - - - - - - - - - - - - - - - - - - -
# load the test data
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
# - - - - - - - - - - - - - - - - - - - - - - -

path <- system.file("exData/", package="calmate"); 
data <- R.utils::loadObject("myASCNData.Rbin", path=path);

# it returns a list where the elements of the list are the SNPs and their allele
# specific copy number values.

dataCalMaTe <- weightedCalMaTeByASCN(data);

#compare freqB results before and after CalMaTe calibration for sample 3 

nSp = 3;                 
plot(data[,1,nSp]/median(data[,1,nSp]),data[,2,nSp]/median(data[,2,nSp]), xlim = c(0,3), ylim = c(0,3))
points(dataCalMaTe[,1,nSp],dataCalMaTe[,2,nSp], col="blue")
