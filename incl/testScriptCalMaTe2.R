#CalMaTe's Test Script
# - - - - - - - - - - - - - - - - - - - - - - -
# In this case the SNPs from chromosome 22 are calibrated by CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - -

workingDir <- "directory";
setwd(workingDir);

library(aroma.affymetrix);
library(MASS);
library(calmate);

# - - - - - - - - - - - - - - - - - - - - - - -
# Setup dataset and chip names
# - - - - - - - - - - - - - - - - - - - - - - -

projectName <- "projectName";
chipType <- "chipType";

# - - - - - - - - - - - - - - - - - - - - - - -
# Define CDF file
# - - - - - - - - - - - - - - - - - - - - - - -
# Get CDF
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);
cs <- AffymetrixCelSet$byName(projectName, cdf=cdf);
print(cs);
gi <- getGenomeInformation(cdf);
print(gi);
si <- getSnpInformation(cdf);
print(si);

# If using cross allele calibration
acc <- AllelicCrosstalkCalibration(cs, model="CRMAv2");
print(acc);
csC <- process(acc, verbose=verbose);

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn);
csN <- process(bpn, verbose=verbose);

plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE);
fit(plm, verbose=verbose);

ces <- getChipEffects(plm);

# Process Chr22 for now
units <- getUnitsOnChromosome(gi, c(22));
dataCRMA <- extractTotalAndFreqB(ces, units = units);

# In case there are na or nan
ind <- is.na(dataCRMA[,2,1]);
dataCRMA <- dataCRMA[!ind,,];

ind <- is.nan(dataCRMA[,2,1]);
dataCRMA <- dataCRMA[!ind,,];

# Allele specific copy numbers
dataAB <- dataCRMA;
colnames(dataAB)<-c("A","B")
dataAB[,1,] <- dataCRMA[,1,]*(1-dataCRMA[,2,]);
dataAB[,2,] <- dataCRMA[,1,]*(dataCRMA[,2,]);

# In the case there are no normal references we do not give Refs as argument
# in other case we can give the references as a logical vector/matrix or even
# an index vector (see CalMaTe code for more information on this option)
# It returns a 3 dimensional matrix (Jx2xI) of allele specific copy numbers, 
# where J is the number of SNPs and the I the number of samples
dataCalMaTe <- calmateByThetaAB(dataAB);

# Comparing allele specific copy number results before and after CalMaTe calibration
# Scaling the initial data to Copy Number values
dataAB <- 2*dataAB/rowMedians(dataAB[,1,]+dataAB[,2,]);

# Plot to "random" arrays
Clim <- c(0,3);
subplots(4, ncol=2, byrow=TRUE);
for (ii in c(1,5)) {
  plot(dataAB[,,ii], cex = .3, pch = 16, xlim=Clim, ylim=Clim);
  title(main=dimnames(C)[[3]][ii]);
  plot(dataCalMaTe[,,ii], cex = .3, pch = 16, xlim=Clim, ylim=Clim);
  title(main=sprintf("%s\ncalibrated", dimnames(C)[[3]][ii]));
}

