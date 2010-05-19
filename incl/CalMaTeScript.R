#CalMaTe's script

workingDir <- "directory"
setwd(workingDir)

library(aroma.affymetrix)
library(MASS)
source("CalMaTeWeighted.R")
source("refineCN_rlmWeighted.R")

# - - - - - - - - - - - - - - - - - - - - - - -
# setup dataset and chip names
# - - - - - - - - - - - - - - - - - - - - - - -

log <- Arguments$getVerbose(-10, timestamp=TRUE);
projectName <- "projectName"
chipType <- "chipType"

# - - - - - - - - - - - - - - - - - - - - - - -
# define CDF file
# - - - - - - - - - - - - - - - - - - - - - - -
# Get CDF
cdf <- AffymetrixCdfFile$byChipType(chipType);
print(cdf);
cs <- AffymetrixCelSet$byName(projectName, cdf=cdf)
print(cs)
gi <- getGenomeInformation(cdf)
print(gi)
si <- getSnpInformation(cdf)
print(si)

#if using cross allele calibration
acc <- AllelicCrosstalkCalibration(cs, model="CRMAv2")
print(acc)
csC <- process(acc, verbose=verbose)

bpn <- BasePositionNormalization(csC, target="zero");
print(bpn)
csN <- process(bpn, verbose=verbose);

plm <- RmaCnPlm(csN, mergeStrands=TRUE, combineAlleles=FALSE)
fit(plm, verbose=verbose);

ces <- getChipEffects(plm);

# Process Chr22 for now
units <- getUnitsOnChromosome(gi, c(22));

dataCRMA <- extractTotalAndFreqB(ces, units = units);

#just in case there are na or nan
ind <- is.na(dataCRMA[,2,1]);
dataCRMA <- dataCRMA[!ind,,];

ind <- is.nan(dataCRMA[,2,1]);
dataCRMA <- dataCRMA[!ind,,];

#allele specific data

dataA <- dataCRMA[,1,]*(1-dataCRMA[,1,])
dataB <- dataCRMA[,1,]*(dataCRMA[,1,])

#in the case there are no normal references we do not give Refs as argument
#in other case we can give the references as a logical vector/matrix or even
#an index vector (see CalMaTe code for more information on this option)

#it returns a list where the elements of the list are the SNPs and their allele
#specific copy number values.
dataCalMaTe <- calMaTeWeighted(dataA, dataB);

#compare freqB results before and after CalMaTe calibration
pos <- getPositions(gi, units = units);
plot(pos[!ind], dataCRMA[,2,1])

#ugly way to do this...

calDataSp1 <- seq(0, nrow(dataCRMA))
for(ii in 1:nwor(dataCRMA)){
  aux <- dataCalMaTe[[1]][,1]
  calDataSp1 <- aux[1]/(aux[1]+aux[2])
}

x11();plot(pos[!ind], calDataSp1,ylim =c(0,1));
