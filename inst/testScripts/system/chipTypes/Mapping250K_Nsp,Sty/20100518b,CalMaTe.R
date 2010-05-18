###########################################################################
# Author: Henrik Bengtsson
# Created on: 2010-05-18
# Last updated: 2010-05-18
#
# Data:
# rawData/GSE12702/Mapping250K_Nsp/
#
# URL: http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE12702
###########################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("calmate");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

dataSet <- "Affymetrix_2006-TumorNormal";
tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
chipType <- "Mapping250K_Nsp";
dsT <- AromaUnitTotalCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
print(dsT);

dsB <- AromaUnitFracBCnBinarySet$byName(dataSet, tags=tags, chipType=chipType);
print(dsB);

stopifnot(identical(getNames(dsT), getNames(dsB)));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Study array #1 and chromosome 22
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
array <- 1;
chr <- 1;

chrTag <- sprintf("Chr%02d", chr);

# Extract the array
dfT <- getFile(dsT, array);
dfB <- getFile(dsB, array);
stopifnot(identical(getName(dfT), getName(dfB)));

# Extract the chromosome
ugp <- getAromaUgpFile(dsT);
units <- getUnitsOnChromosome(ugp, chr);
pos <- getPositions(ugp, units=units);

# Extract (total, fracB) signals
theta <- extractMatrix(dfT, units=units, drop=TRUE);
beta <- extractMatrix(dfB, units=units, drop=TRUE);

# Convert to allele-specific SNP signals (thetaA, thetaB)
thetaB <- theta * beta;
thetaA <- theta - thetaB;
thetaAB <- matrix(c(thetaA, thetaB), ncol=2, byrow=FALSE);
colnames(thetaAB) <- c("A", "B");
rm(thetaA, thetaB);

# Keep non-missing values
idxs <- which(is.finite(thetaAB[,"A"]) & is.finite(thetaAB[,"B"]));
res <- CalMaTeWeighted(thetaAB[idxs,"A"], thetaAB[idxs,"B"]);
res <- lapply(res, FUN=function(x) {
  x <- t(x);
  colnames(x) <- c("A", "B");
  x;
});
naValue <- as.double(NA);
thetaABC <- matrix(naValue, nrow=nrow(thetaAB), ncol=2);
colnames(thetaABC) <- c("A", "B");
thetaABC[idxs,] <- res[[1]];
thetaC <- thetaABC[,"A"] + thetaABC[,"B"];
betaC <- thetaABC[,"B"] / thetaC;


# Plot
ylim <- c(-0.2, 1.2);
subplots(2, ncol=1);

fracB <- RawAlleleBFractions(beta, x=pos);
plot(fracB, pch=".", ylim=ylim);
stext(side=3, pos=0, getName(dfT));
stext(side=3, pos=1, chrTag);

fracBC <- RawAlleleBFractions(betaC, x=pos);
plot(fracBC, pch=".", ylim=ylim);
stext(side=3, pos=0, getName(dfT));
stext(side=3, pos=1, chrTag);


###########################################################################
# HISTORY:
# 2010-05-18
# o Created.
###########################################################################
