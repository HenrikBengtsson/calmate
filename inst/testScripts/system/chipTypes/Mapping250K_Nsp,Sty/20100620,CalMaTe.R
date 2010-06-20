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
library("aroma.core");
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
# CalMaTe
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsList <- list(total=dsT, fracB=dsB);
asn <- CalMaTeNormalization(dsList);
print(asn);

ugp <- getAromaUgpFile(dsT);
chr <- 17;
units <- getUnitsOnChromosome(ugp, chr);

dsNList <- process(asn, units=units, verbose=verbose);
print(dsNList);

ii <- 1;
df <- getFile(dsList$fracB, ii);
dfN <- getFile(dsNList$fracB, ii);

beta <- extractRawAlleleBFractions(df, chromosome=chr);
betaN <- extractRawAlleleBFractions(dfN, chromosome=chr);


subplots(2, ncol=1);
plot(beta);
title(sprintf("%s", getName(beta)));
plot(betaN);
title(sprintf("%s (CalMaTe)", getName(betaN)));


###########################################################################
# HISTORY:
# 2010-06-20
# o Created.
###########################################################################
