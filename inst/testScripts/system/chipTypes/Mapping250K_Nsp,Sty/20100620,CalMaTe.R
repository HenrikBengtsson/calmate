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
asn <- CalMaTeNormalization(list(total=dsT, fracB=dsB));
print(asn);


dsNList <- process(asn, verbose=verbose);
print(dsNList);


###########################################################################
# HISTORY:
# 2010-06-20
# o Created.
###########################################################################
