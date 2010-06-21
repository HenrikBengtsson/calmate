
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CRMAv2 - Preprocess raw Affymetrix data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("aroma.affymetrix");  # Needed for CRMAv2
dataSet <- "Affymetrix_2006-TumorNormal";
chipType <- "Mapping250K_Nsp";
dsList <- doCRMAv2(dataSet, chipType=chipType, combineAlleles=FALSE, 
                                             plm="RmaCnPlm", verbose=-10);
print(dsList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# CalMaTe - Post-normalization of ASCNs estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
asn <- CalMaTeNormalization(dsList);
print(asn);

# For speed issues, we will here only process loci on Chromosome 17.
chr <- 17;
ugp <- getAromaUgpFile(dsList$total);
units <- getUnitsOnChromosome(ugp, chr);

dsNList <- process(asn, units=units, verbose=verbose);
print(dsNList);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Plot allele B fractions (before and after)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Sample #1 and Chromosome 17
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
