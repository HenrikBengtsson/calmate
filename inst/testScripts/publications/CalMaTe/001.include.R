###########################################################################
# 001.include.R
#
# Author: Henrik Bengtsson
# 
# Description:
# This script defines functions and other objects that are needed by
# all scripts in this directory.
###########################################################################
library("aroma.cn");
library("calmate");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
loadSets <- function(..., verbose=FALSE) {
  dsT <- AromaUnitTotalCnBinarySet$byName(..., verbose=verbose);
  dsB <- AromaUnitFracBCnBinarySet$byName(..., verbose=verbose);
  list(total=dsT, fracB=dsB);
} # loadSets()



extractSignals <- function(dsList, sampleName, chromosome, reference=c("none", "median"), verbose=FALSE) {
  # Argument 'reference':
  reference <- match.arg(reference);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting total CN signals");
  
  verbose && enter(verbose, "Extracting sample of interest");
  verbose && cat(verbose, "Sample name: ", sampleName);
  idx <- indexOf(dsList$total, sampleName);
  verbose && cat(verbose, "Index: ", idx);
  idx <- Arguments$getIndex(idx);
  
  dfT <- getFile(dsList$total, idx);
  dfB <- getFile(dsList$fracB, idx);
  
  verbose && print(verbose, list(tumor=dfT, fracB=dfB));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Extracting TCNs");
  tcn <- extractRawCopyNumbers(dfT, logBase=NULL, chromosome=chromosome);
  verbose && print(verbose, tcn);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting BAFs");
  baf <- extractRawAlleleBFractions(dfB, chromosome=chromosome);
  verbose && print(verbose, baf);
  verbose && exit(verbose);
  
  if (reference == "median") {
    dfTR <- getAverageFile(dsList$total, verbose=less(verbose,5));
    
    verbose && enter(verbose, "Extracting TCNs for reference pool");
    tcnR <- extractRawCopyNumbers(dfTR, logBase=NULL, chromosome=chromosome);
    verbose && print(verbose, tcnR);
    verbose && exit(verbose);
    
    tcn <- divideBy(tcn, tcnR);
    tcn$y <- 2*tcn$y;
  }
  
  res <- list(tcn=tcn, baf=baf);
  verbose && exit(verbose);
  
  res;
} # extractSignals()


###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Created.
###########################################################################
