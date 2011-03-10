###########################################################################
# Title:
# Author: Henrik Bengtsson
###########################################################################


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Loading support files
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find the pathname and directory of this script
library("R.utils");
pathname <- names(findSourceTraceback())[1];
path <- dirname(pathname);

# Loading include files
sourceTo(file.path(path, "001.include.R"));
sourceTo(file.path(path, "002.datasets.R"));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose chip type to study
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipTypes <- c("Mapping250K_Nsp", "GenomeWideSNP_6", "Human1M-Duo");
if (interactive() && require("R.menu")) {
  chipType <- textMenu(chipTypes, title="Select chip type:", value=TRUE);
} else {
  chipType <- chipTypes[1];
}


figPath <- file.path("figures", chipType);
figPath <- Arguments$getWritablePath(figPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


segments <- NULL;
if (chipType == "Mapping250K_Nsp") {
  ds <- datasets$geo;
  dataSet <- "GSE12702";
  tags <- "ACC,-XY,BPN,-XY,RMA,FLN,-XY";
  chipType <- "Mapping250K_Nsp";
  
  # Sample of interest
  sampleName <- "GSM318736";
} else if (chipType == "Human1M-Duo") {
  ds <- datasets$tcga;
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.1.1.0";
  tags <- "XY";
  chipType <- "Human1M-Duo";
  sampleName <- "TCGA-23-1027"; 
} else if (chipType == "GenomeWideSNP_6") {
  ds <- datasets$tcga;
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0";
  tags <- "ASCRMAv2";
  chipType <- "GenomeWideSNP_6";
  sampleName <- "TCGA-23-1027"; 
} 

segments <- ds$segments;
getPairs <- ds$getPairs;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrated or not?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
calTags <- c("<none>", "CMTN", "CMTN,refs=N");
calTag <- textMenu(calTags, title="Choose calibration method:", value=TRUE);
if (calTag == "<none>") calTag <- NULL;

anTags <- c("dens", "<none>");
anTags <- anTags[1];
if (length(anTags) > 1) {
  anTags <- textMenu(anTags, title="Choose graphical annotation:", value=TRUE);
  if (anTags == "<none>") anTags <- NULL;
}



dsList <- loadSets(dataSet, tags=c(tags, calTag), chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);

# Sanity check
# stopifnot(all(sapply(dsList, FUN=length) == 40));


dsTags <- getTags(dsList$total, collapse=",");
dsTags <- gsub("ACC,-XY,BPN,-XY,RMA,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("ACC,-XY,BPN,-XY,AVG,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("CMTN", "CalMaTe", dsTags);
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## refs <- indexOf(dsList$total, pairs[,"normal"]);
##  # Sanity check
##  stopifnot(length(idxs) == 20);
## verbose && cat(verbose, "Number of reference samples for CalMaTe: ", length(refs));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extracting one tumor-normal pair
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Extracting tumor-normal pair of interest");

# Identify tumor of interest
pairs <- getPairs(dsList);

# Identify pair
pair <- pairs[sampleName,,drop=TRUE];
verbose && cat(verbose, "Tumor: ", pair["tumor"]);
verbose && cat(verbose, "Normal: ", pair["normal"]);
pairName <- paste(pair, collapse="vs"); 
verbose && cat(verbose, "Pair: ", pairName);

verbose && exit(verbose);
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segments to be studied
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Blim <- c(-0.3,1.3);
xlab <- "Normal BAF";
ylab <- "Tumor BAF";
pch <- 19;

for (kk in seq(along=segments)) {
  segName <- names(segments)[kk];
  segment <- segments[[segName]];
  chr <- segment[1];
  region <- segment[2:3];

  segTag <- sprintf("Chr%02d:%g-%gMb", chr, region[1], region[2]);
  
  verbose && enter(verbose, sprintf("Segment #%d ('%s') of %d", kk, segTag, length(segments)));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dataList <- lapply(pair, FUN=function(name) {
    extractSignals(dsList, sampleName=name, chromosome=chr, region=region*1e6, verbose=verbose);
  });
  names(dataList) <- names(pair);

  # Extract (gammaT, gammaN)
  gammaT <- getSignals(dataList$tumor$tcn);
  gammaN <- getSignals(dataList$normal$tcn);
  
  # Extract (betaT, betaN)
  betaT <- getSignals(dataList$tumor$baf);
  betaN <- getSignals(dataList$normal$baf);
  muN <- callNaiveGenotypes(betaN);
  x <- dataList$normal$tcn$x / 1e6;
  xH <- x[muN == 1/2];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Plot (betaN,betaT)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  for (tbn in c(FALSE, TRUE)) {
    beta <- cbind(N=betaN, T=betaT);

    # Drop 'refs=N' tag
    tagsT <- dropTags(dsTags, drop="refs=N");

    if (tbn) {
      verbose && enter(verbose, "TumorBoosting");
      betaTN <- normalizeTumorBoost(betaT=betaT, betaN=betaN);
      verbose && str(verbose, betaTN);
      verbose && exit(verbose);
      beta[,"T"] <- betaTN;
      tagsT <- fullname(tagsT, "TBN");
    }

    ascn <- cbind(A=(1-beta[,"T"])*gammaT, B=beta[,"T"]*gammaT);

    fullname <- fullname(sampleName, gsub(":","_",segTag), tagsT, "BvsB", anTags);
    devEval("png", name=fullname, width=840, aspectRatio=1, {
      plotBvsB(beta, sampleName=sampleName, segName=segName, segTag=segTag, 
                          dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    }, path=figPath);

    fullname <- fullname(sampleName, gsub(":","_",segTag), tagsT, "ASCN", anTags);
    devEval("png", name=fullname, width=840, aspectRatio=1, {
      plotASCN(ascn, sampleName=sampleName, segName=segName, segTag=segTag, 
                          dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    }, path=figPath);

    fullname <- fullname(sampleName, gsub(":","_",segTag), tagsT, "C1C2", anTags);
    devEval("png", name=fullname, width=840, aspectRatio=0.7, {
      plotC1C2(ascn, x=x, muN=muN, sampleName=sampleName, segName=segName, segTag=segTag, dataSet=dataSet, tagsT=tagsT, chipType=chipType);
    }, path=figPath);
  } # for (tbn ...)

  verbose && exit(verbose);
} # for (kk ...)



###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Created from PN's code in the online vignette.
###########################################################################
