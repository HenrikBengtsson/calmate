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
sourceTo("001.include.R", path=path);
sourceTo("002.datasets.R", path=path);


figPath <- file.path("figures", chipType);
figPath <- Arguments$getWritablePath(figPath);
  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calibrated or not?
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
calTag <- textMenu(calTags, title="Choose calibration method:", value=TRUE);
if (calTag == "<none>") calTag <- NULL;


dsList <- loadSets(dataSet, tags=c(tags, calTag), chipType=chipType, verbose=verbose);
verbose && print(verbose, dsList);

# Sanity check
# stopifnot(all(sapply(dsList, FUN=length) == 40));


dsTags <- getTags(dsList$total, collapse=",");
dsTags <- gsub("ACC,-XY,BPN,-XY,RMA,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("ACC,-XY,BPN,-XY,AVG,FLN,-XY", "ASCRMAv2", dsTags);
dsTags <- gsub("CMTN", "CalMaTe", dsTags);
  
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extracting signals
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
units <- 1000+1:100;
dd <- extractCACB(dsList, units=units);

# Drop 'refs=N' tag
tagsT <- dropTags(dsTags, drop="refs=N");

for (uu in seq(along=units)) {
  unit <- units[uu];
  unitTag <- sprintf("unit%06d", unit);
  unitName <- unit;

  ddT <- dd[uu,,];

  fullname <- fullname(sampleName, tagsT, unitTag, "CACB");
  devEval("png", name=fullname, width=840, aspectRatio=1, {
    plotMultiArrayCACB(ddT, unitName=unitName, tagsT=tagsT);
  }, path=figPath);
} # for (uu ...)


###########################################################################
# HISTORY:
# 2011-03-10 [HB]
# o Created.
###########################################################################
