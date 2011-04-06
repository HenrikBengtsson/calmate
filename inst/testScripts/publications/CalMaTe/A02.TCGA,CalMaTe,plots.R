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



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose chip type to study
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chipTypes <- c("GenomeWideSNP_6", "Human1M-Duo");
if (interactive() && require("R.menu")) {
  chipType <- textMenu(chipTypes, title="Select chip type:", value=TRUE);
} else {
  chipType <- chipTypes[1];
}


figPath <- file.path("figures", chipType);
figPath <- Arguments$getWritablePath(figPath);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data sets
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (chipType == "Human1M-Duo") {
  dataSet <- "hudsonalpha.org_OV.Human1MDuo.1.1.0";
  tags <- "XY";
  chipType <- "Human1M-Duo"
} else if (chipType == "GenomeWideSNP_6") {
  dataSet <- "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0";
  tags <- "ASCRMAv2";
  chipType <- "GenomeWideSNP_6";
}


if (!exists("dsCList", mode="list")) {
  dsList <- loadSets(dataSet, tags=tags, chipType=chipType, verbose=verbose);
  dsCList <- loadSets(dataSet, tags=c(tags, "CMTN"), chipType=chipType, verbose=verbose);
}
verbose && print(verbose, dsList);
verbose && print(verbose, dsCList);


# Sample of interest
sampleName <- "TCGA-23-1027-01";
verbose && cat(verbose, "Sample name: ", sampleName);


# Plot dimensions etc.
xlab <- "Position (Mb)";
ylab <- list(tcn="Copy number", baf="Allele B Fraction");
ylim <- list(tcn=c(0,6), baf=c(-0.1, 1.1));
at <- list(tcn=c(0, 2, 4, 6), baf=c(0, 1/2, 1));

# Default graphical parameters
setOption("devNew/args/par", list(
  mar=c(2.5,2.5,1.1,1)+0.1, mgp=c(1.4,0.4,0), tcl=-0.3,
  cex=2
));

chromosomes <- 1:23;
for (chr in chromosomes) {
  chrTag <- sprintf("Chr%02d", chr);

  verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", chr, chrTag, length(chromosomes)));
  verbose && cat(verbose, "Sample name: ", sampleName);

  verbose && enter(verbose, "Extracting signals");
  data <- extractSignals(dsList, sampleName=sampleName, chromosome=chr, reference="median", verbose=verbose);
  verbose && print(verbose, data);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting CalMaTe signals");
  dataC <- extractSignals(dsCList, sampleName=sampleName, chromosome=chr, verbose=verbose);
  verbose && print(verbose, dataC);
  verbose && exit(verbose);


  verbose && enter(verbose, "Plotting");
  # Plot dimensions
  pos <- data$tcn$x;
  x <- pos/1e6;
  xlim <- range(x, na.rm=TRUE);

  for (figTag in c("tcn", "baf")) {
    for (cal in c(FALSE, TRUE)) {
      if (cal) {
        dat <- dataC;
        tagsT <- c(tags, "CMTN");
      } else {
        dat <- data;
        tagsT <- tags;
      }

      tags <- c(chrTag, chipType, tagsT, toupper(figTag));
      toPNG(name=sampleName, tags=tags, width=840, aspectRatio=1/3, {
        plot(NA, xlim=xlim, ylim=ylim[[figTag]], 
                 xlab=xlab, ylab=ylab[[figTag]], axes=FALSE);
        axis(side=1);
        axis(side=2, at=at[[figTag]]);
        points(dat[[figTag]], pch=".");
        stext(side=3, pos=0, sprintf("%s", sampleName));
        stext(side=3, pos=1, chrTag);
      });
    } # for (cal ...)
  } # for (figTag ...)

  verbose && exit(verbose);
} # for (chr ...)



###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Created from PN's CalMaTe,Illumina.R script from Nov 2010.
###########################################################################
