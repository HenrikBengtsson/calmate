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



extractSignals <- function(dsList, sampleName, ..., reference=c("auto", "none", "median"), verbose=FALSE) {
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
  
  idxT <- indexOf(dsList$total, sampleName);
  dfT <- getFile(dsList$total, idxT);
  idxB <- indexOf(dsList$fracB, sampleName);
  dfB <- getFile(dsList$fracB, idxB);
  
  verbose && print(verbose, list(tumor=dfT, fracB=dfB));
  verbose && exit(verbose);
  
  verbose && enter(verbose, "Extracting TCNs");
  tcn <- extractRawCopyNumbers(dfT, logBase=NULL, ...);
  verbose && print(verbose, tcn);
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting BAFs");
  baf <- extractRawAlleleBFractions(dfB, ...);
  verbose && print(verbose, baf);
  verbose && exit(verbose);

  if (reference == "auto") {
    verbose && enter(verbose, "Inferring from data set name if CN ratios needs to be calculated");
    hasRatios <- hasTag(dsList$total, "CMTN");
    verbose && cat(verbose, "Has 'CMTN' tag: ", hasRatios);
    reference <- ifelse(hasRatios, "none", "median");
    verbose && cat(verbose, "Inferred argument 'reference': ", reference);
    verbose && exit(verbose);
  }
  
  if (reference == "median") {
    dfTR <- getAverageFile(dsList$total, verbose=less(verbose,5));
    
    verbose && enter(verbose, "Extracting TCNs for reference pool");
    tcnR <- extractRawCopyNumbers(dfTR, logBase=NULL, ...);
    verbose && print(verbose, tcnR);
    verbose && exit(verbose);
    
    tcn <- divideBy(tcn, tcnR);
    tcn$y <- 2*tcn$y;
  }
  
  res <- list(tcn=tcn, baf=baf);
  verbose && exit(verbose);
  
  res;
} # extractSignals()


plotBvsB <- function(beta, pch=19, Blim=c(-0.2,+1.2), xlab="Normal BAF", ylab="Tumor BAF", dataSet=NULL, chipType=NULL, sampleName=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Blim, ylim=Blim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=c(0,1/2,1), label=c("0","1/2","1"));
  }
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  points(beta, pch=pch, cex=1);

  if (anTags == "dens") {
    for (dd in 1:2) {
      d <- density(beta[,dd], adjust=0.4, from=Blim[1], to=Blim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotBvsB()



plotASCN <- function(ascn, pch=19, Clim=c(-0.4,4), xlab=expression(C[A]), ylab=expression(C[B]), dataSet=NULL, chipType=NULL, sampleName=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=Clim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  for (ss in 1:2) {
    axis(side=ss, at=0:round(max(Clim)+1L));
  }
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  points(ascn, pch=pch, cex=1);

  if (anTags == "dens") {
    for (dd in 1:2) {
      d <- density(ascn[,dd], adjust=0.4, from=Clim[1], to=Clim[2], na.rm=TRUE);
      draw(d, lwd=3, col="#666666", side=dd, 
                       height=0.05, scale="relative", xpd=FALSE);
    }
  }
} # plotASCN()



plotC1C2 <- function(ascn, x, muN, pch=19, xlim=NULL, Clim=c(-0.4,4), xlab="Position (Mb)", ylab=expression(list(C[1],C[2])), dataSet=NULL, chipType=NULL, sampleName=NULL, segName=NULL, segTag=NULL, tagsT=NULL, ...) {
  if (is.null(xlim)) {
    xlim <- range(x, na.rm=TRUE);
    dx <- diff(xlim);
    xlim[1] <- xlim[1] - 0.05*dx;
    xlim[2] <- xlim[2] + 0.05*dx;
  }

#  par(mar=c(4,5,3,3)+0.1, cex=0.8, cex.lab=2.4, cex.axis=2.2);
  par(mar=c(5.8,7,3,3)+0.1, mgp=c(4,1.4,0), cex=0.8, cex.lab=3.8, cex.axis=2.8);
  plot(NA, xlim=xlim, ylim=Clim, xlab=xlab, ylab=ylab, axes=FALSE);
  axis(side=1);
  axis(side=2, at=0:round(max(Clim)+1L));
  box();
  abline(a=0, b=1, lty=3, lwd=2);

  # Text annotations
  if (!is.null(segName) & !is.null(segTag)) {
    stext(side=3, pos=0, sprintf("%s @ %s", segName, segTag), cex=2.0);
  }
  if (!is.null(tagsT)) {
    tagsT <- fullname(tagsT);
    stext(side=3, pos=1, tagsT, cex=2.0);
  }
  if (!is.null(dataSet) & !is.null(chipType)) {
#    stext(side=4, pos=0, sprintf("%s (%s)", dataSet, chipType), cex=1.5);
  }
  if (!is.null(sampleName)) {
    stext(side=4, pos=1, sampleName,  cex=1.5);
  }

  hets <- which(muN == 1/2);
  xH <- x[hets];
  c1c2 <- ascn[hets,];
  rr <- which(c1c2[,2] < c1c2[,1]);
  c1c2[rr,] <- c1c2[rr,2:1];

  c1 <- c1c2[,1];
  c2 <- c1c2[,2];

  points(xH, c1, pch=pch, cex=1, col="blue");
  points(xH, c2, pch=pch, cex=1, col="red");

  abline(h=median(c1,na.rm=TRUE), lwd=6, col="blue");
  abline(h=median(c2,na.rm=TRUE), lwd=6, col="red");

  if (anTags == "dens") {
    d <- density(c1, adjust=0.8, na.rm=TRUE);
    draw(d, lwd=3, col="blue", side=2, height=0.05, scale="relative", xpd=FALSE);

    d <- density(c2, adjust=0.8, na.rm=TRUE);
    draw(d, lwd=3, col="red", side=2, height=0.05, scale="relative", xpd=FALSE);
  }

} # plotC1C2()


###########################################################################
# HISTORY:
# 2011-03-09 [HB]
# o Added option "auto" for argument 'reference' of extractSignals().
# o Created.
###########################################################################
