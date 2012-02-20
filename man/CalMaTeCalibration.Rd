%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CalMaTeCalibration.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CalMaTeCalibration}
\docType{class}
\alias{CalMaTeCalibration}


 \title{The CalMaTeCalibration class}

 \description{
  Package:  calmate \cr
\bold{Class CalMaTeCalibration}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{CalMaTeCalibration}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CalMaTeCalibration}\cr
extends \link[R.oo]{Object}\cr



  This class represents the CalMaTe method [1], which
  corrects for SNP effects in allele-specific copy-number estimates
  (ASCNs).
 }

 \usage{CalMaTeCalibration(data=NULL, tags="*", references=NULL, flavor=c("v2", "v1"), ...)}

 \arguments{
   \item{data}{A named \code{\link[base]{list}} with data set named \code{"total"} and
     \code{"fracB"} where the former should be of class
     \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}} and the latter of
     class \code{\link[aroma.core]{AromaUnitFracBCnBinarySet}}.  The
     two data sets must be for the same chip type, have the same
     number of samples and the same sample names.}
   \item{tags}{Tags added to the output data sets.}
   \item{references}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying which samples
     should be as reference samples for estimating the model parameters.
     If \code{\link[base]{NULL}}, all samples are used.}
   \item{flavor}{A \code{\link[base]{character}} string specifying which flavor of the
     CalMaTe algorithm to use for fitting the model.
     See \code{\link{fitCalMaTeInternal}} for details.}
   \item{...}{Additional arguments passed to \code{\link{calmateByTotalAndFracB}}().}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findUnitsTodo} \tab  -\cr
 \tab \code{getDataSets} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getOutputDataSets} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getReferences} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \section{Reference samples}{
  In order to estimate the calibration parameters, the model assumes
  that, for any given SNP, there are a majority of samples that are
  diploid at that SNP.  Note that it does not have to be the same set
  of samples for all SNPs.

  By using argument \code{references}, it is possible so specify which
  samples should be used when estimating the calibration parameters.
  This is useful when for instance there are several tumor samples with
  unknown properties as well as a set of normal samples that can be
  assumed to be diploid.

  Theoretical, a minimum of three reference samples are needed in order
  for the model to be identifiable.  If less, an error is thrown.
  However, in practice more reference samples should be used, that is,
  in the order of at least 6-10 reference samples with a diverse set
  of genotypes.
 }

 \examples{\dontrun{
  
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
# CalMaTe - Post-calibration of ASCNs estimates
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
asn <- CalMaTeCalibration(dsList);
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

# Extract raw (TCN,BAF)
df <- getFile(dsList$total, ii);
dfR <- getAverageFile(dsList$total, verbose=verbose);
gamma <- extractRawCopyNumbers(df, logBase=NULL, chromosome=chr);
gammaR <- extractRawCopyNumbers(dfR, logBase=NULL, chromosome=chr);
gamma <- 2*divideBy(gamma, gammaR);
df <- getFile(dsList$fracB, ii);
beta <- extractRawAlleleBFractions(df, chromosome=chr);

# Extract calibrated (TCN,BAF)
dfN <- getFile(dsNList$fracB, ii);
betaN <- extractRawAlleleBFractions(dfN, chromosome=chr);
dfN <- getFile(dsNList$total, ii);
gammaN <- extractRawCopyNumbers(dfN, logBase=NULL, chromosome=chr);

# Plot
subplots(4, ncol=2, byrow=FALSE);
plot(beta);
title(sprintf("\%s", getName(beta)));
plot(gamma);
plot(betaN);
title(sprintf("\%s (CalMaTe)", getName(betaN)));
plot(gammaN);

 }}

 \references{
   [1] ...
 }

 \seealso{
   Low-level versions of the CalMaTe method is available
   via the \code{\link{calmateByThetaAB}}() and
   \code{\link{calmateByTotalAndFracB}}() methods.

  For further information on the internal fit functions, see
  \code{\link{fitCalMaTeInternal}}.
 }
\keyword{classes}
