%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CalMaTeNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{CalMaTeNormalization}
\docType{class}
\alias{CalMaTeNormalization}


 \title{The CalMaTeNormalization class}

 \description{
  Package:  calmate \cr
\bold{Class CalMaTeNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{CalMaTeNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{CalMaTeNormalization}\cr
extends \link[R.oo]{Object}\cr



  This class represents the CalMaTe normalization method [1], which
  corrects for SNP effects in allele-specific copy-number estimates
  (ASCNs).
 }

 \usage{CalMaTeNormalization(data=NULL, references=NULL, TRUNCATE=FALSE, method="FracB", tags="*", ...)}

 \arguments{
   \item{data}{A named \code{\link[base]{list}} with data set named \code{"total"} and
     \code{"fracB"} where the former should be of class
     \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}} and the latter of
     class \code{\link[aroma.core]{AromaUnitFracBCnBinarySet}}.  The
     two data sets must be for the same chip type, have the same
     number of samples and the same sample names.}
   \item{references}{A \code{\link[base]{vector}} specifying which samples should be used
     as the reference set.
     By default, all samples are considered.}
   \item{TRUNCATE}{Logical variable to indicate if a truncation is applied.
     By default TRUNCATE=FALSE.}
   \item{method}{Name of the truncation method to use.
     By default method="FracB".}
   \item{tags}{Tags added to the output data sets.}
   \item{...}{Not used.}
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
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \details{
   ...
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
title(sprintf("\%s", getName(beta)));
plot(betaN);
title(sprintf("\%s (CalMaTe)", getName(betaN)));

 }}

 \references{
   [1] ...
 }

 \seealso{
   Low-level versions of the CalMaTe normalization method is available
   via the \code{\link{calmateByThetaAB.array}}() and
   \code{\link{calmateByTotalAndFracB.array}}() methods.
 }
\keyword{classes}
