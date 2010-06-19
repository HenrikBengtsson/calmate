###########################################################################/**
# @set "class=array"
# @RdocMethod calmateByThetaAB
# @alias calmateByThetaAB
# 
# @title "Normalize allele-specific copy numbers (CA,CB)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#  \item{references}{A @logical @vector of length I, or an index @vector 
#          with values in [0,I]. If @NULL, all samples are considered to
#          be reference samples (==1:I).}
#  \item{...}{Additional arguments passed to internal 
#             \code{calmate().}}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric array.
# }                                   
#
# @examples "../incl/calmateByThetaAB.Rex"
#
# \seealso{
#  To calibrate (total,fracB) data, 
#  see @seemethod "calmateByTotalAndFracB".
# }
#*/###########################################################################
setMethodS3("calmateByThetaAB", "array", function(data, references=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.array(data)) {
    throw("Argument 'data' is not an array: ", class(data)[1]);
  }
  dim <- dim(data);
  dimnames <- dimnames(data);
  if (length(dim) != 3) {
    throw("Argument 'data' is not a 3-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (dim[2] != 2) {
    throw("Argument 'data' is not a Jx2xI-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (!is.null(dimnames[[2]])) {
    if (!identical(dimnames[[2]], c("A", "B"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'A' & 'B': ", paste(dimnames[[2]], collapse=", "));
    }
  }

  # Argument 'references':
  if (is.null(references)) {
    # The default is that all samples are used to calculate the reference.
    references <- rep(TRUE, times=dim[3]);
  } else if (is.logical(references)) {
    if (length(references) != dim[3]) {
      throw("Length of argument 'references' does not match the number of samples in argument 'data': ", length(references), " != ", dim[3]);
    }
  } else if (is.numeric(references)) {
    references <- as.integer(references);
    if (any(references < 1 | references > dim[3])) {
      throw(sprintf("Argument 'references' is out of range [1,%d]", dim[3]));
    }
    idxs <- references;
    references <- rep(FALSE, times=dim[3]);
    references[idxs] <- TRUE;
  }



  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);    
  
  verbose && enter(verbose, "calmateByThetaAB()");
  verbose && cat(verbose, "ASCN signals:");
  verbose && str(verbose, data);


  verbose && enter(verbose, "Identifying non-finite data points");
  # Keep finite values
  ok <- (is.finite(data[,"A",]) & is.finite(data[,"B",]));
  ok <- rowAlls(ok);
  verbose && summary(verbose, ok);
  hasNonFinite <- any(!ok);
  if (hasNonFinite) {
    verbose && enter(verbose, "Excluding non-finite data points");
    dataS <- data[ok,,,drop=FALSE];
    verbose && str(verbose, data);
    verbose && exit(verbose);
    dim <- dim(dataS);
  } else {
    verbose && cat(verbose, "All data points are finite.");
    dataS <- data;
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Fitting CalMaTe");
  nbrOfSNPs <- dim(dataS)[1];
  verbose && cat(verbose, "Number of SNPs: ", nbrOfSNPs);
  verbose && printf(verbose, "Number of SNPs left:");
  # Drop dimnames for faster processing
  dimnames(dataS) <- NULL;
  for (jj in seq(length=nbrOfSNPs)) {
    if (verbose && (jj %% 100 == 1)) printf(verbose, "%d,", nbrOfSNPs-jj+1);
    Cjj <- dataS[jj,,,drop=TRUE];  # An 2xI matrix
    CCjj <- fitCalMaTe(Cjj, refs=references, ...);
    # Sanity check
    stopifnot(identical(dim(CCjj), dim(Cjj)));
    dataS[jj,,] <- CCjj;
  } # for (jj ...)
  if (verbose) cat(verbose, "done.");
#  verbose && str(verbose, dataS);
  verbose && exit(verbose);

  if (hasNonFinite) {
    verbose && enter(verbose, "Expanding to array with non-finite");
    dataC <- data;
    dataC[ok,,] <- dataS;
    verbose && str(verbose, dataC);
    verbose && exit(verbose);
  } else {
    dataC <- dataS;
    dimnames(dataC) <- dimnames(data);
  }
  rm(dataS);

  # Sanity check
  stopifnot(identical(dim(dataC), dim(data)));

  verbose && cat(verbose, "Calibrated ASCN signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # calmateByThetaAB()


###########################################################################
# HISTORY:
# 2010-06-18 [HB]
# o Now calmateByThetaAB() calls internal refineCN2().
# o Added argument 'references' to calmateByThetaAB().
# 2010-06-04 [MO]
# o Created.
###########################################################################
