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
#  \item{data}{An Jx2xI @numeric @array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#  \item{references}{An index @vector in [1,I] or a @logical @vector 
#     of length I specifying which samples are used when calculating the
#     reference signals.  If @NULL, all samples are used.}
#  \item{...}{Additional arguments passed to internal @seemethod "fitCalMaTe.matrix"}
#  \item{truncate}{If @TRUE, final ASCNs are forced to be non-negative
#     while preserving the total CNs.}
#  \item{refAvgFcn}{(optional) A @function that takes a JxI @numeric @matrix
#     an argument \code{na.rm} and returns a @numeric @vector of length J.
#     It should calculate some type of average for each of the J rows, e.g.
#     @see "matrixStats::rowMedians".  
#     If specified, then the total copy numbers of the calibrated ASCNs
#     are standardized toward (twice) the average of the total copy numbers
#     of the calibrated reference ASCNs.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric @array
#   with the same dimension names as argument \code{data}.
# }                                   
#
# @examples "../incl/calmateByThetaAB.Rex"
#
# \seealso{
#  To calibrate (total,fracB) data, 
#  see @seemethod "calmateByTotalAndFracB".
# }
#*/###########################################################################
setMethodS3("calmateByThetaAB", "array", function(data, references=NULL, ..., truncate=FALSE, refAvgFcn=NULL, verbose=FALSE) {
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

  nbrOfSamples <- dim[3];
  if (nbrOfSamples <= 2) {
    throw("Argument 'data' contains less than two samples: ", nbrOfSamples);
  }

  # Argument 'references':
  if (is.null(references)) {
    # The default is that all samples are used to calculate the reference.
    references <- seq(length=nbrOfSamples);
  } else if (is.logical(references)) {
    if (length(references) != nbrOfSamples) {
      throw("Length of argument 'references' does not match the number of samples in argument 'data': ", length(references), " != ", nbrOfSamples);
    }
    references <- which(references);
    if (length(references) == 0) {
      throw("No references samples.");
    }
  } else if (is.numeric(references)) {
    references <- as.integer(references);
    if (any(references < 1 | references > nbrOfSamples)) {
      throw(sprintf("Argument 'references' is out of range [1,%d]", nbrOfSamples));
    }
    if (length(references) == 0) {
      throw("No references samples.");
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);    



  # From here on we force dimension names on the 2nd dimension
  dimnames(data)[[2]] <- c("A", "B");

  
  verbose && enter(verbose, "calmateByThetaAB()");
  verbose && cat(verbose, "ASCN signals:");
  verbose && str(verbose, data);
  verbose && cat(verbose, "Reference samples:");
  verbose && str(verbose, references);

  verbose && enter(verbose, "Identifying non-finite data points");
  # Keep finite values
  ok <- (is.finite(data[,"A",,drop=FALSE]) & is.finite(data[,"B",,drop=FALSE]));
  dim(ok) <- dim(ok)[-2]; # Drop 2nd dimension
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
  verbose && printf(verbose, "Number of SNPs left: ");
  # Drop dimnames for faster processing
  dimnames(dataS) <- NULL;
  for (jj in seq(length=nbrOfSNPs)) {
    if (verbose && (jj %% 100 == 1)) {
      printf(verbose, "%d,", nbrOfSNPs-jj+1, timestamp=FALSE);
    }
    Cjj <- dataS[jj,,,drop=FALSE];  # An 1x2xI array
    dim(Cjj) <- dim(Cjj)[-1]; # A 2xI matrix
    CCjj <- fitCalMaTe(Cjj, references=references, ...);
    # Sanity check
    stopifnot(identical(dim(CCjj), dim(Cjj)));
    dataS[jj,,] <- CCjj;
  } # for (jj ...)
  if (verbose) cat(verbose, "done.");
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

  if (truncate){
    dataC <- truncateThetaAB(dataC);
    verbose && cat(verbose, "Truncated ASCN signals:");
    verbose && str(verbose, dataC);
  } 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Standardize toward a custom average of the references?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(refAvgFcn)) {
    verbose && enter(verbose, "Standardize total copy numbers toward the average reference signals");
    # Extract reference signals
    dataCR <- dataC[,,references,drop=FALSE];
    # Calculate total copy number signals
    yCR <- dataCR[,1,,drop=FALSE]+dataCR[,2,,drop=FALSE];
    dim(yCR) <- dim(yCR)[-2]; # Drop 2nd dimension
    # Calculate the average
    yCR <- refAvgFcn(yCR, na.rm=TRUE);
    # Standardize ASCNs to this average
    dataC <- 2 * dataC / yCR;
    verbose && exit(verbose);
  }

  # Enforce the same dimension names as the input data
  dimnames(dataC) <- dimnames;

  verbose && cat(verbose, "Calibrated (A,B) signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # calmateByThetaAB()


###########################################################################
# HISTORY:
# 2011-03-18 [HB]
# o BUG FIX: calmateByThetaAB() required that the 2nd dimension
#   of argument 'data' had names "A" and "B".
# 2010-08-05 [HB]
# o ROBUSTNESS: Now calmateByThetaAB() asserts that there is at least
#   two samples.
# o BUG FIX: calmateByThetaAB() would not work with only one unit or only
#   one sample.
# 2010-08-02 [HB]
# o Added argument 'refAvgFcn' to calmateByThetaAB().
# 2010-06-19 [HB]
# o Now calmateByThetaAB() uses internal truncateThetaAB() to truncate
#   (CA,CB) values.  Since it operates on arrays, it is much faster.
# 2010-06-18 [HB]
# o Now calmateByThetaAB() calls internal fitCalMaTe().
# o Added argument 'references' to calmateByThetaAB().
# 2010-06-04 [MO]
# o Created.
###########################################################################
