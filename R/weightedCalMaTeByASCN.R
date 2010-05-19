# @title "Internal CalMaTe fit function"
#
# \description{
#  @get "title".
# }
#
# @usage
#
# \arguments{
#  \item{data}{An Jx2xI @numeric array, where J is the number of SNPs,
#          2 is the number of alleles, and I is the number of samples.}
#  \item{...}{Additional arguments passed to internal 
#          \code{CalMaTeWeighted().}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric array.
# }
#
setMethodS3("weightedCalMaTeByASCN", "array", function(data, ..., verbose=FALSE) {
  # Argument 'data':
  if (!is.array(data)) {
    throw("Argument 'data' is not an array: ", class(data)[1]);
  }
  dim <- dim(data);
  if (length(dim) != 3) {
    throw("Argument 'data' is not a 3-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (dim[2] != 2) {
    throw("Argument 'data' is not a Jx2xI-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);

  verbose && enter(verbose, "weightedCalMaTeByASCN()");
  verbose && cat(verbose, "ASCN signals:");
  verbose && str(verbose, data);

  dim <- dim(data);
  dimnames <- dimnames(data);

  verbose && enter(verbose, "Identifying non-finite data points");
  # Keep finite values
  ok <- (is.finite(data[,1,]) & is.finite(data[,2,]));
  ok <- rowAlls(ok);
  verbose && summary(verbose, ok);
  hasNonFinite <- any(!ok);
  if (hasNonFinite) {
    verbose && enter(verbose, "Excluding non-finite data points");
    dataS <- dataS[ok,,,drop=FALSE];
    verbose && str(verbose, data);
    verbose && exit(verbose);
  } else {
    verbose && cat(verbose, "All data points are finite.");
    dataS <- data;
  }
  verbose && exit(verbose);
  verbose && enter(verbose, "Fitting CalMaTe");
  dataA <- dataS[,1,,drop=FALSE];
  dataB <- dataS[,2,,drop=FALSE];
  dim(dataA) <- dim[-2];
  dim(dataB) <- dim[-2];
  rm(dataS);
  dataT <- weightedCalMaTe(dataA, dataB, ...);
  verbose && str(verbose, head(dataT));
  rm(dataA, dataB);
  verbose && exit(verbose);

  verbose && enter(verbose, "Restructuring into an array");
  # There is one element per SNP
  dimT <- c(dim(dataT[[1]]), length(dataT));
  dataT <- unlist(dataT, use.names=FALSE);
  dim(dataT) <- dimT;
  dataT <- aperm(dataT, perm=c(3,1,2));
  verbose && str(verbose, dataT);
  verbose && exit(verbose);

  if (hasNonFinite) {
    verbose && enter(verbose, "Expanding to array with non-finite");
    dataC <- data;
    dataC[ok,,] <- dataT;
    verbose && str(verbose, dataC);
    verbose && exit(verbose);
  } else {
    dataC <- dataT;
  }
  rm(dataT);

  # Sanity check
  stopifnot(identical(dim(dataC), dim(data)));

  verbose && cat(verbose, "Calibrated ASCN signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # weightedCalMaTeByASCN()


###########################################################################
# HISTORY:
# 2010-05-18
# o Created.
###########################################################################
