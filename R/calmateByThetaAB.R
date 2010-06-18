###########################################################################/**
# @set "class=array"
# @RdocMethod calmateByThetaAB
# @alias calmateByThetaAB
# 
# @title "Internal CalMaTe fit function"
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
setMethodS3("calmateByThetaAB", "array", function(data, ..., verbose=FALSE) {
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
    dimnames <- dimnames(dataS);    
  } else {
    verbose && cat(verbose, "All data points are finite.");
    dataS <- data;
  }
  verbose && exit(verbose);

  verbose && enter(verbose, "Fitting CalMaTe");
  dataA <- dataS[,"A",,drop=FALSE];
  dataB <- dataS[,"B",,drop=FALSE];
  dim(dataA) <- dim[-2];
  dim(dataB) <- dim[-2];
  rm(dataS);
  dataT <- calmate(dataA, dataB, ...);
  save(dataT, file="dataT.Rdata");
  verbose && str(verbose, head(dataT));
  rm(dataA, dataB);
  verbose && exit(verbose);

  verbose && enter(verbose, "Restructuring into an array");
  # There is one element per SNP
  dimT <- c(dim(dataT[[1]]), length(dataT));
  dataT <- unlist(dataT, use.names=FALSE);
  dim(dataT) <- dimT;
  dataT <- aperm(dataT, perm=c(3,1,2));
  dimnames(dataT) <- dimnames;
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
}) # calmateByThetaAB()


###########################################################################
# HISTORY:
# 2010-06-04 [MO]
# o Created.
###########################################################################
