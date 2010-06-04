###########################################################################/**
# @set "class=array"
# @RdocMethod calmateByTotalAndFracB
# @alias calmateByTotalAndFracB
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
#                      2 is total and fracB, and I is the number of samples.}
#  \item{...}{Additional arguments passed to 
#         @seemethod "calmateByThetaAB".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric array.
# }
#
# @examples "../incl/calmateByTotalAndFracB.Rex"
#
# \seealso{
#  To calibrate (thetaA,thetaB) data, 
#  see @seemethod "calmateByThetaAB".
# }
#*/###########################################################################
setMethodS3("calmateByTotalAndFracB", "array", function(data, ..., verbose=FALSE) {
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
    if (!identical(dimnames[[2]], c("total", "fracB"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'total' & 'fracB': ", paste(dimnames[[2]], collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "calmateByTotalAndFracB()");
  verbose && cat(verbose, "(total,fracB) signals:");
  verbose && str(verbose, data);

  verbose && enter(verbose, "Transforming to (thetaA, thetaB)");
  data <- totalAndFracB2ThetaAB(data, verbose=less(verbose, 5));
  verbose && str(verbose, data);
  verbose && exit(verbose);

  dataC <- calmateByThetaAB(data, ..., verbose=verbose);  

  verbose && enter(verbose, "Backtransforming to (total, fracB)");
  dataC <- thetaAB2TotalAndFracB(dataC, verbose=less(verbose, 5));
  verbose && str(verbose, dataC);
  verbose && exit(verbose);

  # Sanity check
  stopifnot(identical(dim(dataC), dim(data)));

  verbose && cat(verbose, "Calibrated (total,fracB) signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # calmateByTotalAndFracB()


###########################################################################
# HISTORY:
# 2010-06-04 [MO]
# o Created.
###########################################################################
