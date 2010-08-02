###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitCalMaTeCNprobes
# 
# @title "Normalizes non-polymorphic copy number loci according to the CalMaTe method"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{T}{A JxI @numeric @matrix, where J is the number of loci
#                      and I is the number of samples.}
#  \item{references}{A @logical or @numeric @vector specifying which
#     samples should be used as the reference set.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @numeric @vector of length J.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fitCalMaTeCNprobes", "matrix", function(T, references, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.

  # Extract the reference samples
  Tref <- T[,references, drop=FALSE];
  
  TR <- rowMedians(Tref, na.rm=TRUE);

  res <- 2 * T/TR;
  
  #nbrOfReferences <- length(references);
#  H <- matrix(2, nrow=nbrOfReferences, ncol=1,byrow=FALSE);
#  matSum <- matrix(0, nrow=1, ncol=nrow(T));
#
#  ind <- which(rowSums(T)!=0);
#  fit <- matrix(0, nrow=1, ncol=length(ind));
#
#  for(ii in ind){
#    fit[ii] <- rlm(T[ii,references], H, maxit=50);
#  }
#  
#  if(length(ind)!=0 && ind != 0)
#    matSum[ind] <- fit$coefficients;
#  res <- matSum %*% T;
#
  res;
}, protected=TRUE) # fitCalMaTeCNprobes()


###########################################################################
# HISTORY:
# 2010-08-02 [HB]
# o ROBUSTNESS: fitCalMaTeCNprobes() can now also handle missing values.
# o Made into an S3 method for matrix:es.
# 2010-06-22 [MO]
# o Created.
###########################################################################
