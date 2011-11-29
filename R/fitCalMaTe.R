###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitCalMaTe
# @alias fitCalMaTe
#
# @title "Calibrates SNP loci according to the CalMaTe method"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{dataT}{A 2xI @numeric @matrix of allele specific copy numbers (ASCNs), 
#     where 2 is the number alleles and I is the number of samples.}
#  \item{references}{A @logical or @numeric @vector specifying which
#     samples should be used as the reference set.}
#  \item{fB1,fB2}{Thresholds for calling genotypes AA, AB, BB from the
#     allele B fractions.}
#  \item{maxIter}{The maximum number of iterations without converging
#     before the algorithm quits.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a 2xI @numeric @matrix of calibrated ASCNs.
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fitCalMaTe", "matrix", function(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  nbrOfSNPs <- nrow(dataT);
  nbrOfReferences <- length(references);
  
  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  dataT[dataT < eps] <- eps;
  
  a <- max(max(dataT[2,] / (pmax(dataT[1,],0) + 1e-4)), max(dataT[1,] / (pmax(dataT[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  dataT <- Giro %*% dataT;

  # Extract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking if all the samples are homozygous
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Twist half of the samples in case there is only one allele?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  onlyOneAllele <- (abs(sum(naiveGenoDiff)/2) == length(naiveGenoDiff));
  if (onlyOneAllele) {
    idxs <- seq(length=ncol(dataT)/2);
    dataT[1:2,idxs] <- dataT[2:1,idxs, drop=FALSE];

    # Update precalcalculated signals
    TR <- dataT[,references, drop=FALSE];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy numbers must be close to 2 for the reference samples or
  # (if there are not control samples) for most of the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  H <- matrix(2, nrow=nbrOfReferences, ncol=1, byrow=FALSE);
  fit <- rlm(t(TR), H, maxit=maxIter);
  matSum <- fit$coefficients;
  coeffs <- fit$w;
  dataT <- diag(matSum) %*% dataT;

  # Reextract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  fit <- rlm(t(TR), naiveGenoDiff, maxit=maxIter, weights=coeffs);
  matDiff <- fit$coefficients;

  # T matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ T ]   [MatDiff[1] MatDiff[2]]
  U <- matrix(c(0.5, 0.5, 0.5, -0.5), nrow=2, ncol=2, byrow=FALSE);
  V <- matrix(c(c(1,1), matDiff), nrow=2, ncol=2, byrow=TRUE);
  T <- U %*% V;
  
  res <- T %*% dataT;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Undo the previous change applied to the data in case there is 
  # only one allele    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (onlyOneAllele) {
    idxs <- seq(length=ncol(res)/2);
    res[1:2,idxs] <- res[2:1,idxs, drop=FALSE];
  }

  # Return parameter estimates(?)
  ## attr(res, "modelFit") <- list(fit=fit);

  res;
}, protected=TRUE) # fitCalMaTe()


###########################################################################
# HISTORY:
# 2011-11-29 [MO]
# o Change matrix "T" by "dataT" and "P" by "T"
# 2010-08-02 [HB]
# o ROBUSTNESS: Now fitCalMaTe() also works (technically) when there is
#   only one reference.
# o Made into an S3 method for matrix:es.
# 2010-06-18 [HB]
# o Created from refineCN.list().
###########################################################################
