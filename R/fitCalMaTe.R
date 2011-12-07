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
  
  # Argument "references"
  if(nbrOfReferences < 3){
    throw("At least 3 reference samples or if it is null, at least 3 samples in the dataset");
  }

  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  dataT[dataT < eps] <- eps;

  a <- max(max(dataT[2,] / (pmax(dataT[1,],0) + 1e-4)), max(dataT[1,] / (pmax(dataT[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  dataT <- Giro %*% dataT;

  # Extract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # Set some weights based on the median
  S <- colSums(TR);
  S[S<0] <- 0;
  CN <- 2* S / median(S);
  w1 <- 0.1  + 0.9 * sqrt(sqrt(2*(1- pnorm(abs(CN -2 ) / median(abs((CN -2)))))))

  if (sum(is.nan(w1)) > 0){
    w1 = rep(1, length(w1));
  }


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
    idxs <- references[seq(length=ncol(TR)/2)];  # Changed from previous version
    dataT[1:2,idxs] <- dataT[2:1,idxs, drop=FALSE];

    # Update precalcalculated signals
    TR <- dataT[,references, drop=FALSE];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy numbers must be close to 2 for the reference samples or
  # (if there are not control samples) for most of the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  H <- matrix(2, nrow=nbrOfReferences, ncol=1, byrow=FALSE);
  fit <- rlm(H ~0 + t(TR), maxit=maxIter, weights=w1); # alternatively can be set method = "MM"
  matSum <- fit$coefficients;
  coeffs <- fit$w;
  coeffs[coeffs<0] <- 1e-8;
  coeffs <- coeffs * w1;
  dataT <- diag(matSum) %*% dataT;

  # Reextract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  fit <- rlm(naiveGenoDiff ~ 0 + t(TR), maxit=maxIter, weights=coeffs);
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
# 2011-12-07 [MO]
# o At least 3 reference samples.
# 2011-04-12 [AR]
#   · Bug fixed: there was a bug for SNPs with a single allele
#     when using a set of references
#   · Set some initial weights based on the median that improves
#     the breakdown point if no reference samples are provided
#     when using a set of references
# 2011-11-29 [MO]
# o Change matrix "T" by "dataT" and "P" by "T"
# 2010-08-02 [HB]
# o ROBUSTNESS: Now fitCalMaTe() also works (technically) when there is
#   only one reference.
# o Made into an S3 method for matrix:es.
# 2010-06-18 [HB]
# o Created from refineCN.list().
###########################################################################
