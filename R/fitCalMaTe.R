###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitCalMaTe
# @alias fitCalMaTe
#
# @title "Calibrates SNP loci according to the CalMaTe method"
#
# \description{
#  @get "title".
#  \emph{Note: This is an internal function of the package. It should not
#   be used elsewhere. Its API may change at any time.}
# }
#
# @synopsis
#
# \arguments{
#  \item{dataT}{A 2xI @numeric @matrix of allele specific copy numbers (ASCNs),
#     where 2 is the number alleles and I is the number of samples.}
#  \item{references}{A @integer @vector with elements in [1,I] specifying
#     which samples should be used as the reference set.}
#  \item{fB1, fB2}{Thresholds for calling genotypes AA, AB, BB from the
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

  eps2 <- 1e-4;
  a <- max(max(dataT[2,] / (pmax(dataT[1,],0) + eps2)), max(dataT[1,] / (pmax(dataT[2,],0) + eps2)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  dataT <- Giro %*% dataT;

  # Extract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # Set some weights based on the median
  S <- colSums(TR);
  S[S < 0] <- 0;
  CN <- 2*S/median(S);
  # Q: Is it really supposed to be sqrt(sqrt(...)) here? /HB 2012-01-31
  w1 <- 0.1 + 0.9 * sqrt(sqrt(2*(1-pnorm(abs(CN-2) / median(abs(CN-2))))));

  if (sum(is.nan(w1)) > 0) {
    w1 = rep(1, times=length(w1));
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
    idxsSwap <- references[seq(length=ncol(TR)/2)];
    dataT[1:2,idxsSwap] <- dataT[2:1,idxsSwap, drop=FALSE];

    # Update precalcalculated signals
    TR <- dataT[,references, drop=FALSE];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy numbers must be close to 2 for the reference samples or
  # (if there are not control samples) for most of the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  H <- matrix(2, nrow=nbrOfReferences, ncol=1, byrow=FALSE);
  # Alternatively to the below, on can use method = "MM". /AR 2011-12-04
  fit <- rlm(H ~0 + t(TR), maxit=maxIter, weights=w1);
  # Why not the following? /HB 2011-12-15
  ##  fit <- rlm(x=t(TR), y=H, weights=w1, maxit=maxIter);
  matSum <- fit$coefficients;
  coeffs <- fit$w;
  eps3 <- 1e-8;
  coeffs[coeffs < eps3] <- eps3;
  coeffs <- coeffs * w1;
  dataT <- diag(matSum) %*% dataT;

  # Reextract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  fit <- rlm(naiveGenoDiff ~ 0 + t(TR), maxit=maxIter, weights=coeffs);
  # Why not the following? /HB 2011-12-15
  ##  fit <- rlm(x=t(TR), y=naiveGenoDiff, weights=coeffs, maxit=maxIter);
  matDiff <- fit$coefficients;

  # T matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ T ]   [MatDiff[1] MatDiff[2]]
  U <- matrix(c(0.5, 0.5, 0.5, -0.5), nrow=2, ncol=2, byrow=FALSE);
  V <- matrix(c(c(1, 1), matDiff), nrow=2, ncol=2, byrow=TRUE);
  T <- U %*% V;

  res <- T %*% dataT;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Undo the previous change applied to the data in case there is
  # only one allele
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (onlyOneAllele) {
    res[1:2,idxsSwap] <- res[2:1,idxsSwap, drop=FALSE];
  }

  # Return parameter estimates(?)
  ## attr(res, "modelFit") <- list(fit=fit);

  res;
}, protected=TRUE) # fitCalMaTe()


###########################################################################
# HISTORY:
# 2012-01-31 [HB]
# o CLEANUP: Pruning up the code, e.g. dropping extraneous spaces and
#   parenthesis and adding missing ones. Dropping ambigous comments 
#   using terminologies such as "previous", "latest" and "final".
# o DOCUMENTATION: Argument "references" have to be an index vector
#   (and cannot be a logical vector as previously said).
# 2012-01-31 [MO]
# o BUG FIX: the index "idxs" was recalculated to undo the change when 
#   there is only one allele, and it was done as the previous version,
#   taking into account all the samples, not only the references.     
# 2011-12-15 [HB]
# o CORRECTION: Now doing x[x < eps] <- eps instead of x[x < 0] <- eps.
# o Dropped the validation of argument 'references'.  This is an internal
#   function and validation should already have been done elsewhere.
# 2011-12-07 [MO]
# o At least 3 reference samples.
# 2011-12-04 [AR]
# o BUG FIX: there was a bug for SNPs with a single allele when using a
#   set of references. Set some initial weights based on the median that
#   improves the breakdown point if no reference samples are provided
#   when using a set of references.
# 2011-11-29 [MO]
# o Change matrix "T" by "dataT" and "P" by "T"
# 2010-08-02 [HB]
# o ROBUSTNESS: Now fitCalMaTe() also works (technically) when there is
#   only one reference.
# o Made into an S3 method for matrix:es.
# 2010-06-18 [HB]
# o Created from refineCN.list().
###########################################################################
