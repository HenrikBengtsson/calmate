fitCalMaTe <- function(T, refs, fB1=1/3, fB2=2/3, maxIter=50, truncate=TRUE, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  nbrOfSNPs <- nrow(T);
  nSamples <- ncol(T);
  nbrOfRefs <- length(refs);

  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  T[T < 1e-6] <- 1e-6;
  
  a <- max(max(T[2,] / (pmax(T[1,],0) + 1e-4)), max(T[1,] / (pmax(T[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  T <- Giro %*% T;

  # Extract the signals for the reference set
  TR <- T[,refs];

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
    idxs <- seq(length=ncol(T)/2);
    T[1:2,idxs] <- T[2:1,idxs];

    # Update precalcalculated signals
    TR <- T[,refs];
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Total copy numbers must be close to 2 for the reference samples or
  # (if there are not control samples) for most of the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  H <- matrix(2, nrow=nbrOfRefs, ncol=1, byrow=FALSE);
  fit <- rlm(t(TR), H, maxit=maxIter);
  matSum <- fit$coefficients;
  coeffs <- fit$w;
  T <- diag(matSum) %*% T;

  # Reextract the signals for the reference set
  TR <- T[,refs];

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  fit <- rlm(t(TR), naiveGenoDiff, maxit=maxIter, weights=coeffs);
  matDiff <- fit$coefficients;

  # P matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ P ]   [MatDiff[1] MatDiff[2]]
  P <- matrix(c(0.5, 0.5, 0.5, -0.5), nrow=2, ncol=2, byrow=FALSE) %*% matrix(c(c(1,1), matDiff), nrow=2, ncol=2, byrow=TRUE);
  
  res <- P %*% T;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Truncate ASCNs to avoid non-positives?
  # Truncation is done such that TCN is preserved regardlessly.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (truncate) {
    delta <- matrix(0, nrow=nrow(res), ncol=ncol(res));

    # (CA*,CB*) = (CA,CB) + (-delta,+delta) for all SNPs where CA < 0.
    idxs <- which(res[1,] < 0);
    d <- res[1,idxs];
    delta[1,idxs] <- -d;
    delta[2,idxs] <- +d;

    # (CA*,CB*) = (CA,CB) + (+delta,-delta) for all SNPs where CB < 0.
    idxs <- which(res[2,] < 0);
    d <- res[2,idxs];
    delta[1,idxs] <- +d;
    delta[2,idxs] <- -d;

    # Note, here all(colSums(delta) == 0) is TRUE.
    res <- res + delta;
  } # if (!allowNegative)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Undo the previous change applied to the data in case there is 
  # only one allele    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (onlyOneAllele) {
    idxs <- seq(length=ncol(res)/2);
    res[1:2,idxs] <- res[2:1,idxs];
  }

  res;
} # fitCalMaTe()


###########################################################################
# HISTORY:
# 2010-06-19 [HB]
# o Added argument 'truncate' for optional truncating of (CA,CB).
# 2010-06-18 [HB]
# o Created from refineCN.list().
###########################################################################
