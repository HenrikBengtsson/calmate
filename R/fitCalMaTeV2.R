fitCalMaTeV2 <- function(dataT, references, fB1=1/3, fB2=2/3, maxIter=50, ...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  dataInit <- dataT;
  nbrOfSNPs <- nrow(dataT);
  nbrOfReferences <- length(references);

  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  dataT[dataT < eps] <- eps;

  a <- max(max(dataT[2,references] / (pmax(dataT[1,references],0) + 1e-4)), max(dataT[1,references] / (pmax(dataT[2,references],0) + 1e-4)));
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
  # Twist half of the samples in case there is only one allele
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
  # Alternatively to the below, on can use method = "MM". /AR 2011-12-04
  suppressWarnings(fit <- rlm(H ~0 + t(TR), maxit=maxIter, weights=w1));
  if (fit$converged == FALSE){ 
    return(fitCalMaTeMedians(dataInit, references, fB1=1/3, fB2=2/3));
  }

  matSum <- fit$coefficients;
  coeffs <- fit$w;
  eps2 <- 1e-8;
  coeffs[coeffs < eps2] <- eps2;
  coeffs <- coeffs * w1;
  dataT <- diag(matSum) %*% dataT;

  #Reextract the signals for the reference set
  TR <- dataT[,references, drop=FALSE];

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  if(abs(sum(naiveGenoDiff)/2) == length(naiveGenoDiff)){
    return(fitCalMaTeMedians(dataInit, references, fB1=1/3, fB2=2/3)); 
  }
  suppressWarnings(fit <- rlm(naiveGenoDiff ~ 0 + t(TR), maxit=maxIter, weights=coeffs));
  if (fit$converged == FALSE){
    return(fitCalMaTeMedians(dataInit, references, fB1=1/3, fB2=2/3));
  }

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
    res[1:2,idxs] <- res[2:1,idxs, drop=FALSE];
  }

  res;
} # fitCalMaTeV2()


###########################################################################
# HISTORY:
# 2012-02-21 [MO]
# o Created from fitCalMaTeV1
# o Use of weigths on "rlm" function
# o Warnings suppresed
# o Use of fitCalMaTeMedians when rlm does not converge
###########################################################################
