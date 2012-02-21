###########################################################################/**
# @set "class=matrix"
# @RdocMethod fitCalMaTeMedians
# @alias fitCalMaTeMedians
#
# @title "Calibrates SNP loci according to the CalMaTe method, not using "rlm" function" .
#
# \description{
#  @get "title".
#  \emph{Note: This is an internal function of the package, which is kept
#   only kept to provide easy access to the internal fit functions.
#   It it actually not elsewhere in the package, and should nor by others.}
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
#  \item{...}{Additional arguments passed to the internal fit functions.}
# }
#
# \value{
#   Returns a 2xI @numeric @matrix of calibrated ASCNs.
# }
#*/###########################################################################

setMethodS3("fitCalMaTeMedians", "matrix", function(T, references, fB1=1/3, fB2=2/3,...) {
  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  data <- T;
  nbrOfSNPs <- nrow(T);
  nbrOfReferences <- length(references);

  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  T[T < eps] <- eps;
  
  a <- max(max(T[2,references] / (pmax(T[1,references],0) + 1e-4)), max(T[1,references] / (pmax(T[2,references],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  T <- Giro %*% T;

  # Extract the signals for the reference set
  TR <- T[,references, drop=FALSE];
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Genotyping
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Group the three possibilities
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  wgths <- c(sum(naiveGenoDiff==2),
             sum(naiveGenoDiff==0),
             sum(naiveGenoDiff==-2))
  
  TR <- cbind(rowMedians(TR[,naiveGenoDiff==2, drop=FALSE]),
              rowMedians(TR[,naiveGenoDiff==0, drop=FALSE]),
              rowMedians(TR[,naiveGenoDiff==-2,drop=FALSE]));
  # Remove possible NaNs
  TR[is.nan(TR)] <- 0;
  
  # Handle cases a single genotype is present in the samples

  # A very unlikely case: all the samples are heterozygous
  if ((wgths[1] == 0)&(wgths[3] == 0)) {
    TR[1,1] <- 2*TR[1,2];
    TR[2,1] <- .1*TR[2,2];

    wgths[1] <- wgths[2]
  }
  # Only BB
  if ((wgths[1] == 0)&(wgths[2] == 0)) {
    TR[,1] <- TR[2:1,3];
    wgths[1] <- wgths[3]
  }
  if ((wgths[3] == 0)&(wgths[2] == 0)) {
    TR[,3] <- TR[2:1,1];
    wgths[3] <- wgths[1]
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Callibration
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  Genotypes <- cbind(c(2,0), c(1,1), c(0,2));
  P <- qr.solve(t(TR) * wgths^2, t(Genotypes) * wgths^2);
  res <- t(P) %*% T;
  
  T <- res;

  TR <- T[,references, drop=FALSE];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Have the genotypes changed?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fracB <- TR[2,] / (TR[1,] + TR[2,]);
  naiveGenoDiff_1 <- 2*(fracB < fB1) - 2*(fracB > fB2);
  
  if (!identical(naiveGenoDiff_1, naiveGenoDiff)) {
        # Repeat the process once
        naiveGenoDiff <- naiveGenoDiff_1;
        wgths <- c(sum(naiveGenoDiff==2),
                   sum(naiveGenoDiff==0),
                   sum(naiveGenoDiff==-2))

        TR <- cbind(rowMedians(TR[,naiveGenoDiff==2, drop=FALSE]),
                    rowMedians(TR[,naiveGenoDiff==0, drop=FALSE]),
                    rowMedians(TR[,naiveGenoDiff==-2,drop=FALSE]));
        # Remove possible NaNs
        TR[is.nan(TR)] <- 0;

        # Handle cases in which not all the genotypes are present in the samples

        # A very unlikely case: all the samples are heterozygous
        if ((wgths[1] == 0)&(wgths[3] == 0)) {
          TR[1,1] <- 2*TR[1,2];
          TR[2,1] <- .1*TR[2,2];

          wgths[1] <- wgths[2]
        }
        # Only BB
        if ((wgths[1] == 0)&(wgths[2] == 0)) {
          TR[,1] <- TR[2:1,3];
          wgths[1] <- wgths[3]
        }
        if ((wgths[3] == 0)&(wgths[2] == 0)) {
          TR[,3] <- TR[2:1,1];
          wgths[3] <- wgths[1]
        }

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Callibration
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Genotypes <- cbind(c(2,0), c(1,1), c(0,2));
        P <- qr.solve(t(TR) * wgths^2, t(Genotypes) * wgths^2);
        res <- t(P) %*% T;
  }
  res;
}, protected=TRUE) # fitCalMaTeMedians()

###########################################################################
# HISTORY:
# 2011-04-12 [AR]
# o Created from fitCalMaTe().
#   · the previous version used rlm. This one does not.
#   · Robustness is obtained by using medians among different genotypes
#       and weighting according to the number or samples in each genotype
###########################################################################
