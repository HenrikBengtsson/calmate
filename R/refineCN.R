###########################################################################/**
# @set "class=list"
# @RdocMethod refineCN
#
# @title "CalMaTe calibration function"
#
# \description{
#  @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{input}{A @list with two elements. 
#      The first element is a 2xI @numeric array, where 2 is the number 
#      of alleles, and I is the number of samples. 
#      The second element is an array pointed out the reference samples.}
#  \item{fB1}{Lower heterozygous threshold in [0,1]. Initially set to 1/3.}
#  \item{fB2}{Higher heterozygous threshold in [0,1]. Initially set to 2/3.}
#  \item{maxIter}{Maximum number of iterations for @see "MASS::rlm" function.}
#  \item{...}{Not used.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an 2xI @numeric array.
# }
#
# @examples "../incl/refineCN.Rex"
#
# \seealso{
#  To calibrate (total,fracB) data,
#  see @seemethod "calmateByTotalAndFracB".
# }
#*/###########################################################################
setMethodS3("refineCN", "list", function(input, fB1=1/3, fB2=2/3, maxIter=50, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'input':
  if (!is.list(input)) {
    throw("Argument 'data' is not a list: ", class(input)[1]);
  }

  # Argument 'fB1' & 'fB2':
  fB1 <- as.double(fB1);
  if (fB1 < 0 || fB1 > 1) {
    throw("Argument 'fB1' is out of range [0,1]: ", fB1);
  }

  fB2 <- as.double(fB2);
  if (fB2 < 0 || fB2 > 1) {
    throw("Argument 'fB2' is out of range [0,1]: ", fB2);
  }

  # Argument 'maxIter':
  maxIter <- as.integer(maxIter);
  if (length(maxIter) != 1) {
    throw("Argument 'maxIter' has to be a single value: ", length(maxIter));
  }
  if (!is.finite(maxIter)) {
    throw("Argument 'maxIter' is non-finite: ", maxIter);
  }
  if (maxIter <= 0) {
    throw("Argument 'maxIter' is non-positive: ", maxIter);
  }


  # Organizing the arguments  
  input <- input[[1]];
  inputData <- input$inputData[[1]];    
  # Adding a small value so there are "non" 0 values
  eps <- 1e-6;
  inputData[inputData < eps] = eps;
         
  refs <- input$inputRefs[[1]];  
  
  nSamples <- length(inputData) / 2;
  dataA <- inputData[1:nSamples];
  dataB <- inputData[(nSamples+1):(2*nSamples)];

  # Sanity check  
  if (length(dataA) != length(dataB) || length(refs) != length(dataB)) {
    stop("Wrong input to refineCN() function");
  }

  # Axis change
  data <- c(dataA, dataB);
  Tinput <- matrix(data, nrow=2, ncol=nSamples, byrow=TRUE);
  
  a <- max(max(Tinput[2,] / (pmax(Tinput[1,],0) + 1e-4)), max(Tinput[1,] / (pmax(Tinput[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2, byrow=FALSE);
  Giro <- solve(Giro);
  Tinput <- Giro %*% Tinput;

  # Checking if all the samples are homozygous
  fracB <- Tinput[2,refs] / (Tinput[1,refs] + Tinput[2,refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);         
  oneAllele <- (abs(sum(naiveGenoDiff)/2) == length(naiveGenoDiff));

  # Twist half of the samples in case there is only one allele
  if (oneAllele) {
    Tinput[1:2,1:(ncol(Tinput)/2)] <- Tinput[2:1,1:(ncol(Tinput)/2)];
  }

  # Total copy numbers must be close to 2 for the reference samples or (if there are not control samples) for most of the samples
  fit <- rlm(t(Tinput[,refs]), matrix(data=2, nrow=ncol(Tinput[,refs]), ncol=1, byrow=FALSE), maxit=maxIter);
  matSum <- fit$coefficients;
  coeffs <- fit$w;
  Tinput <- diag(matSum) %*% Tinput;

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- Tinput[2,refs] / (Tinput[1,refs] + Tinput[2,refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  fit <- rlm(t(Tinput[,refs]), naiveGenoDiff, maxit=maxIter, weights=coeffs);
  matDiff <- fit$coefficients;

  # P matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ P ]   [MatDiff[1] MatDiff[2]]
  P <- matrix(c(0.5, 0.5, 0.5, -0.5), nrow=2, ncol=2, byrow=FALSE) %*% matrix(c(c(1,1), matDiff), nrow=2, ncol=2, byrow=TRUE);
  
  res <- P %*% Tinput;

  # Calculate total CNs
  C <- colSums(res);

  # Truncate fracB values to 0 and 1.
  fracB <-  res[2,] / C;
  eps <- 1e-5;
  fracB[(fracB < eps)] <- eps;  # Why?!? /HB
  fracB[(fracB > 1)] <- 1;      # Why not here?!? /HB

  res[1,] <- C*(1-fracB);
  res[2,] <- C*(fracB);

  # Correct the previous change applied to the data in case there is 
  # only one allele    
  if (oneAllele) {
    res[1:2,1:(ncol(res)/2)] <- res[2:1,1:(ncol(res)/2)];
  }

  res;
}) # refineCN()

###########################################################################
# HISTORY:
# 2010-06-18 [HB]
# o CLEAN UP: Renamed variables to non-Spanish etc.
# o BUG FIX: Now "truncating" by x[x < eps] <- eps (was x[x == 0] <- eps).
# 2010-06-04 [MO]
# o Created.
###########################################################################
