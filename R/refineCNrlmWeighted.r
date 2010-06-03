###########################################################################/**
# @set "class=list"
# @RdocMethod refineCNrlmWeighted
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
#  \item{input}{A list with two elements. First a 2xI @numeric array, where 2 is the number of alleles,
#               and I is the number of samples. Second an array pointed out the reference samples.}
#  \item{fB1}{Lower heterozygous threshold. Initially set to .33.}
#  \item{fB2}{Higher heterozygous threshold. Initially set to .66.}
# }
#
# \value{
#   Returns an 2xI @numeric array.
# }
#
# @examples "../incl/refineCNrlmWeighted.Rex"
#
# \seealso{
#  To calibrate (total,fracB) data,
#  see @seemethod "weightedCalMaTeByTotalAndFracB".
# }
#*/###########################################################################
setMethodS3("refineCNrlmWeighted", "list", function(input, fB1=0.33, fB2=0.66,..., verbose=FALSE) {
  require("MASS") || stop("Package not loaded: MASS");

  if (!is.list(input)) {
    throw("Argument 'data' is not a list: ", class(input)[1]);
  }

  # Organizing the arguments  
  input <- input[[1]];
  inputData <- input$inputData[[1]];           
  refs <- input$inputRefs[[1]];  
  
  nSamples <- length(inputData)/2;
  dataA <- inputData[1:nSamples];
  dataB <- inputData[(nSamples+1):(2*nSamples)];
  
  if (length(dataA) != length(dataB) || length(refs) != length(dataB)) {
    stop("Wrong input to refineCN function");
  }

  # Axis change
  Tinput <- matrix(data=0, nrow=2, ncol=nSamples);
  Tinput[1,] <- dataA;
  Tinput[2,] <- dataB;  
  
  a <- max(max(Tinput[2,] / (pmax(Tinput[1,],0) + 1e-4)), max(Tinput[1,] / (pmax(Tinput[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2);
  Giro <- solve(Giro);
  Tinput <- Giro %*% Tinput;

  # Checking if all the samples are homozygous
  fracB <- Tinput[2,refs] / (Tinput[1,refs] + Tinput[2,refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  allHomo <- sum(abs(naiveGenoDiff))/2 == length(naiveGenoDiff);

  # Twist Tinput in case only one allele appears
  if (allHomo) {
    n <- round(ncol(Tinput)/2);
    Tinput[c(2,1),1:n] <- Tinput[,1:n];
  }

  # Total copy numbers must be close to 2 for the reference samples or (if there are not control samples) for most of the samples
  Weight <- mean(Tinput) * dim(Tinput)[2];
  outputRlm <- rlm(rbind(t(Tinput[,refs]),c(Weight,-Weight)), rbind(matrix(data=2, nrow=ncol(Tinput[,refs]), ncol=1),0));
  matSum <- outputRlm$coefficients;
  coeffs <- outputRlm$w;
  Tinput <- diag(matSum) %*% Tinput;

  # The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- Tinput[2,refs] / (Tinput[1,refs] + Tinput[2,refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
# matDiff <- rlm(t(Tinput[,refs]), naiveGenoDiff,weights = coeffs[1:(length(coeffs)-1)]);
# matDiff <- rlm(rbind(t(Tinput[,refs]),c(Weight,Weight)), c(naiveGenoDiff,0), weights=coeffs);
  matDiff <- rlm(t(Tinput[,refs]), naiveGenoDiff);

  matDiff <- matDiff$coefficients;

  # P matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ P ]   [MatDiff[1] MatDiff[2]]
  # solve(matrix(c(1,1,1,-1),2,2)) gives matrix(c(.5, .5, .5, -.5),2,2)
  P <- matrix(c(0.5, 0.5, 0.5, -0.5), nrow=2, ncol=2) %*% matrix(c(c(1,1), matDiff), nrow=2, ncol=2, byrow=TRUE);
  Salida <- P %*% Tinput;
  
  # Truncate freqB values to 0 and 1.
  freqB <-  Salida[2,]/colSums(Salida);
  ind <- freqB<0;
  freqB[ind] <- 0;
  ind <- freqB>1;
  freqB[ind] <- 1;
  
  SalidaAux <- Salida;
  SalidaAux[1,] <- colSums(Salida)*(1-freqB);
  SalidaAux[2,] <- colSums(Salida)*(freqB);
  Salida <- SalidaAux;
    
  # Setting Tinput as it was
  if (allHomo) {
    Salida[c(2,1),1:n] <- Salida[,1:n];
  }

  Salida;
}) # refineCN_rlmWeighted()

###########################################################################
# HISTORY:
# 2010-06-3 [MO]
# o Test to check "input" is a list added.
# 2010-06-2 [MO]
# o Comments added.
###########################################################################
