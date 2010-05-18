refineCN_rlmWeighted <- function (input, fB1=0.33,fB2=0.66)
{

  require("MASS") || stop("Package not loaded: MASS");

#  save(input, file="input.Rdata")
  input <- input[[1]];
  inputData <- input$inputData[[1]];
  Refs <- input$inputRefs[[1]];  
  
  nSamples <- length(inputData)/2;
  DataA <- inputData[1:nSamples];
  DataB <- inputData[(nSamples+1):(2*nSamples)];
  
  if(length(DataA)!= length(DataB) || length(Refs) != length(DataB)){
    stop("Wrong input to refineCN function")
  }
  #Axis change
    
  Tinput <- matrix(data=0, nrow = 2, ncol = nSamples);
  Tinput[1,] <- DataA;
  Tinput[2,] <- DataB;  
  
  a <- max(max(Tinput[2,] / (pmax(Tinput[1,],0) + 1e-4)), max(Tinput[1,] / (pmax(Tinput[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2);
  Giro <- solve(Giro);
  Tinput <- Giro%*%Tinput;

  #Check if all the samples are homozygous
  fracB <- Tinput[2,Refs] / (Tinput[1,Refs] + Tinput[2,Refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  AllHomo <- sum(abs(naiveGenoDiff))/2 == length(naiveGenoDiff);

  #Twist Tinput in case only one allele appears
  if(AllHomo){
    n=round(ncol(Tinput)/2);
    Tinput[c(2,1),1:n] <- Tinput[,1:n];
  }

  #Total copy numbers must be close to 2 for the reference samples or (if there are not control samples) for most of the samples
  Weight <- mean(Tinput) * dim(Tinput)[2];
  Outputrlm <- rlm(rbind(t(Tinput[,Refs]),c(Weight,-Weight)),rbind(matrix(data=2, nrow=ncol(Tinput[,Refs]), ncol=1),0));
  matSum <- Outputrlm$coefficients;
  coeffs <- Outputrlm$w;
  Tinput <- diag(matSum) %*% Tinput;

  #The difference of the copy numbers must be 2, 0 or -2 depending genotyping
  fracB <- Tinput[2,Refs] / (Tinput[1,Refs] + Tinput[2,Refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
# matDiff <- rlm(t(Tinput[,Refs]), naiveGenoDiff,weights = coeffs[1:(length(coeffs)-1)]);
  matDiff <- rlm(rbind(t(Tinput[,Refs]),c(Weight,Weight)), c(naiveGenoDiff,0),weights = coeffs);

  matDiff <- matDiff$coefficients;

  # P matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ P ]   [MatDiff[1] MatDiff[2]]
  # solve(matrix(c(1,1,1,-1),2,2)) gives matrix(c(.5, .5, .5, -.5),2,2)
  P <- matrix(c(.5, .5, .5, -.5), nrow=2, ncol=2) %*% matrix(c(c(1,1), matDiff), nrow=2, ncol=2, byrow=TRUE)
  Salida <- P%*%Tinput;

  # Setting Tinput as it was
  if (AllHomo){Salida[c(2,1),1:n] <- Salida[,1:n];}
  return(Salida);

}