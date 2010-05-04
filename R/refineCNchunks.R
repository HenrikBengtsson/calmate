indexGiro <- function(dataA, dataB){

  Tinput <- matrix(data=0, nrow = 2, ncol = nSamples);
  Tinput[1,] <- dataA;
  Tinput[2,] <- dataB;  

  a <-  max(max(Tinput[2,] / (pmax(Tinput[1,],0) + 1e-4)), max(Tinput[1,] / (pmax(Tinput[2,],0) + 1e-4)));  
}

refineCN <- function (DataA, DataB, Refs, fB1=0.33,fB2=0.66)
{
  require("MASS") || stop("Package not loaded: MASS");

  save(DataA, DataB,Refs, file="CalMaTeTry.Rdata")
  if(is.null(dim(DataA)) || is.null(dim(DataB)) || length(Refs) != ncol(DataB)){
    stop("Wrong input to refineCN function")
  }
  #Axis change
  nSamples <- length(Refs)
    
#  Tinput <- matrix(data=0, nrow = 2, ncol = nSamples);
#  Tinput[1,] <- DataA;
#  Tinput[2,] <- DataB;  
  
  a <- apply(DataA, 1, FUN=indexGiro(DataB));      
  a <- max(max(Tinput[2,] / (pmax(Tinput[1,],0) + 1e-4)), max(Tinput[1,] / (pmax(Tinput[2,],0) + 1e-4)));
  Giro <- matrix(c(1, 1/a, 1/a, 1), nrow=2, ncol=2);
  
  Giro <- solve(Giro);
  Tinput <- Giro%*%Tinput;

  #Check if all the samples are homozygous
  fracB <- Tinput[2,Refs] / (Tinput[1,Refs] + Tinput[2,Refs]);
  naiveGenoDiff <- 2*(fracB < fB1) - 2*(fracB > fB2);
  AllHomo <- sum(abs(naiveGenoDiff))/2 == length(naiveGenoDiff);
  #Dar la vuelta a Tinput
  if (AllHomo==TRUE)
  {
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
  matDiff <- rlm(t(Tinput[,Refs]), naiveGenoDiff);
  matDiff <- matDiff$coefficients;

  # P matrix is:
  #  [1  1] [   ] = [MatSum[1]   MatSum[2]] (We have already applied it) MatSum is 1,1
  #  [1 -1] [ P ]   [MatDiff[1] MatDiff[2]]
  # solve(matrix(c(1,1,1,-1),2,2)) gives matrix(c(.5, .5, .5, -.5),2,2)
  P <- matrix(c(.5, .5, .5, -.5),2,2) %*% matrix(c(c(1,1), matDiff),2,2,byrow=TRUE)
  Salida <- P%*%Tinput;

  # Volver a dar la vuelta
  if (AllHomo==TRUE)
  {
    Salida[c(2,1),1:n] <- Salida[,1:n];
  }
  return(Salida);
}