fitCalMaTeCNprobes <- function(T, references,...) {

  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  
  #refCN <- rowMedians(T[,references]);
  #res <- 2*T/refCN;
  
  nbrOfReferences <- length(references);
  H <- matrix(2, nrow=nbrOfReferences, ncol=1,byrow=FALSE);
  matSum <- matrix(0, nrow=1, ncol=nrow(T));

  ind <- which(rowSums(T)!=0);
  fit <- matrix(0, nrow=1, ncol=length(ind));

  for(ii in ind){
    fit[ii] <- rlm(T[ii,references], H, maxit=50);
  }
  
  if(ind != 0)
    matSum[ind] <- fit$coefficients;
  res <- matSum %*% T;

  res;
} # fitCalMaTeCNprobes()


###########################################################################
# HISTORY:
# 2010-06-22 [MO]
# o Created.
###########################################################################
