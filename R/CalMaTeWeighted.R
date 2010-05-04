CalMaTeWeighted <- function (DataA, DataB, Refs = 0)
{

  #Checking the reference information
  createdRefs <- FALSE;
  if(is.null(dim(DataA))){
    nSamples <- length(DataA);
    nSNPs <- 1;
  }else{
    nSamples <- ncol(DataA);
    nSNPs <- nrow(DataA);
  }
  #there are not given references
  if(is.null(dim(Refs)) && length(Refs)==1) 
  {
    Refs <- rep(TRUE, nSamples)
    createdRefs <- TRUE;    
  }
  #in the case it is only one SNP
  if(is.null(dim(DataA)) && is.null(dim(Refs))){
    if(is.logical(Refs) && length(Refs) == nSamples){
      createdRefs <- TRUE;
    }
    if(!is.logical(Refs)){
      aux <- rep(FALSE, nSamples);
      aux[Refs] <- TRUE;
      Refs <- aux;
      createdRefs <- TRUE;
    }
    if(createdRefs == FALSE){
      stop("Wrong reference information")
    }else{
      Refs <- t(as.matrix(Refs));
      DataA <- t(as.matrix(DataA));
      DataB <- t(as.matrix(DataB));
    }    
  }else{
    #it is not a matrix as the one that contains the data
    if(is.null(dim(Refs)) || dim(Refs)!=dim(DataA) || !is.logical(Refs)){ 
      #it is an index vector
      if(length(Refs) != nSamples && !is.logical(Refs)){
        aux <- matrix(data=FALSE, ncol=ncol(DataB), nrow=nrow(DataA))
        aux[,Refs] <- TRUE;
        Refs <- aux;
      }
      #it is a logical reference vector with more than one sample as reference
      if(is.logical(Refs) && length(Refs) == ncol(DataA) && sum(Refs) > 1){
        createdRefs <- TRUE;            
      }
      #non of these cases
      if(createdRefs == FALSE){
        stop("Wrong reference information")
      }else{
        #generate the reference matrix
        pr <- rep(Refs, nrow(DataA));
        dim(pr) <- c(length(Refs), nrow(DataA));
        Refs <- t(pr);
      }  
    }
  }


#  inputData <- cbind(DataA, DataB)
  
  inputData <- apply(cbind(DataA,DataB),1,list)
  inputRefs <- apply(Refs,1,list)
  input <- cbind(inputData, inputRefs);   
  input <- apply(input,1,list)
  refineData <- lapply(X=input, FUN = refineCN_rlmWeighted);
#  refineData <- lapplyInChunks(c(1:nSNPs), function(rr) {
#  refineCN_rlmWeighted(input[rr]);}, chunkSize=500e3)

  return(refineData)
}