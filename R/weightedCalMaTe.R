###########################################################################/**
# @set "class=matrix"
# @RdocMethod weightedCalMaTe
#
# @title "Internal Calmate function"
#
# \description{
#  @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{dataA}{Matrix with the number of copies of A allele of J SNPs in the I samples.}
#  \item{dataB}{Matrix with the number of copies of B allele of J SNPs in the I samples.}
#  \item{refs}{Vector indicating the reference samples.}
# }
#
# \value{
#   Returns a list with J elements (number of SNPs) an each containing a 2xI @numeric array.
# }
#
#
# \seealso{
#  To calibrate (total,fracB) data,
#  see @seemethod "weightedCalMaTeByTotalAndFracB".
#  see @seemethod "refineCNrlmWeighted".
# }
#*/###########################################################################
setMethodS3("weightedCalMaTe", "matrix", function(dataA, dataB, refs=0,..., verbose=FALSE) {

 if (!is.matrix(dataA)) {
    throw("Argument 'data' is not a matrix: ", class(dataA)[1]);
  }

  #Checking the reference information
  createdrefs <- FALSE;
  nSamples <- 0;
  if(is.null(dim(dataA))){
    nSamples <- length(dataA);
    nSNPs <- 1;
  }else{
    nSamples <- ncol(dataA);
    nSNPs <- nrow(dataA);
  }
  #there are not given references
  if(is.null(dim(refs)) && length(refs)==1) 
  {
    refs <- rep(TRUE, nSamples)
    createdrefs <- TRUE;    
  }
  #in the case it is only one SNP
  if(is.null(dim(dataA)) && is.null(dim(refs))){
    if(is.logical(refs) && length(refs) == nSamples){
      createdrefs <- TRUE;
    }
    if(!is.logical(refs)){
      aux <- rep(FALSE, nSamples);
      aux[refs] <- TRUE;
      refs <- aux;
      createdrefs <- TRUE;
    }
    if(createdrefs == FALSE){
      stop("Wrong reference information")
    }else{
      refs <- t(as.matrix(refs));
      dataA <- t(as.matrix(dataA));
      dataB <- t(as.matrix(dataB));
    }    
  }else{
    #it is not a matrix as the one that contains the data
    if(is.null(dim(refs)) || dim(refs)!=dim(dataA) || !is.logical(refs)){ 
      #it is an index vector
      if(length(refs) != nSamples && !is.logical(refs)){
        aux <- matrix(data=FALSE, ncol=ncol(dataB), nrow=nrow(dataA))
        aux[,refs] <- TRUE;
        refs <- aux;
      }
      #it is a logical reference vector with more than one sample as reference
      if(is.logical(refs) && length(refs) == ncol(dataA) && sum(refs) > 1){
        createdrefs <- TRUE;            
      }
      #non of these cases
      if(!createdrefs){
        stop("Wrong reference information")
      }else{
        #generate the reference matrix
        pr <- rep(refs, nrow(dataA));
        dim(pr) <- c(length(refs), nrow(dataA));
        refs <- t(pr);
      }  
    }
  }

#  inputData <- cbind(dataA, dataB)
  
  inputData <- apply(cbind(dataA,dataB),1,list)
  inputRefs <- apply(refs,1,list)
  input <- cbind(inputData, inputRefs);   
  input <- apply(input,1,list)
  refineData <- lapply(X=input, FUN = refineCNrlmWeighted);
#  refineData <- lapplyInChunks(c(1:nSNPs), function(rr) {
#  refineCN_rlmWeighted(input[rr]);}, chunkSize=500e3)

  return(refineData)
}) # weightedCalMaTe() 

###########################################################################
# HISTORY:
# 2010-06-3 [MO]
# o Comments added.
###########################################################################
