calMaTeWeighted <- function(dataA, dataB, refs=0) {
  #Checking the reference information
  createdRefs <- FALSE;
  nSamples <- 0;
  if (is.null(dim(dataA))) {
    nSamples <- length(dataA);
    nSNPs <- 1;
  } else {
    nSamples <- ncol(dataA);
    nSNPs <- nrow(dataA);
  }

  #there are not given references
  if (is.null(dim(refs)) && length(refs) == 1) {
    refs <- rep(TRUE, times=nSamples);
    createdRefs <- TRUE;
  }

  #in the case it is only one SNP
  if (is.null(dim(dataA)) && is.null(dim(refs))) {
    if (is.logical(refs) && length(refs) == nSamples) {
      createdRefs <- TRUE;
    }

    if (!is.logical(refs)) {
      aux <- rep(FALSE, times=nSamples);
      aux[refs] <- TRUE;
      refs <- aux;
      createdRefs <- TRUE;
    }

    if (!createdRefs) {
      stop("Wrong reference information");
    } else {
      refs <- t(as.matrix(refs));
      dataA <- t(as.matrix(dataA));
      dataB <- t(as.matrix(dataB));
    }
  } else {
    #it is not a matrix as the one that contains the data
    if (is.null(dim(refs)) || dim(refs) != dim(dataA) || !is.logical(refs)) { 
      #it is an index vector
      if (length(refs) != nSamples && !is.logical(refs)) {
        aux <- matrix(data=FALSE, ncol=ncol(dataB), nrow=nrow(dataA));
        aux[,refs] <- TRUE;
        refs <- aux;
      }

      #it is a logical reference vector with more than one sample as reference
      if (is.logical(refs) && length(refs) == ncol(dataA) && sum(refs) > 1) {
        createdRefs <- TRUE;
      }

      #non of these cases
      if (!createdRefs) {
        stop("Wrong reference information");
      } else {
        #generate the reference matrix
        pr <- rep(refs, times=nrow(dataA));
        dim(pr) <- c(length(refs), nrow(dataA));
        refs <- t(pr);
      }  
    }
  }

  inputData <- apply(cbind(dataA,dataB), MARGIN=1, FUN=list);
  inputRefs <- apply(refs, MARGIN=1, FUN=list);
  input <- cbind(inputData, inputRefs);
  input <- apply(input, MARGIN=1, FUN=list);
  refineData <- lapply(X=input, FUN=refineCN_rlmWeighted);
# refineData <- lapplyInChunks(c(1:nSNPs), function(rr) {
#   refineCN_rlmWeighted(input[rr]);
# }, chunkSize=500e3);

  refineData;
} # calMaTeWeighted()
