###########################################################################/**
# @set "class=matrix"
# @RdocMethod calmate
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
#  \item{dataA}{A @numeric JxI copy-number @matrix, 
#        where J is the number of SNPs and I is the number samples.}
#  \item{dataB}{A @numeric JxI copy-number @matrix.}
#  \item{refs}{A @logical @vector indicating the reference samples.}
#  \item{...}{Additional arguments passed to 
#        @seemethod "refineCN".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @list with J elements (SNPs) each containing 
#   a 2xI @numeric @matrix.
# }
#
# @examples "../incl/calmate.Rex"
#
# \seealso{
#  To calibrate (total,fracB) data,
#  see @seemethod "calmateByTotalAndFracB".
#  see @seemethod "refineCN".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("calmate", "matrix", function(dataA, dataB, refs=0, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'dataA':
  if (!is.matrix(dataA)) {
    throw("Argument 'dataA' is not a matrix: ", class(dataA)[1]);
  }

  # Argument 'dataB':
  if (!is.matrix(dataB)) {
    throw("Argument 'dataB' is not a matrix: ", class(dataB)[1]);
  }
  if (!identical(dim(dataB), dim(dataA))) {
    dimAStr <- paste(dim(dataA), collapse="x");
    dimBStr <- paste(dim(dataB), collapse="x");
    throw(sprintf("Arguments 'dataA' and 'dataB' have different dimensions: (%s) != (%s)", dimAStr, dimBStr));
  }

  # Argument 'refs':




  nSamples <- ncol(dataA);
  nSNPs <- nrow(dataA);

  #Checking the reference information
  createdRefs <- FALSE;

  #if there are not given references     
  if (is.null(dim(refs)) && is.numeric(refs) && 
                            length(refs) == 1 && refs == 0) {
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
    if (is.null(dim(refs)) || !identical(dim(refs), dim(dataA)) || !is.logical(refs)) { 
      #it is an index vector
      if (length(refs) != nSamples && !is.logical(refs)) {
        if (max(refs) > ncol(dataB)) {
          stop("Wrong reference information");  
        }
        aux <- matrix(data=FALSE, ncol=ncol(dataB), nrow=1, byrow=FALSE);
        aux[1,refs] <- TRUE;
        refs <- aux;
        createdRefs <- TRUE;
      }

      #it is a logical reference vector with more than one sample as reference
      if (is.logical(refs) && length(refs) == ncol(dataA) && sum(refs) >= 1) {
        createdRefs <- TRUE;            
      }

      #non of these cases
      if (!createdRefs) {
        stop("Wrong reference information");
      } else {
        #generate the reference matrix
        aux <- rep(refs, times=nrow(dataA));
        dim(aux) <- c(length(refs), nrow(dataA));
        refs <- t(aux);
      }  
    }
  }

  inputData <- apply(cbind(dataA,dataB), MARGIN=1, FUN=list);
  inputRefs <- apply(refs, MARGIN=1, FUN=list);
  input <- cbind(inputData, inputRefs);   
  input <- apply(input, MARGIN=1, FUN=list);
  refineData <- base::lapply(input, FUN=refineCN, ...);

  return(refineData);
}) # calmate() 

###########################################################################
# HISTORY:
# 2010-06-18 [HB]
# o BUG FIX: Argument '...' was not passed to refineCN().
# o BUG FIX: Invalid test (dim(x) != dim(y)) was used.
# o CLEAN UP: Dropped argument 'maxIter', which is now passed in '...'.
# o CLEAN UP: Tidying up code.  More is needed.
# 2010-06-04 [MO]
# o Created.
###########################################################################
