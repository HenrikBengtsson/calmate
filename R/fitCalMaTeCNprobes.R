fitCalMaTeCNprobes <- function(T, references,...) {

  # This is an internal function. Because of this, we will assume that
  # all arguments are valid and correct.  No validation will be done.
  
  refCN <- rowMedians(T[,references]);
  dataCNC <- 2*T/refCN;

  res;
} # fitCalMaTeCNprobes()


###########################################################################
# HISTORY:
# 2010-06-22 [MO]
# o Created.
###########################################################################
