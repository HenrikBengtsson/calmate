###########################################################################/**
# @set "class=array"
# @RdocMethod calmateByTotalAndFracB
# @alias calmateByTotalAndFracB
# 
# @title "Normalize allele-specific copy numbers (total,fracB)"
#
# \description{
#  @get "title", where total is the total (non-polymorphic) signal and
#  fracB is the allele B fraction.
#  It is only loci with a non-missing (@NA) fracB value that are
#  considered to be SNPs and normalized by CalMaTe.  The other loci
#  are left untouched.
# }
#
# @synopsis
#
# \arguments{
#  \item{data}{An Jx2xI @numeric @array, where J is the number of loci,
#                      2 is total and fracB, and I is the number of samples.}
#  \item{references}{A @logical or @numeric @vector specifying which
#     samples should be used as the reference set.  
#     By default, all samples are considered.}
#  \item{...}{Additional arguments passed to 
#         @seemethod "calmateByThetaAB".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns an Jx2xI @numeric @array.
# }
#
# @examples "../incl/calmateByTotalAndFracB.Rex"
#
# \seealso{
#  To calibrate (thetaA,thetaB) or (CA,CB) signals, 
#  see @seemethod "calmateByThetaAB".
# }
#*/###########################################################################
setMethodS3("calmateByTotalAndFracB", "array", function(data, references=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'data':
  if (!is.array(data)) {
    throw("Argument 'data' is not an array: ", class(data)[1]);
  }
  dim <- dim(data);
  dimnames <- dimnames(data);
  if (length(dim) != 3) {
    throw("Argument 'data' is not a 3-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (dim[2] != 2) {
    throw("Argument 'data' is not a Jx2xI-dimensional array: ", 
                                                paste(dim, collapse="x"));
  }
  if (!is.null(dimnames[[2]])) {
    if (dimnames[[2]][2]=="freqB"){
      dimnames(data)[[2]][2] = "fracB";  
      dimnames[[2]][2]="fracB";
    }    
    if (!identical(dimnames[[2]], c("total", "fracB"))) {
      throw("If given, the names of the allele (2nd) dimension of the Jx2xI-dimensional array (argument 'data') have to be 'total' & 'fracB': ", paste(dimnames[[2]], collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);


  verbose && enter(verbose, "calmateByTotalAndFracB()");
  verbose && cat(verbose, "(total,fracB) signals:");
  verbose && str(verbose, data);

  verbose && enter(verbose, "Identifying SNPs (non-missing fracB)");
  nok <- is.na(data[,"fracB",]);
  nok <- rowAlls(nok);
  snps <- which(!nok);
  verbose && printf(verbose, "Number of SNPs: %d (%.2f%%)\n",
                            length(snps), 100*length(snps)/dim(data)[1]);
  verbose && exit(verbose);

  verbose && enter(verbose, "Transforming SNPs to (thetaA, thetaB)");
  theta <- data[snps,,,drop=FALSE];
  theta <- totalAndFracB2ThetaAB(theta, verbose=less(verbose, 5));
  verbose && str(verbose, theta);
  verbose && exit(verbose);

  thetaC <- calmateByThetaAB(theta, references=references, ..., verbose=verbose);  
  rm(theta); # Not needed anymore

  verbose && enter(verbose, "Backtransforming SNPs to (total, fracB)");
  dataC <- data;
  dataC[snps,,] <- thetaAB2TotalAndFracB(thetaC, verbose=less(verbose, 5));
  verbose && str(verbose, dataC);
 
  rm(snps); # Not needed anymore
  verbose && exit(verbose);

  verbose && enter(verbose, "Calibrating non-polymorphic probes");
  dataC[nok,"total",] <- fitCalMaTeCNprobes(data[nok,"total",], references=references);
  verbose && str(verbose, dataC[nok,,]);
  verbose && exit(verbose);
  
  
  verbose && cat(verbose, "Calibrated (total,fracB) signals:");
  verbose && str(verbose, dataC);

  verbose && exit(verbose);

  dataC;
}) # calmateByTotalAndFracB()


###########################################################################
# HISTORY:
# 2010-06-22 [MO]
# o Now calmateByTotalAndFracB() calibrates also non-polymorphic loci.
# 2010-06-18 [HB]
# o Now calmateByTotalAndFracB() handles also non-polymorphic loci.
# 2010-06-04 [MO]
# o Created.
###########################################################################
