###########################################################################/**
# @RdocClass CalMaTeNormalization
#
# @title "The CalMaTeNormalization class"
#
# \description{
#  @classhierarchy
#
#  This class represents the CalMaTe normalization method [1], which 
#  corrects for SNP effects in allele-specific copy-number estimates
#  (ASCNs).
# }
# 
# @synopsis 
#
# \arguments{
#   \item{data}{A named @list with data set named \code{"total"} and
#     \code{"fracB"} where the former should be of class
#     @see "aroma.core::AromaUnitTotalCnBinarySet" and the latter of
#     class @see "aroma.core::AromaUnitFracBCnBinarySet".  The
#     two data sets must be for the same chip type, have the same
#     number of samples and the same sample names.}
#   \item{tags}{Tags added to the output data sets.}
#   \item{...}{Arguments passed to calmateByTotalAndFracB or calmateByThetaAB.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# \details{
#   ...
# }
#
# \examples{\dontrun{
#   @include "../incl/CalMaTeNormalization.Rex"
# }}
#
# \references{
#   [1] ...
# }
#
# \seealso{
#   Low-level versions of the CalMaTe normalization method is available
#   via the @see "calmateByThetaAB.array" and 
#   @see "calmateByTotalAndFracB.array" methods.
# }
#
#*/###########################################################################
setConstructorS3("CalMaTeNormalization", function(data=NULL, tags="*", ...) {
  # Validate arguments
  if (!is.null(data)) {
    if (!is.list(data)) {
      throw("Argument 'data' is not a list: ", class(data)[1]);
    }
    reqNames <- c("total", "fracB");
    ok <- is.element(reqNames, names(data));
    if (!all(ok)) {
      throw(sprintf("Argument 'data' does not have all required elements (%s): %s", paste(reqNames, collapse=", "), paste(reqNames[!ok], collapse=", ")));
    }
    data <- data[reqNames];

    # Assert correct classes
    className <- "AromaUnitTotalCnBinarySet";
    ds <- data$total;
    if (!inherits(ds, className)) {
      throw(sprintf("The 'total' data set is not of class %s: %s", className, class(ds)[1]));
    }

    className <- "AromaUnitFracBCnBinarySet";
    ds <- data$fracB;
    if (!inherits(ds, className)) {
      throw(sprintf("The 'fracB' data set is not of class %s: %s", className, class(ds)[1]));
    }

    # Assert that the chip types are compatile
    if (getChipType(data$total) != getChipType(data$fracB)) {
      throw("The 'total' and 'fracB' data sets have different chip types: ", 
            getChipType(data$total), " != ", getChipType(data$fracB));
    }

    # Assert that the data sets have the same number data files
    nbrOfFiles <- nbrOfFiles(data$total)
    if (nbrOfFiles != nbrOfFiles(data$fracB)) {
      throw("The number of samples in 'total' and 'fracB' differ: ", 
            nbrOfFiles, " != ", nbrOfFiles(data$fracB));
    }

    # Assert that the data sets have the same samples
    if (!identical(getNames(data$total), getNames(data$fracB))) {
      throw("The samples in 'total' and 'fracB' have different names.");
    }

  }


  # Arguments '...':
  args <- list(...);
  if (length(args) > 0) {
    argsStr <- paste(names(args), collapse=", ");
    throw("Unknown arguments: ", argsStr);
  }

  this <- extend(Object(...), "CalMaTeNormalization",
    .data = data
  );

  setTags(this, tags);

  this; 
})


setMethodS3("as.character", "CalMaTeNormalization", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  dsList <- getDataSets(this);
  s <- c(s, sprintf("Data sets (%d):", length(dsList)));
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    s <- c(s, sprintf("<%s>:", capitalize(names(dsList)[kk])));
    s <- c(s, as.character(ds));
  } 
 
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



setMethodS3("getAsteriskTags", "CalMaTeNormalization", function(this, collapse=NULL, ...) {
  tags <- "CMTN";

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }
  
  tags;
}, private=TRUE) 


setMethodS3("getName", "CalMaTeNormalization", function(this, ...) {
  dsList <- getDataSets(this);
  ds <- dsList$total;
  getName(ds);
}) 



setMethodS3("getTags", "CalMaTeNormalization", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  dsList <- getDataSets(this);
  ds <- dsList$total;
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("setTags", "CalMaTeNormalization", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})


setMethodS3("getFullName", "CalMaTeNormalization", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getDataSets", "CalMaTeNormalization", function(this, ...) {
  this$.data;
})
 

setMethodS3("getRootPath", "CalMaTeNormalization", function(this, ...) {
  "totalAndFracBData";
})


setMethodS3("getPath", "CalMaTeNormalization", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  dsList <- getDataSets(this);
  ds <- dsList$total;
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})


setMethodS3("nbrOfFiles", "CalMaTeNormalization", function(this, ...) {
  dsList <- getDataSets(this);
  ds <- dsList$total;
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSets", "CalMaTeNormalization", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  res <- this$.outputDataSets;
  if (is.null(res)) {
    res <- allocateOutputDataSets(this, ..., verbose=less(verbose, 10));
    this$.outputDataSets <- res;
  }
  res;
}) 


setMethodS3("allocateOutputDataSets", "CalMaTeNormalization", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  verbose && enter(verbose, "Retrieve/allocation output data sets");

  dsList <- getDataSets(this);
  path <- getPath(this);

  res <- list();
  for (kk in seq(along=dsList)) {
    ds <- dsList[[kk]];
    verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                 kk, getName(ds), length(dsList)));

    for (ii in seq(ds)) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                        ii, getName(df), nbrOfFiles(ds)));

      filename <- getFilename(df);
      pathname <- Arguments$getWritablePathname(filename, path=path, 
                                                       mustNotExist=FALSE);
      # Skip?
      if (isFile(pathname)) {
        verbose && cat(verbose, "Already exists. Skipping.");
        verbose && exit(verbose);
        next;
      }

      # Create temporary file
      pathnameT <- sprintf("%s.tmp", pathname);
      pathnameT <- Arguments$getWritablePathname(pathnameT, mustNotExist=TRUE);

      # Copy source file
      copyFile(getPathname(df), pathnameT);

      # Make it empty by filling it will missing values
      # AD HOC: We should really allocate from scratch here. /HB 2010-06-21
      dfT <- newInstance(df, pathnameT);
      dfT[,1] <- NA;

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Renaming temporary file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Renaming temporary output file");
      file.rename(pathnameT, pathname);
      if (!isFile(pathname)) {
        throw("Failed to rename temporary file ('", pathnameT, "') to final file ('", pathname, "')");
      }
      verbose && exit(verbose);

      verbose && cat(verbose, "Copied: ", pathname);
      verbose && exit(verbose);
    } # for (ii ...)

    dsOut <- byPath(ds, path=path, ...);

    # AD HOC: The above byPath() grabs all *.asb files. /HB 2010-06-20
    keep <- is.element(sapply(dsOut, getFilename), sapply(ds, getFilename));
    dsOut <- extract(dsOut, keep);

    res[[kk]] <- dsOut;
    rm(ds, dsOut);
    
    verbose && exit(verbose);
  } # for (kk ...)

  names(res) <- names(dsList);

  this$.outputDataSets <- res;

  verbose && exit(verbose);

  res;
}, protected=TRUE)




setMethodS3("findUnitsTodo", "CalMaTeNormalization", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


  verbose && enter(verbose, "Finding units to do");

  dsList <- getOutputDataSets(this);

  # The last data set is always updated last
  ds <- dsList[[length(dsList)]];
  verbose && print(verbose, ds);

  # The last file (in lexicographic ordering) is always updated last
  fullnames <- getFullNames(ds);
  verbose && str(verbose, fullnames);

  o <- order(fullnames, decreasing=TRUE);
  idx <- o[1];
  df <- getFile(ds, idx);
  verbose && print(verbose, df);

  # Read all values
  values <- df[,1,drop=TRUE];

  verbose && cat(verbose, "Number of units: ", length(values));


  # Identify all missing values
  nok <- is.na(values);
  units <- which(nok);

  verbose && printf(verbose, "Number of units to do: %d (%.2f%%)\n", 
                       length(units), 100*length(units)/length(values));

  verbose && cat(verbose, "Units to do (with missing values):");
  verbose && str(verbose, units);

  verbose && exit(verbose);

  units;
})


setMethodS3("process", "CalMaTeNormalization", function(this, units="remaining", ..., force=FALSE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsList <- getDataSets(this);
  dsTCN <- dsList$total;
  dsBAF <- dsList$fracB;
  
  # Argument 'units':
  if (is.null(units)) {
    units <- "remaining";
  }
  df <- getFile(dsTCN, 1);
  nbrOfUnits <- nbrOfUnits(df);
  if (identical(units, "remaining")) {
    units <- seq(length=nbrOfUnits);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'ram':
  ram <- getRam(aromaSettings, ram);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "CalMaTe normalization of ASCNs");
  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of arrays: ", nbrOfFiles);

  verbose && cat(verbose, "Units:");
  verbose && str(verbose, units);

  # Skip already processed units or not?
  if (!force) {
    verbose && enter(verbose, "Skipping already processed units");
    unitsTodo <- findUnitsTodo(this, verbose=less(verbose, 25));
    verbose && cat(verbose, "Units remaining:");
    verbose && str(verbose, unitsTodo);

    units <- intersect(units, unitsTodo);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);

    rm(unitsTodo);
    verbose && exit(verbose);
  }

  nbrOfUnits <- length(units);
  verbose && printf(verbose, "Number of units to do: %d (%.2f%%)\n", 
                      length(units), 100*length(units)/nbrOfUnits(df));


  chipType <- getChipType(dsTCN, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  rm(dsList);


  sampleNames <- getNames(dsTCN);
  dimnames <- list(NULL, sampleNames, c("total", "fracB"));

  outPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate output data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- getOutputDataSets(this, verbose=less(verbose, 5));
  if (nbrOfUnits == 0) {
    verbose && cat(verbose, "No more units to process. Skipping.");
    verbose && exit(verbose);
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Process in chunks
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calculating number of units to fit per chunk");
  verbose && cat(verbose, "RAM scale factor: ", ram);

  bytesPerChunk <- 100e6;       # 100Mb
  verbose && cat(verbose, "Bytes per chunk: ", bytesPerChunk);

  bytesPerUnitAndArray <- 2*8;  # Just a rough number; good enough?
  verbose && cat(verbose, "Bytes per unit and array: ", bytesPerUnitAndArray);

  bytesPerUnit <- nbrOfFiles * bytesPerUnitAndArray;
  verbose && cat(verbose, "Bytes per unit: ", bytesPerUnit);

  unitsPerChunk <- ram * bytesPerChunk / bytesPerUnit;
  unitsPerChunk <- as.integer(max(unitsPerChunk, 1));
  unitsPerChunk <- min(unitsPerChunk, nbrOfUnits);
  verbose && cat(verbose, "Number of units per chunk: ", unitsPerChunk);

  nbrOfChunks <- ceiling(nbrOfUnits / unitsPerChunk);
  verbose && cat(verbose, "Number of chunks: ", nbrOfChunks);
  verbose && exit(verbose);

  idxs <- 1:nbrOfUnits;
  head <- 1:unitsPerChunk;

  count <- 1L;
  while (length(idxs) > 0) {
    tTotal <- processTime();

    verbose && enter(verbose, sprintf("Chunk #%d of %d", count, nbrOfChunks));
    if (length(idxs) < unitsPerChunk) {
      head <- 1:length(idxs);
    }
    uu <- idxs[head];

    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units[uu]);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Reading (total,fracB) data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Reading (total,fracB) data");
    total <- extractMatrix(dsTCN, units=units[uu], verbose=less(verbose,25));
    verbose && str(verbose, total);
    fracB <- extractMatrix(dsBAF, units=units[uu], verbose=less(verbose,25));
    verbose && str(verbose, fracB);
    verbose && exit(verbose);

    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Normalizing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Combining into an (total,fracB) array");
    dim <- c(nrow(total), ncol(total), 2);
    data <- c(total, fracB);
    rm(total, fracB);
    data <- array(data, dim=dim, dimnames=dimnames);
    data <- aperm(data, perm=c(1,3,2));
    verbose && str(verbose, data);
    verbose && exit(verbose);

    verbose && enter(verbose, "Normalizing");
    dataN <- calmateByTotalAndFracB(data, ..., verbose=verbose);
    fit <- attr(dataN, "modelFit");
    verbose && str(verbose, fit);
    verbose && str(verbose, dataN);
    verbose && exit(verbose);

    rm(data);  # Not needed anymore
    gc <- gc();
    verbose && print(verbose, gc);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing model fit
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing model fit");
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Storing normalized data
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Storing normalized data");
    for (kk in seq(along=res)) {
      ds <- res[[kk]];
      verbose && enter(verbose, sprintf("Data set #%d ('%s') of %d", 
                                           kk, getName(ds), length(res)));

      # Store in lexicograph ordering
      fullnames <- getFullNames(ds);
      idxs <- order(fullnames, decreasing=FALSE);
      
      for (ii in idxs) {
        df <- getFile(ds, ii);
        verbose && enter(verbose, sprintf("Data file #%d ('%s') of %d", 
                                        ii, getName(df), nbrOfFiles(ds)));

        signals <- dataN[,kk,ii];
        verbose && cat(verbose, "Signals:");
        verbose && str(verbose, signals);
        df[units[uu],1] <- signals;
        rm(signals);

        verbose && exit(verbose);
      } # for (ii ...)

      verbose && exit(verbose);
    } # for (kk ...)
    verbose && exit(verbose);

    # Next chunk
    idxs <- idxs[-head];
    count <- count + 1;

    # Garbage collection
    tGc <- processTime();
    gc <- gc();
    verbose && print(verbose, gc);

    verbose && exit(verbose);
  } # while(length(idxs) > 0)

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2010-07-22
# o BUG FIX: Now process() for CalMaTeNormalization returns immediately 
#   if there are no units left.
# o BUG FIX: process(..., verbose=TRUE) would give "Error in sprintf("Chunk
#   #%d of %d", count, nbrOfChunks) :  invalid format '%d'; use format %f, 
#   %e, %g or %a for numeric objects".
# 2010-06-29
# o Added support for process(..., force=TRUE) in CalMaTeNormalization.
# 2010-06-23
# o Added support for argument 'references' and 'truncate'.
#   TODO: The asterisk tags should probably reflect these settings.
# o Code cleanup.
# 2010-06-21
# o Added minor Rdoc comments.
# o ROBUSTNESS: Added more assertions to the CalMaTe constructor.
# o ROBUSTNESS: Now process() stores results in lexicograph ordering to
#   assure that the lexicographicly last file is updated last, which is 
#   the file that findUnitsTodo() is querying.
# o Added getOutputDataSets().
# o Added findUnitsTodo().
# o Now allocateOutputDataSets() "clears" the files.
# 2010-06-20
# o First test shows that it seems to work.
# o Created.
############################################################################
