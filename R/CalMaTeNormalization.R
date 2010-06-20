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

    # Assert that the chip types are compatile
    if (getChipType(data$total) != getChipType(data$fracB)) {
      throw("The 'total' and 'fracB' data sets have different chip types: ", 
            getChipType(data$total), " != ", getChipType(data$fracB));
    }

    # Assert that the data sets have the same number data files
    if (nbrOfFiles(data$total) != nbrOfFiles(data$fracB)) {
      throw("The number of samples in 'total' and 'fracB' differ: ", 
            nbrOfFiles(data$total), " != ", nbrOfFiles(data$fracB));
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
  ds <- getInputDataSet(this);
  getName(ds);
}) 



setMethodS3("getTags", "CalMaTeNormalization", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
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
 

setMethodS3("getInputDataSet", "CalMaTeNormalization", function(this, ...) {
  data <- getDataSets(this, ...);
  data$total;
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
  ds <- getInputDataSet(this);
  chipType <- getChipType(ds, fullname=FALSE);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
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
  ds <- getInputDataSet(this);
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSet", "CalMaTeNormalization", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- byPath(ds, path=path, ...);
  res;
}) 


setMethodS3("allocateOutputDataSets", "CalMaTeNormalization", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 


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

      # Copy source file
      copyFile(getPathname(df), pathname);

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

  res;
}, protected=TRUE)


setMethodS3("process", "CalMaTeNormalization", function(this, units=NULL, ..., force=FALSE, ram=NULL, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsList <- getDataSets(this);
  dsTCN <- dsList$total;
  dsBAF <- dsList$fracB;

  # Argument 'units':
  df <- getFile(dsTCN, 1);
  nbrOfUnits <- nbrOfUnits(df);
  if (is.null(units)) {
    units <- seq(length=nbrOfUnits);
  } else {
    units <- Arguments$getIndices(units, max=nbrOfUnits);
  }
  nbrOfUnits <- length(units);

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

  chipType <- getChipType(dsTCN, fullname=FALSE);
  verbose && cat(verbose, "Chip type: ", chipType);
  rm(dsList);

  sampleNames <- getNames(dsTCN);
  dimnames <- list(NULL, sampleNames, c("total", "fracB"));

  outPath <- getPath(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate output data sets
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- allocateOutputDataSets(this, verbose=less(verbose, 1)); 


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

  count <- 1;
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
    total <- extractMatrix(dsTCN, units=units[uu], verbose=less(verbose,5));
    verbose && str(verbose, total);
    fracB <- extractMatrix(dsBAF, units=units[uu], verbose=less(verbose,5));
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
    dataN <- calmateByTotalAndFracB(data, verbose=less(verbose,5));
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

      for (ii in seq(ds)) {
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

#  res <- getOutputDataSet(this, verbose=less(verbose, 1)); 

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2010-06-20
# o First test shows that it seems to work.
# o Created.
############################################################################
