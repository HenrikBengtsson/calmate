datasets <- list(geo=list(), tcga=list());

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# GSE12702
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
datasets$geo$sampleName <- "GSM318736";
  
datasets$geo$segments <- list(
  "normal (1,1)" = c(8,  0.0, 16.0),
  "loss (0,1)"   = c(8, 16.5, 45.0), 
  "gain (1,2)"   = c(8, 45.0,150.0)
);

datasets$geo$getPairs <- function(dsList, ...) {
  # For the GSE12702 data set, it happens to be that samples, when sorted
  # lexicographically, are ordered in pairs of tumors and normals.
  sampleNames <- getNames(dsList$total);
  sampleNames <- matrix(sampleNames, nrow=20, ncol=2, byrow=TRUE);
  colnames(sampleNames) <- c("tumor", "normal");

  pairs <- array(sampleNames, dim=dim(sampleNames));
  colnames(pairs) <- c("tumor", "normal");
  rownames(pairs) <- sampleNames[,"tumor"];

  ## When could add patient ID annotation too.
  ## rownames(pairs) <- c(24, 25, 27, 31, 45, 52, 58, 60, 75, 110, 115, 122,
  ##                      128, 137, 138, 140, 154, 167, 80, 96);

  pairs;
} # getPairs()
  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# TCGA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
datasets$tcga$sampleName <- "TCGA-23-1027";

datasets$tcga$segments <- list(
  "normal (1,1)" = c(2, 105.0, 124.0),
  "gain (1,2)"   = c(2, 125.0, 140.0), 
  "CN-LOH (0,2)" = c(2, 150.0, 170.0)
);

datasets$tcga$getPairs <- function(dsList, ...) {
    fullnames <- getFullNames(dsList$total);
    sampleNames <- gsub(",(total|fracB).*", "", fullnames);
    sampleNames <- gsub("-[0-9]{2}[A-Z]-[0-9]{4}-[0-9]{2}.*", "", sampleNames);

    names <- gsub("-[0-9]{2}[A-Z]", "", sampleNames);
    uNames <- unique(names);
    pairs <- matrix(NA, nrow=length(uNames), ncol=2);
    rownames(pairs) <- uNames;
    colnames(pairs) <- c("tumor", "normal");
    for (name in uNames) {
      patternT <- sprintf("%s-(01)[A-Z]$", name);
      patternN <- sprintf("%s-(10|11)[A-Z]$", name);
      idxT <- grep(patternT, sampleNames);
      idxN <- grep(patternN, sampleNames);
      pair <- c(idxT[1], idxN[1]);
      pair <- sampleNames[pair];
      pairs[name,] <- pair;
    }
    pairs;
} # getPairs()

