datasets <- list();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# GSE12702
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
datasets$geo <- env({
"Mapping250K_Nsp" <- list(
  dataSet  = "GSE12702",
  tags     = "ACC,-XY,BPN,-XY,RMA,FLN,-XY",
  chipType = "Mapping250K_Nsp"
);

sampleName <- "GSM318736";
  
segments <- list(
  "normal (1,1)" = c(8,  0.0, 16.0),
  "loss (0,1)"   = c(8, 16.5, 45.0), 
  "gain (1,2)"   = c(8, 45.0,150.0)
);

getPairs <- function(dsList, ...) {
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
}); # geo  

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# TCGA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
datasets$tcga <- env({
"Human1M-Duo" <- list(
  dataSet  = "hudsonalpha.org_OV.Human1MDuo.1.1.0",
  tags     = "XY",
  chipType = "Human1M-Duo"
);

"GenomeWideSNP_6" <- list(
  dataSet  = "broad.mit.edu_OV.Genome_Wide_SNP_6.12.6.0",
  tags     = "ASCRMAv2",
  chipType = "GenomeWideSNP_6"
);

sampleName <- "TCGA-23-1027";

segments <- list(
  "normal (1,1)" = c(2, 105.0, 124.0),
  "gain (1,2)"   = c(2, 125.0, 140.0), 
  "CN-LOH (0,2)" = c(2, 150.0, 170.0)
);

getPairs <- function(dsList, ...) {
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
}) # tcga


datasets <- lapply(datasets, FUN=function(env) as.list(env));

# Identify available chip types
chipTypes <- unlist(datasets);
chipTypes <- chipTypes[grep("chipType", names(chipTypes))];
chipTypes <- unlist(chipTypes, use.names=FALSE);

calTags <- c("<none>", "CMTN", "CMTN,refs=N")[-2];
anTags <- c("dens", "hets", "2Mb");

hetCol <- c("#999999", "red", "blue", "orange")[4];


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Choose chip type to study
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (interactive() && require("R.menu")) {
  chipType <- textMenu(chipTypes, title="Select chip type:", value=TRUE);
} else {
  chipType <- chipTypes[1];
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (chipType == "Mapping250K_Nsp") {
  ds <- datasets$geo;
} else if (chipType == "Human1M-Duo") {
  ds <- datasets$tcga;
} else if (chipType == "GenomeWideSNP_6") {
  ds <- datasets$tcga;
} 
attachLocally(ds);
attachLocally(ds[[chipType]]);
