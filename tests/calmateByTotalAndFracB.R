library("calmate")

message("*** calmateByTotalAndFracB() ...")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

# Transform to (total,fracB) signals
data <- thetaAB2TotalAndFracB(theta)

# Calibrate (total,fracB) by CalMaTe
dataC <- calmateByTotalAndFracB(data)

# Calculate copy-number ratios
theta <- data[,"total",]
thetaR <- matrixStats::rowMedians(theta, na.rm=TRUE)
data[,"total",] <- 2*theta/thetaR

# Plot two "random" arrays
Clim <- c(0,4)
Blim <- c(0,1)
subplots(4, ncol=2, byrow=FALSE)
for (ii in c(1,5)) {
  sampleName <- dimnames(data)[[3]][ii]
  sampleLabel <- sprintf("Sample #%d ('%s')", ii, sampleName)
  plot(data[,,ii], xlim=Clim, ylim=Blim)
  title(main=sampleLabel)
  plot(dataC[,,ii], xlim=Clim, ylim=Blim)
  title(main=sprintf("%s\ncalibrated", sampleLabel))
}


# Assert that it also works with a single unit
dummy <- calmateByTotalAndFracB(data[1,,,drop=FALSE])
stopifnot(length(dim(dummy)) == 3)


message("*** calmateByTotalAndFracB() - misc ...")

dataT <- data[1:2,,]
dataU <- truncateFracB(dataT)
str(dataU)
stopifnot(all(dataU[,2,] >= 0))

dataT <- data[1:2,,1]
dataV <- truncateFracB(dataT)
str(dataV)
stopifnot(all(dataV[,2] >= 0))
stopifnot(all.equal(dataV, dataU[,,1]))

message("*** calmateByTotalAndFracB() - misc ... DONE")


message("*** calmateByTotalAndFracB(..., references) ...")

references <- rep(TRUE, length=dim(data)[3])
fit1 <- calmateByTotalAndFracB(data, references=references)
str(fit1)

references <- 1:dim(data)[3]
fit2 <- calmateByTotalAndFracB(data, references=references)
str(fit2)
stopifnot(all.equal(fit2, fit1))

fit3 <- calmateByTotalAndFracB(data, refAvgFcn=matrixStats::rowMedians)
str(fit3)

message("*** calmateByTotalAndFracB(..., references) ... DONE")


message("*** calmateByTotalAndFracB() - exceptions ...")

## Exceptions
dataT <- data[1,,,drop=FALSE]
dim(dataT) <- c(dim(dataT), 1L)
res <- try(fit <- calmateByTotalAndFracB(dataT))
stopifnot(inherits(res, "try-error"))

dataT <- data[1,1,,drop=FALSE]
res <- try(fit <- calmateByTotalAndFracB(dataT))
stopifnot(inherits(res, "try-error"))

references <- rep(TRUE, length=dim(data)[3]-1L)
res <- try(fit <- calmateByTotalAndFracB(data, references=references))
stopifnot(inherits(res, "try-error"))

references <- seq_len(dim(data)[3]+1L)
res <- try(fit <- calmateByTotalAndFracB(data, references=references))
stopifnot(inherits(res, "try-error"))

message("*** calmateByTotalAndFracB() - exceptions ... DONE")

message("*** calmateByTotalAndFracB() ... DONE")
