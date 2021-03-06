library("calmate")

message("*** calmateByThetaAB() ...")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

# Calculate (CA,CB)
thetaR <- matrixStats::rowMedians(theta[,"A",] + theta[,"B",], na.rm=TRUE)
C <- 2*theta/thetaR

# For each available CalMaTe fitting algorithm...
flavors <- sort(eval(formals(calmateByThetaAB.array)$flavor))

CList <- list(raw=C)
for (flavor in flavors) {
  # Calibrate (CA,CB) by CalMaTe
  CList[[flavor]] <- calmateByThetaAB(theta, flavor=flavor)

  # Assert that it also works with a single unit
  dummy <- calmateByThetaAB(theta[1,,,drop=FALSE], flavor=flavor)
  stopifnot(length(dim(dummy)) == 3)
}

# Create plot
Clim <- c(-0.2,4)

if (interactive()) {
  devNew(type="x11", aspectRatio=1.9)
} else {
  devNew(type="png", "test-calmateByTheta.png", aspectRatio=1.9)
}

subplots(2*length(CList)+2, ncol=2, byrow=FALSE)
par(mar=c(3,3,1,1)+0.1, mgp=c(1.8,0.7,0))
# Plot two "random" arrays
for (ii in c(1,5)) {
  for (kk in 1:length(CList)) {
    key <- names(CList)[kk]
    C <- CList[[key]]
    plot(C[,,ii], col=kk, xlim=Clim, ylim=Clim)

    if (kk == 1L) {
      sampleName <- dimnames(C)[[3]][ii]
      label <- sprintf("Sample #%d ('%s')", ii, sampleName)
    } else {
      label <- sprintf("calibrated (%s)", key)
    }
    stext(side=3, pos=0.5, label)
  } # for (kk ...)

  for (kk in 1:length(CList)) {
    key <- names(CList)[kk]
    C <- CList[[key]]

    if (kk == 1L) {
      plot(C[,,ii], col=kk, xlim=Clim, ylim=Clim)
    } else {
      points(C[,,ii], col=kk)
    }
    stext(side=3, pos=0.5, "All together")
  } # for (kk ...)
} # for (ii ...)

devDone()


("*** calmateByThetaAB() - NAs ...")

thetaT <- theta
thetaT[,2,4] <- NA_real_
fit <- calmateByThetaAB(thetaT)
str(fit)

("*** calmateByThetaAB() - NAs ... DONE")


("*** calmateByThetaAB() - misc ...")

thetaT <- theta[1:2,,]
fit <- calmateByThetaAB(thetaT, truncate=TRUE)
str(fit)

fit <- calmateByThetaAB(thetaT, refAvgFcn=matrixStats::rowMedians)
str(fit)

thetaU <- truncateThetaAB(thetaT)
str(thetaU)
stopifnot(all(thetaU >= 0))

thetaV <- theta[1:2,,1]
thetaV <- truncateThetaAB(thetaV)
str(thetaV)
stopifnot(all(thetaV >= 0))
stopifnot(all.equal(thetaV, thetaU[,,1]))

references <- rep(TRUE, length=dim(theta)[3])
fit <- calmateByThetaAB(theta, references=references)
str(fit)

("*** calmateByThetaAB() - misc ... DONE")


message("*** calmateByThetaAB() - exceptions ...")

thetaT <- theta[1,,,drop=FALSE]
dim(thetaT) <- c(dim(thetaT), 1L)
res <- try(fit <- calmateByThetaAB(thetaT))
stopifnot(inherits(res, "try-error"))

thetaT <- theta[1,1,,drop=FALSE]
res <- try(fit <- calmateByThetaAB(thetaT))
stopifnot(inherits(res, "try-error"))

references <- rep(TRUE, length=dim(theta)[3]-1L)
res <- try(fit <- calmateByThetaAB(theta, references=references))
stopifnot(inherits(res, "try-error"))

references <- seq_len(dim(theta)[3]+1L)
res <- try(fit <- calmateByThetaAB(theta, references=references))
stopifnot(inherits(res, "try-error"))

message("*** calmateByThetaAB() - exceptions ... DONE")

message("*** calmateByThetaAB() ... DONE")
