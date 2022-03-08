library("calmate")

message("*** thetaAB2TotalAndFracB() ...")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

data <- thetaAB2TotalAndFracB(theta)
str(data)

theta2 <- totalAndFracB2ThetaAB(data)
str(theta2)

stopifnot(all.equal(theta2, theta))

## Exceptions
thetaT <- theta
dim(thetaT) <- c(dim(theta), 1L)
res <- try(dataT <- thetaAB2TotalAndFracB(thetaT))
stopifnot(inherits(res, "try-error"))

thetaT <- theta[,1,,drop=FALSE]
res <- try(dataT <- thetaAB2TotalAndFracB(thetaT))
stopifnot(inherits(res, "try-error"))

dataT <- data
dim(dataT) <- c(dim(data), 1L)
res <- try(thetaT <- totalAndFracB2ThetaAB(dataT))
stopifnot(inherits(res, "try-error"))

dataT <- data[,1,,drop=FALSE]
res <- try(thetaT <- totalAndFracB2ThetaAB(dataT))
stopifnot(inherits(res, "try-error"))

message("*** thetaAB2TotalAndFracB() ... DONE")
