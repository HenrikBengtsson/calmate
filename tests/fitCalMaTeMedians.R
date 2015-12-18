library("calmate")

message("*** fitCalMaTeMedians() ...")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

references <- 1:10

## Typical
theta1 <- theta[1,,]
fit <- fitCalMaTeMedians(theta1, references=references)
str(fit)
stopifnot(all.equal(dim(fit), dim(theta1)))

message("- fitCalMaTeMedians() - all AA, BB, or heterozygous ...")

## All samples are AA
theta1 <- theta[1,,]
theta1[2,] <- 0
fit <- fitCalMaTeMedians(theta1, references=references)
str(fit)
stopifnot(all.equal(dim(fit), dim(theta1)))

## All samples are BB
theta1 <- theta[1,,]
theta1[1,] <- 0
fit <- fitCalMaTeMedians(theta1, references=references)
str(fit)
stopifnot(all.equal(dim(fit), dim(theta1)))

## All samples are heterozygous
theta1 <- theta[1,,]
theta1[2,] <- theta1[1,]
fit <- fitCalMaTeMedians(theta1, references=references)
str(fit)
stopifnot(all.equal(dim(fit), dim(theta1)))

message("- fitCalMaTeMedians() - all AA, BB, or heterozygous ... DONE")

message("*** fitCalMaTeMedians() ... DONE")
