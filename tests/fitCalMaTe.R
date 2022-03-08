library("calmate")

message("*** fitCalMaTe() ...")

# Load example (thetaA,thetaB) signals
path <- system.file("exData", package="calmate")
theta <- loadObject("thetaAB,100x2x40.Rbin", path=path)

theta1 <- theta[1,,]
references <- 1:10

for (flavor in c("v1", "v2")) {
  message(sprintf(" - fitCalMaTe(..., flavor='%s') ...", flavor))
  fit <- fitCalMaTe(theta1, references=references, flavor=flavor)
  str(fit)
  stopifnot(all.equal(dim(fit), dim(theta1)))
  message(sprintf(" - fitCalMaTe(..., flavor='%s') ... DONE", flavor))
}

message("*** fitCalMaTe() ... DONE")
