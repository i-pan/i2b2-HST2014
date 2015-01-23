## Load libraries
library(MASS)

## Generate parameters
beta0 <- runif(806, min = -1, max = 1)
beta1 <- rnorm(806, mean = 0, sd = 0.05)
beta2 <- runif(806, min = 0, max = 0.02)

## Generate number of facts ~ Poisson
Fcases <- rpois(5000, 50)
Fcontrols <- rpois(100000, 30)

## Generate multivariate normal variable (susceptibility to disease)
# Generate covariance matrix
d <- 806 
covmat <- matrix(rnorm(d^2, mean = 0, sd = 0.02), d, d)
covmat <- covmat + t(covmat)
covmat <- cor(covmat)
# Generate random variable
susceptX <- mvrnorm(n = 1, mu = rep(0, length = 806), Sigma = covmat)

## Develop logistic regression model: logit(ODDS) = beta0 + beta1*[1 if case, 0 if control] + beta2*[# of facts] + susceptX
# Define expit function
expit <- function(x) exp(x)/(1+exp(x)) # solve for probability of getting disease
# Generate cases
diseaseCASE <- matrix(0, nrow = length(Fcases), ncol = 807)
indexrow <- 1
for(eachcase in Fcases) {
	logitCASE <- beta0 + beta1 + beta2*eachcase + 0.25*susceptX
	Pcases <- expit(logitCASE)
	Pcases <- Pcases/10
	indexcol <- 1
	for(eachP in Pcases) {
		diseaseCASE[indexrow, indexcol] <- rbinom(1, 1, eachP)
		indexcol <- indexcol + 1
	}
	diseaseCASE[indexrow, 807] <- eachcase
	indexrow <- indexrow + 1
}
# Generate controls
diseaseCONTROL <- matrix(0, nrow = length(Fcontrols), ncol = 807)
indexrow <- 1
for(eachcont in Fcontrols) {
	logitCONT <- beta0 + beta2*eachcont + 0.25*susceptX
	Pconts <- expit(logitCONT)
	Pconts <- Pconts/10
	indexcol <- 1
	for(eachP in Pconts) {
		diseaseCONTROL[indexrow, indexcol] <- rbinom(1, 1, eachP)
		indexcol <- indexcol + 1
	}
	diseaseCONTROL[indexrow, 807] <- eachcont
	indexrow <- indexrow + 1
}