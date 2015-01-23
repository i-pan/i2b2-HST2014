## Set starting time
ptm <- proc.time()

## Set seed
set.seed(11)

## Load libraries
library(MASS)

## Generate parameters
beta0 <- runif(806, min = -3, max = -2)
beta1 <- rnorm(806, mean = 0.5, sd = 0.1)
beta2 <- runif(806, min = 0, max = 0.02)

## Generate number of facts ~ Poisson
Fcases <- rpois(5000, 50)
Fcontrols <- rpois(100000, 30)

### Develop logistic regression model: logit(ODDS) = beta0 + beta1*[1 if case, 0 if control] + beta2*[# of facts] + susceptX
# Define expit function
expit <- function(x) exp(x)/(1+exp(x)) # solve for probability of getting disease

## Generate multivariate normal variable (susceptibility to disease)
print('Generating disease susceptibility for cases...')
# Generate covariance matrix
d <- 806 
covmat <- matrix(rnorm(d^2, mean = 0, sd = 0.02), d, d)
covmat <- covmat + t(covmat)
covmat <- cor(covmat)
# Generate random variable
susceptX <- mvrnorm(n = 5000, mu = rep(0, length = 806), Sigma = covmat)

# Generate cases
print('Generating cases...')
diseaseCASE <- matrix(0, nrow = length(Fcases), ncol = 807)
indexrow <- 1
for(eachcase in Fcases) {
	print(paste('Generating case', as.character(indexrow), '...'))
	logitCASE <- beta0 + beta1 + beta2*eachcase + 0.25*susceptX[indexrow, ]
	Pcases <- expit(logitCASE)
	indexcol <- 1
	for(eachP in Pcases) {
		diseaseCASE[indexrow, indexcol] <- rbinom(1, 1, eachP)
		indexcol <- indexcol + 1
	}
	diseaseCASE[indexrow, 807] <- eachcase
	indexrow <- indexrow + 1
}
print('Finished generating cases')
# Generate controls - I actually ended up using the exact same distribution of facts as I did for the cases
# so I wouldn't have to do a matching algorithm to compute relative risk later on.
# Originally, I planned to generate a larger distribution of facts for controls so that the matching
# algorithm could be applied, but this is easily adjustable

## Generate multivariate normal variable (susceptibility to disease)
print('Generating disease susceptibility for controls...')
# Generate covariance matrix
d <- 806 
covmat <- matrix(rnorm(d^2, mean = 0, sd = 0.02), d, d)
covmat <- covmat + t(covmat)
covmat <- cor(covmat)
# Generate random variable
susceptX <- mvrnorm(n = 5000, mu = rep(0, length = 806), Sigma = covmat)

print('Generating controls...')
diseaseCONTROL <- matrix(0, nrow = length(Fcontrols), ncol = 807)
indexrow <- 1
for(eachcont in Fcases) {
	print(paste('Generating control', as.character(indexrow), '...'))
	logitCONT <- beta0 + beta2*eachcont + 0.25*susceptX[indexrow, ]
	Pconts <- expit(logitCONT)
	indexcol <- 1
	for(eachP in Pconts) {
		diseaseCONTROL[indexrow, indexcol] <- rbinom(1, 1, eachP)
		indexcol <- indexcol + 1
	}
	diseaseCONTROL[indexrow, 807] <- eachcont
	indexrow <- indexrow + 1
}
print('Finished generating controls')

## Compute relative risks
print('Computing relative risks...')
relrisk <- numeric(806)
for(eachdisease in 1:806) {
	print(paste('Computing relative risk for disease', as.character(eachdisease), '...'))
	eachrisk <- sum(diseaseCASE[, eachdisease])/sum(diseaseCONTROL[, eachdisease]) # not sure if this approach is valid
	relrisk[eachdisease] <- eachrisk
}
print('Finished computing relative risks')
riskdensity <- density(relrisk)
hist(relrisk, prob=T, col='grey', ylim=c(0,max(riskdensity$y*1.2)))
lines(density(relrisk), col='blue', lwd=2)
print('Printing summary of relative risks...')
print(summary(relrisk))

## Calculate runtime
runtime <- proc.time() - ptm
print('Calculating runtime...')
print(runtime)