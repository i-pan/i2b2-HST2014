## Runtime calculation
ptm <- proc.time()

## Load packages and set seed
library(MASS)
library(plyr)
library(survival)
library(spcov)
set.seed(10)

## Generate indicator variables to indicate which predictors will apply to 
## which diseases

p <- 0.05  # insert probability here
numpre <- 100  # insert number of predictors
numdis <- 806  # insert number of diseases
indicators <- matrix(0, nrow = numdis, ncol = numpre)
for(eachrow in 1:numdis) {
	indicators[eachrow, ] = rbinom(numpre, 1, p)
}

## Generate beta parameters for each predictor for each disease 

numcomb <- 50  # number of comorbidities
# for cases
betas.cases <- matrix(rnorm(numdis*numpre, mean = 0, sd = log(1.5)/2), nrow = 
	numdis, ncol = numpre)
for(eachbeta in 1:numdis) {
	betas.cases[eachbeta, ] <- betas.cases[eachbeta, ]*indicators[eachbeta, ]
}
# for controls
betas.controls <- betas.cases
for(eachdis in 1:numcomb) {
	for(eachbeta in 1:numpre) {
		if(betas.controls[eachdis, eachbeta] != 0) {
			betas.controls[eachdis, eachbeta] <- 
				betas.controls[eachdis, eachbeta] - log(2)/7
		}
	}
}

## Generate predictors for each patient

beta0 <- runif(806, -2, 0)
numcases <- 500  # insert number of cases
numconts <- 5000  # insert number of controls
covmat.cases <- GenerateCliquesCovariance(numpre/2, 2, log(1.5)/2)$Sigma
pre.cases <- mvrnorm(numcases, mu = rep(0.5, length = numpre), 
	Sigma = covmat.cases)
covmat.conts <- GenerateCliquesCovariance(numpre/2, 2, log(1.5)/2)$Sigma
pre.conts <- mvrnorm(numconts, mu = rep(0.5, length = numpre), 
	Sigma = covmat.conts)

## Generate facts for cases and controls and its parameter

facts.cases <- round(rnorm(numcases, 200, 60))
facts.conts <- round(rweibull(numconts, 1, 1)*100)
a0 <- log(1.5)/100

## Generate logit for each patient

# for cases
logitcases <- matrix(0, nrow = numcases, ncol = numdis)
for(eachcase in 1:numcases) {
	print(paste('Generating logit for case', eachcase, '...'))
	for(eachdis in 1:numdis) {
		logitcases[eachcase, eachdis] <- sum(betas.cases[eachdis, ]*
			pre.cases[eachcase, ]) + a0*facts.cases[eachcase] + beta0[eachdis]
	}
}
# for controls
logitconts <- matrix(0, nrow = numconts, ncol = numdis)
for(eachcont in 1:numconts) {
	print(paste('Generating logit for control', eachcont, '...'))
	for(eachdis in 1:numdis) {
		logitconts[eachcont, eachdis] <- sum(betas.controls[eachdis, ]*
			pre.conts[eachcont, ]) + a0*facts.conts[eachcont] + beta0[eachdis]
	}
}

## Calculate probabilities of disease from logit

expit <- function(x) 1/(1+exp(-x))
Pcases <- expit(logitcases)
Pconts <- expit(logitconts)

## Generate disease variables (binary) from probabilities

diCases <- matrix(0, nrow = numcases, ncol = numdis)
diConts <- matrix(0, nrow = numconts, ncol = numdis)

for(row in 1:nrow(diCases)) {
	print(paste('Generating disease profile for case', row, '...'))
	for(col in 1:ncol(diCases)) {
		pDis <- Pcases[row, col]
		diCases[row, col] <- rbinom(1, 1, pDis)
	}
}

for(row in 1:nrow(diConts)) {
	print(paste('Generating disease profile for control', row, '...'))
	for(col in 1:ncol(diConts)) {
		pDis <- Pconts[row, col]
		diConts[row, col] <- rbinom(1, 1, pDis)
	}
}

diCases <- cbind(rep(1, length=nrow(diCases)), diCases, facts.cases)
diConts <- cbind(rep(0, length=nrow(diConts)), diConts, facts.conts)

varNames <- character(808)
varNames[808] <- "FACTS"
for(eachVar in 1:806) {
	varNames[eachVar+1] <- paste("D", eachVar, sep="")
}
varNames[1] <- "EXPOSURE"

colnames(diCases) <- varNames
colnames(diConts) <- varNames

diCases <- as.data.frame(diCases)
diConts <- as.data.frame(diConts)
sampleConts <- diConts[sample(nrow(diConts), numcases), ]
mergedData <- rbind(diCases, sampleConts)

coefs.nofacts <- numeric(806)
coefs.facts <- numeric(806)
for(i in 1:806) {
	fit <- glm(mergedData[, i+1] ~ mergedData[, 1]  + mergedData[, 808], 
		family=binomial)
	fit2 <- glm(mergedData[, i+1] ~ mergedData[, 1], family=binomial)
	coefs.facts[i] <- coef(fit)[2]
	coefs.nofacts[i] <- coef(fit2)[2]
}

odds.nofacts <- exp(coefs.nofacts)
odds.facts <- exp(coefs.facts)
trueOdds.nofacts <- odds.nofacts[odds.nofacts < 5]
trueOdds.facts <- odds.facts[odds.facts < 5]

hist(trueOdds.facts, col=rgb(1,0,0,1/4), prob=T, breaks=20, main="Distribution 
	of Odds Ratios", xlab="Odds Ratios")
lines(density(trueOdds.facts), col=rgb(1,0,0,1), lwd=3)
hist(trueOdds.nofacts, col=rgb(0,0,1,1/4), prob=T, breaks=20, add=T)
lines(density(trueOdds.nofacts), col=rgb(0,0,1,1), lwd=3)

diCases.sorted <- diCases[order(diCases$FACTS), ]
diConts.sorted <- diConts[order(diConts$FACTS), ]

## Conditional logistic regression

## First, determine quantiles for bins
quantBins <- quantile(diCases$FACTS, prob=seq(0,1,.2))

## Bin 1 will be patients with facts between 0th and 20th percentile
print("Generating bin 1 ...")
logicalBin1 <- diConts$FACTS >= quantBins[1] & diConts$FACTS <= quantBins[2]
indexBin1 <- which(logicalBin1)
indexBin1 <- sample(indexBin1, numcases/5)
Bin1 <- data.frame(stringsAsFactors=FALSE) 
for(eachIndex1 in indexBin1) {
	Bin1 <- rbind.fill(Bin1, diConts[eachIndex1, ])
}

## Bin 2 will be patients with facts between 20th and 40th percentile
print("Generating bin 2 ...")
logicalBin2 <- diConts$FACTS >= quantBins[2] & diConts$FACTS <= quantBins[3]
indexBin2 <- which(logicalBin2)
indexBin2 <- sample(indexBin2, numcases/5)
Bin2 <- data.frame(stringsAsFactors=FALSE)
for(eachIndex2 in indexBin2) {
	Bin2 <- rbind.fill(Bin2, diConts[eachIndex2, ])
}

## Bin 3 will be patients with facts between 40th and 60th percentile
print("Generating bin 3 ...")
logicalBin3 <- diConts$FACTS >= quantBins[3] & diConts$FACTS <= quantBins[4]
indexBin3 <- which(logicalBin3)
indexBin3 <- sample(indexBin3, numcases/5)
Bin3 <- data.frame(stringsAsFactors=FALSE) 
for(eachIndex3 in indexBin3) {
	Bin3 <- rbind.fill(Bin3, diConts[eachIndex3, ])
}

## Bin 4 will be patients with facts between 60th and 80th percentile
print("Generating bin 4 ...")
logicalBin4 <- diConts$FACTS >= quantBins[4] & diConts$FACTS <= quantBins[5]
indexBin4 <- which(logicalBin4)
indexBin4 <- sample(indexBin4, numcases/5)
Bin4 <- data.frame(stringsAsFactors=FALSE) 
for(eachIndex4 in indexBin4) {
	Bin4 <- rbind.fill(Bin4, diConts[eachIndex4, ])
}

## Bin 5 will be patients with facts between 80th and 100th percentile
print("Generating bin 5 ...")
logicalBin5 <- diConts$FACTS >= quantBins[5] & diConts$FACTS <= quantBins[6]
indexBin5 <- which(logicalBin5)
indexBin5 <- sample(indexBin5, numcases/5)
Bin5 <- data.frame(stringsAsFactors=FALSE) 
for(eachIndex5 in indexBin5) {
	Bin5 <- rbind.fill(Bin5, diConts[eachIndex5, ])
}

binConts <- rbind(Bin1, Bin2, Bin3, Bin4, Bin5)
binConts <- binConts[order(binConts$FACTS), ]

## Combine binned controls with binned cases
binSet <- rbind(diCases.sorted, binConts)

## Use conditional logistic regression
condCoefs <- numeric(806)
for(i in 1:806) {
	condFit <- clogit(binSet[, i+1] ~ binSet[, 1] + binSet[, 808] + strata(rep(1:numcases, 2)))
	condCoefs[i] <- coef(condFit)[1]
}
condOdds <- exp(condCoefs)

hist(condOdds, col=rgb(0,1,0,1/4), prob=T, breaks=20, add=T)
lines(density(condOdds), col=rgb(0,1,0,1), lwd=3)

## Runtime calculation
runtime <- proc.time()-ptm
print('Calculating runtime...')
print(runtime)