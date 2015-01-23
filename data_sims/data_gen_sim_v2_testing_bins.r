# runtime calculation
ptm <- proc.time()

# load packages and set seed
library(MASS)
library(plyr)
library(survival)
library(spcov)
set.seed(10)

## generate indicator variables to indicate which predictors will apply to 
## which diseases

p <- 0.05  # insert probability here
numpre <- 100  # insert number of predictors
numdis <- 806  # insert number of diseases
indicators <- matrix(rbinom(numpre*numdis, 1, p), nrow = numdis, ncol = numpre)
# for(eachrow in 1:numdis) {
	# indicators[eachrow, ] = rbinom(numpre, 1, p)
# }

## generate beta parameters for each predictor for each disease 

numcomb <- 50  # number of comorbidities

# for cases
betas.cases <- matrix(rnorm(numdis*numpre, mean = 0, sd = log(1.5)/2), nrow = 
	numpre, ncol = numdis)
for(eachbeta in 1:numdis) {
	betas.cases[, eachbeta] <- betas.cases[, eachbeta]*indicators[eachbeta, ]
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

## generate predictors for each patient

beta0 <- runif(806, -2, 0)
numcases <- 500  # insert number of cases
numconts <- 5000  # insert number of controls
covmat.cases <- GenerateCliquesCovariance(numpre/2, 2, log(1.5)/2)$Sigma
pre.cases <- mvrnorm(numcases, mu = rep(0.5, length = numpre), 
	Sigma = covmat.cases)
covmat.conts <- GenerateCliquesCovariance(numpre/2, 2, log(1.5)/2)$Sigma
pre.conts <- mvrnorm(numconts, mu = rep(0.5, length = numpre), 
	Sigma = covmat.conts)

# generate facts for cases and controls and its parameter
facts.cases <- round(rnorm(numcases, 200, 60))
facts.conts <- round(rweibull(numconts, 1, 1)*100)
a0 <- log(1.5)/100

# combine facts with predictors
pre.cases <- cbind(pre.cases, facts.cases)
pre.conts <- cbind(pre.conts, facts.conts)

# combine a0 with betas
betas.cases <- rbind(betas.cases, rep(a0, length=numdis))
betas.controls <- rbind(betas.controls, rep(a0, length=numdis))

## generate logit for each patient

# for cases
logitcases <- pre.cases%*%betas.cases
# for(eachcase in 1:numcases) {
	# print(paste('Generating logit for case', eachcase, '...'))
	# for(eachdis in 1:numdis) {
		# logitcases[eachcase, eachdis] <- sum(betas.cases[eachdis, ]*
			# pre.cases[eachcase, ]) + a0*facts.cases[eachcase] + beta0[eachdis]
	# }
# }
# for controls
logitconts <- pre.conts%*%betas.controls
# for(eachcont in 1:numconts) {
	# print(paste('Generating logit for control', eachcont, '...'))
	# for(eachdis in 1:numdis) {
		# logitconts[eachcont, eachdis] <- sum(betas.controls[eachdis, ]*
			# pre.conts[eachcont, ]) + a0*facts.conts[eachcont] + beta0[eachdis]
	# }
# }

## calculate probabilities of disease from logit

expit <- function(x) 1/(1+exp(-x))
Pcases <- expit(logitcases)
Pconts <- expit(logitconts)

## generate disease variables (binary) from probabilities

# diCases <- matrix(0, nrow = numcases, ncol = numdis)
# diConts <- matrix(0, nrow = numconts, ncol = numdis)

# for(row in 1:nrow(diCases)) {
	# print(paste('Generating disease profile for case', row, '...'))
	# for(col in 1:ncol(diCases)) {
		# pDis <- Pcases[row, col]
		# diCases[row, col] <- rbinom(1, 1, pDis)
	# }
# }

# for(row in 1:nrow(diConts)) {
	# print(paste('Generating disease profile for control', row, '...'))
	# for(col in 1:ncol(diConts)) {
		# pDis <- Pconts[row, col]
		# diConts[row, col] <- rbinom(1, 1, pDis)
	# }
# }

GetBinomMat2 <- function(prob.mat)
    matrix(sweep(prob.mat, 1:2, prob.mat, rbinom, size = 1),
        ncol = ncol(prob.mat))

print("Generating disease profile for cases ...")		
diCases <- GetBinomMat2(Pcases)
print("Generating disease profile for controls ...")
diConts <- GetBinomMat2(Pconts)

# add exposure variable and facts to the dataset
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

## logistic regression model

# randomly sample controls
sampleConts <- diConts[sample(nrow(diConts), numcases), ]
mergedData <- rbind(diCases, sampleConts)

coefs.nofacts <- numeric(806)
coefs.facts <- numeric(806)
pvalues.nofacts <- numeric(806)
pvalues.facts <- numeric(806)

# compute regression fits
for(i in 1:806) {
	print(paste("Fitting logistic regression for disease", i, "..."))
	fit <- glm(mergedData[, i+1] ~ mergedData[, 1]  + mergedData[, 808], 
		family=binomial) # facts in model 
	fit2 <- glm(mergedData[, i+1] ~ mergedData[, 1], family=binomial) # no facts
	pvalues.facts[i] <- summary(fit)$coefficients[, 4][2]
	pvalues.nofacts[i] <- summary(fit2)$coefficients[, 4][2]
	coefs.facts[i] <- coef(fit)[2]
	coefs.nofacts[i] <- coef(fit2)[2]
}

# compute odds ratios
odds.nofacts <- exp(coefs.nofacts)
odds.facts <- exp(coefs.facts)
trueOdds.nofacts <- odds.nofacts[odds.nofacts < 5]
trueOdds.facts <- odds.facts[odds.facts < 5]

# graph plots 
hist(trueOdds.facts, col=rgb(1,0,0,1/4), prob=T, breaks=20, main="Distribution 
	of Odds Ratios", xlab="Odds Ratios", xlim=c(0.5,2.5), ylim=c(0,3.5))
lines(density(trueOdds.facts), col=rgb(1,0,0,1), lwd=3)
hist(trueOdds.nofacts, col=rgb(0,0,1,1/4), prob=T, breaks=20, add=T)
lines(density(trueOdds.nofacts), col=rgb(0,0,1,1), lwd=3)

# sort datasets
diCases.sorted <- diCases[order(diCases$FACTS), ]
diConts.sorted <- diConts[order(diConts$FACTS), ]

## conditional logistic regression

# generate bins
numBins <- 5 # insert number of bins
quantBins <- quantile(diCases$FACTS, prob=seq(0,1,(1/numBins)))
logicalBins <- vector("list", numBins)
indexBins <- vector("list", numBins)

binConts <- matrix(0, nrow=numcases, ncol=ncol(diCases))

for(i in 1:numBins) {
	print(paste("Generating bin", i, "..."))
	logicalBins[[i]] <- diConts$FACTS >= quantBins[i] & diConts$FACTS <= 
		quantBins[i+1]
	indexBins[[i]] <- which(logicalBins[[i]])
	indexBins[[i]] <- sample(indexBins[[i]], numcases/numBins)
	for(eachIndex in indexBins[[i]]) {
		binConts <- rbind(binConts, diConts[eachIndex, ])
	}
}

# sort bins
binConts <- as.data.frame(binConts)
binConts <- binConts[order(binConts$FACTS), ]

# combine binned controls with binned cases
binSet <- rbind(diCases.sorted, binConts)

# use conditional logistic regression
condCoefs <- numeric(806)
condPvalues <- numeric(806)
for(i in 1:806) {
	print(paste("Fitting conditional logistic regression for disease", i, 
		"..."))
	condFit <- clogit(binSet[, i+1] ~ binSet[, 1] + binSet[, 808] + 
		strata(rep(1:numcases, 2)))
	condPvalues[i] <- summary(condFit)$logtest[3]
	condCoefs[i] <- coef(condFit)[1]
}
# compute odds ratios
condOdds <- exp(condCoefs)

# plot odds ratios
hist(condOdds, col=rgb(0,1,0,1/4), prob=T, breaks=20, add=T)
lines(density(condOdds), col=rgb(0,1,0,1), lwd=3)

# runtime calculation
runtime <- proc.time()-ptm
print('Calculating runtime...')
print(runtime)