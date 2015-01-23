# runtime calculation
ptm <- proc.time()

# load packages and set seed
library(MASS)
library(plyr)
library(survival)
library(spcov)
library(locfdr)
set.seed(10)

## generate indicator variables to indicate which predictors will apply to 
## which diseases

p <- 0.05  # insert probability here
numpre <- 100  # insert number of predictors
numdis <- 806  # insert number of diseases
indicators <- matrix(rbinom(numpre*numdis, 1, p), nrow = numpre, ncol = numdis)

## generate beta parameters for each predictor for each disease 

numcomb <- 112  # number of comorbidities

betas.cases <- matrix(rnorm(numdis*numpre, mean = 0, sd = log(1.5)/2), nrow = 
	numpre, ncol = numdis)
betas.controls <- betas.cases
betas.controls[, 1:numcomb] <- betas.controls[, 1:numcomb]-0.2

betas.cases <- betas.cases*indicators
betas.controls <- betas.controls*indicators

## generate predictors for each patient

beta0 <- runif(806, -2, 0)
numcases <- 2500  # insert number of cases
numconts <- 25000  # insert number of controls
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

# for controls
logitconts <- pre.conts%*%betas.controls


## calculate probabilities of disease from logit

expit <- function(x) 1/(1+exp(-x))
Pcases <- expit(logitcases)
Pconts <- expit(logitconts)

## generate disease variables (binary) from probabilities

getDisease <- function(prob.mat) 
    matrix(rbinom(prob.mat, 1, prob.mat), ncol = ncol(prob.mat))

print("Generating disease profile for cases ...")
diCases <- getDisease(Pcases)
print("Generating disease profile for controls ...")
diConts <- getDisease(Pconts)

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

## logistic regression model

# randomly sample controls
sampleConts <- diConts[sample(nrow(diConts), numcases), ]
mergedData <- rbind(diCases, sampleConts)

coefs.nofacts <- numeric(806)
coefs.facts <- numeric(806)
pvalues.nofacts <- numeric(806)
pvalues.facts <- numeric(806)
zvalues.facts <- numeric(806)
zvalues.nofacts <- numeric(806)

# compute regression fits
for(i in 1:806) {
	print(paste("Fitting logistic regression for disease", i, "..."))
	fit <- glm(mergedData[, i+1] ~ mergedData[, 1]  + mergedData[, 808], 
		family=binomial) # facts in model 
	fit2 <- glm(mergedData[, i+1] ~ mergedData[, 1], family=binomial) # no facts
	pvalues.facts[i] <- summary(fit)$coefficients[, 4][2]
	pvalues.nofacts[i] <- summary(fit2)$coefficients[, 4][2]
	zvalues.facts[i] <- coef(summary(fit))[, "z value"][2]
	zvalues.nofacts[i] <- coef(summary(fit2))[, "z value"][2]
	coefs.facts[i] <- coef(fit)[2]
	coefs.nofacts[i] <- coef(fit2)[2]
}

# compute odds ratios
odds.nofacts <- exp(coefs.nofacts)
odds.facts <- exp(coefs.facts)
trueOdds.nofacts <- odds.nofacts[odds.nofacts < 5]
trueOdds.facts <- odds.facts[odds.facts < 5]

# graph plots
tOfDensity <- density(trueOdds.facts)
hist(trueOdds.facts, col=rgb(1,0,0,1/4), prob=T, breaks=20, main="Distribution 
	of Odds Ratios", xlab="Odds Ratios", xlim=c(0.5,2.5), ylim=c(0,
	max(tOfDensity$y)*1.5))
lines(tOfDensity, col=rgb(1,0,0,1), lwd=3)
hist(trueOdds.nofacts, col=rgb(0,0,1,1/4), prob=T, breaks=20, add=T)
lines(density(trueOdds.nofacts), col=rgb(0,0,1,1), lwd=3)

# sort datasets
diCases.sorted <- diCases[order(diCases[, 808]), ]
diConts.sorted <- diConts[order(diConts[, 808]), ]

## conditional logistic regression

# generate bins
numBins <- 10 # insert number of bins
quantBins <- quantile(diCases[, 808], prob=seq(0,1,(1/numBins)))
logicalBins <- vector("list", numBins)
indexBins <- vector("list", numBins)

binConts <- matrix(, nrow=0, ncol=808)

# for(i in 1:numBins) {
	# print(paste("Generating bin", i, "..."))
	# logicalBins[[i]] <- diConts[, 808] >= quantBins[i] & diConts[, 808] <= 
		# quantBins[i+1]
	# indexBins[[i]] <- which(logicalBins[[i]])
	# indexBins[[i]] <- sample(indexBins[[i]], numcases/numBins)
	# for(eachIndex in indexBins[[i]]) {
		# binConts <- rbind(binConts, diConts[eachIndex, ])
	# }
# }

for(i in 1:numBins) {
	print(paste("Generating bin", i, "..."))
	logicalBins[[i]] <- diConts[, 808] >= quantBins[i] & diConts[, 808] <= 
		quantBins[i+1]
	indexBins[[i]] <- which(logicalBins[[i]])
	indexBins[[i]] <- sample(indexBins[[i]], numcases/numBins)
	binConts <- rbind(binConts, diConts[indexBins[[i]], ])
}
# sort bins
binConts <- binConts[order(binConts[, 808]), ]

# combine binned controls with binned cases
binSet <- rbind(diCases.sorted, binConts)

# use conditional logistic regression
zvalues.cond <- numeric(806)
condCoefs <- numeric(806)
condPvalues <- numeric(806)
for(i in 1:806) {
	print(paste("Fitting conditional logistic regression for disease", i, 
		"..."))
	condFit <- clogit(binSet[, i+1] ~ binSet[, 1] + binSet[, 808] + 
		strata(rep(1:numcases, 2)))
	condPvalues[i] <- summary(condFit)$logtest[3]
	condCoefs[i] <- coef(condFit)[1]
	zvalues.cond[i] <- coef(summary(condFit))[, 'z'][1]
}
# compute odds ratios
condOdds <- exp(condCoefs)

# estimate empirical null distribution
locfdr.facts <- locfdr(zvalues.facts, bre=50, pct0=0.2)
locfdr.nofacts <- locfdr(zvalues.nofacts, bre=50, pct0=0.2, nulltype=2)
locfdr.cond <- locfdr(zvalues.cond, bre=50, pct0=0.2)

MLEmean.facts <- locfdr.facts$fp0[3, 1]
MLEsd.facts <- locfdr.facts$fp0[3, 2]
MLEmean.nofacts <- locfdr.facts$fp0[3, 1]
MLEsd.nofacts <- locfdr.facts$fp0[3, 2]
MLEmean.cond <- locfdr.facts$fp0[3, 1]
MLEsd.cond <- locfdr.facts$fp0[3, 2]

# adjust p-values
adj.pvalues.facts <- 2*(1-pnorm(zvalues.facts, mean=MLEmean.facts, sd=MLEsd.facts))
adj.pvalues.nofacts <- 2*(1-pnorm(zvalues.nofacts, mean=MLEmean.nofacts, sd=MLEsd.nofacts))
adj.pvalues.cond <- 2*(1-pnorm(zvalues.cond, mean=MLEmean.cond, sd=MLEsd.cond))

# turn adj pvalues into fdr
adj.fdr.facts <- p.adjust(adj.pvalues.facts, method="fdr")
adj.fdr.nofacts <- p.adjust(adj.pvalues.nofacts, method="fdr")
adj.fdr.cond <- p.adjust(adj.pvalues.cond, method="fdr")

temp.df <- data.frame(disease = 1:806, odds_facts = odds.facts, pvalues_facts = pvalues.facts,
                      pvalues_nofacts = pvalues.nofacts, pvalues_cond = condPvalues,
											adj_pvalues_facts = adj.pvalues.facts, adj_pvalues_nofacts = adj.pvalues.nofacts,
											adj_pvalues_cond = adj.pvalues.cond, fdr_facts = adj.fdr.facts,
											fdr_nofacts = adj.fdr.nofacts, fdr_cond = adj.fdr.cond)

# plot odds ratios
hist(condOdds, col=rgb(0,1,0,1/4), prob=T, breaks=20, add=T)
lines(density(condOdds), col=rgb(0,1,0,1), lwd=3)

# runtime calculation
runtime <- proc.time()-ptm
print('Calculating runtime...')
print(runtime)