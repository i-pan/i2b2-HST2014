# Load packages
library(spcov)
library(MASS)

# Generate indicator variables
indicators <- matrix(0, nrow = 806, ncol = 100)
for(eachrow in 1:806) {
	indicators[eachrow, ] <- rbinom(100, 1, 0.1)
}

# Generate all betas
all_betas <- rnorm(100, mean = 0, sd = log(1.5)/2)

# Generate case betas
case_betas <- matrix(0, nrow = 806, ncol = 100) 
for(eachrow in 1:806) {
	case_betas[eachrow, ] <- indicators[eachrow, ]*all_betas
}

# Generate control betas 
control_betas <- matrix(0, nrow = 806, ncol = 100)
for(eachrow in 1:50) {
	for(eachelement in 1:100) {
		if(case_betas[eachrow, eachelement] != 0) {
			control_betas[eachrow, eachelement] <- case_betas[eachrow, eachelement] - 0.05
		}
	}
}

# Generate values for random variables X_1 through X_100 for cases
howmanycases <- 100
cov_matrix <- GenerateCliquesCovariance(howmanycases/2, 2, 0.5)$Sigma
X_mat <- mvrnorm(howmanycases, mu = rep(0, length = 100), Sigma = cov_matrix)

# Generate logits for cases
logitCASES <- matrix(0, nrow = howmanycases, ncol = 806)
for(eachcase in 1:howmanycases) {
	for(eachrow in 1:806) {
		logitCASES[eachcase, eachrow] <- sum(case_betas[eachrow, ]*X_mat[eachcase, ]) 
	}
}

# Generate probabilities for cases
expit <- function(x) exp(x)/(1+exp(x))
Pcases <- expit(logitCASES)