##### RUNTIME #####
ptm <- proc.time()

##### LOAD LIBRARIES #####

library(survival)
library(locfdr)
set.seed(10)

##### LOAD DATA #####

dem <- read.csv('dem.csv') # demographics file
dx <- read.csv('dx.csv') # diagnostics file 

# phewas groupings, obtained from:
# http://knowledgemap.mc.vanderbilt.edu/research/content/phewas
phewas.code <- read.table("phewas_code.txt", header = T) 
phewas.code <- as.matrix(phewas.code)
phewas.code[, 'phewas_code'] <- gsub(" ", "", phewas.code[, 'phewas_code'])


##### ORGANIZE DATA #####

# split into cases and controls
case <- dem[dem$grp == "Case", ]
cont <- dem[dem$grp == "Control", ]
case <- cbind(case, status = rep(1, length = nrow(case)))
cont <- cbind(cont, status = rep(0, length = nrow(cont)))

# create age bins (X-year bins) i.e. 10-year bins = 0-10, 10-20, etc.
year.bins <- 10
case <- cbind(case, 
              age_bin = as.numeric(cut(case$age, 
							                         breaks = seq(0, 11*year.bins, 
																			              by = year.bins), 
																			 labels = c(seq(0, year.bins, by = 1)), 
	                                     right = TRUE)))
cont <- cbind(cont, 
              age_bin = as.numeric(cut(cont$age, 
							                         breaks = seq(0, 11*year.bins, 
																						        by = year.bins), 
																			 labels = c(seq(0, year.bins, by = 1)), 
	                                     right = TRUE)))

# create facts bins 
fact.bins <- 20 
fact.qtile <- c(quantile(case$num_facts, probs = seq(0, 1, by = 1/fact.bins)))
case <- cbind(case, 
              fact_bin = as.numeric(cut(case$num_facts, 
							                          breaks = fact.qtile, 
																				labels = seq(1, fact.bins, 
																							       by = 1), 
																				include.lowest = TRUE)))
cont <- cbind(cont, 
							fact_bin = as.numeric(cut(cont$num_facts, 
							                          breaks = fact.qtile, 
																				labels = seq(1, fact.bins, 
																							       by = 1), 
																				include.lowest = TRUE)))
	
# sort by gender, race, age, facts
case <- case[order(case$gender, case$race, case$age_bin, case$fact_bin, 
	                 case$num_facts), ]
cont <- cont[order(cont$gender, cont$race, cont$age_bin, cont$fact_bin, 
	                 cont$num_facts), ]


##### MATCH CASES TO CONTROLS (1-TO-1) #####

index <- vector("list", nrow(case))
index.sample <- vector("numeric", nrow(case))
for(each.case in 1:nrow(case)) {
	print(paste("Matching case", each.case, 
							"based on gender, race, exact age, exact facts ..."))
  logic.case <- cont$gender == case[each.case, ]$gender &
                cont$race == case[each.case, ]$race &
                cont$age == case[each.case, ]$age & 
		            cont$num_facts == case[each.case, ]$num_facts
  if(!(TRUE %in% logic.case)) {
	  print("Could not find exact match ...")
		print(paste("Matching case", each.case, 
								"based on gender, race, age bin, exact facts ..."))
    logic.case <- cont$gender == case[each.case, ]$gender &
                  cont$race == case[each.case, ]$race &
                  cont$age_bin == case[each.case, ]$age_bin &
                  cont$num_facts == case[each.case, ]$num_facts
    if(!(TRUE %in% logic.case)) {
		  print("Could not find exact match ...")
			print(paste("Matching case", each.case, 
									"based on gender, race, age bin, fact bin ..."))
      logic.case <- cont$gender == case[each.case, ]$gender &
                    cont$race == case[each.case, ]$race &
                    cont$age_bin == case[each.case, ]$age_bin &
                    cont$fact_bin == case[each.case, ]$fact_bin
      if(!(TRUE %in% logic.case)) {
			  print("Could not find exact match ...")
				print(paste("Matching case", each.case, 
										"based on gender, age bin, fact bin ..."))
        logic.case <- cont$gender == case[each.case, ]$gender &
                      cont$age_bin == case[each.case, ]$race &
                      cont$fact_bin == case[each.case, ]$fact_bin
        if(!(TRUE %in% logic.case)) {
				print("Could not find exact match ...")
				print(paste("Matching case", each.case, 
										"based on age bin, fact bin ..."))
          logic.case <- cont$age_bin == case[each.case, ]$age_bin &
                        cont$fact_bin == case[each.case, ]$fact_bin
          if(!(TRUE %in% logic.case)) {
					  print("Could not find exact match ...")
						print(paste("Matching case", each.case, "based on fact bin ..."))
            logic.case <- cont$fact_bin == case[each.case, ]$fact_bin
          }
        } 
      }
    }
  }
  index[[each.case]] <- which(logic.case)
  if(length(index[[each.case]]) == 1) {
    index.sample[each.case] <- index[[each.case]]
  } else {
  index.sample[each.case] <- sample(index[[each.case]], 1)
  }
}

cont.sample <- cont[index.sample, ]

# combine cases and controls
case.cont <- rbind(case, cont.sample)

##### OPTIONAL CODE TO PheWAS CODES BASED ON ROLLUP BOOLEAN VALUE #####

# for(each.code in 1:nrow(phewas.code)) {
	# if(phewas.code[each.code, ]$rollup_bool == 1) {
		# phewas.code[each.code, ]$phewas_code <- 
			# trunc(phewas.code[each.code, ]$phewas_code)
	# }
# }

##### ASSIGN PheWAS CODES TO RESPECTIVE ICD-9 CODES #####

case.cont.ID <- case.cont$patient_num
case.cont.dx <- dx[which(dx$patient_num %in% case.cont.ID), ]
case.cont.dx <- as.matrix(case.cont.dx)
case.cont.dx <- cbind(case.cont.dx, 
											phewas_code = rep(0, length = nrow(case.cont.dx)),
											exclude_range = rep(0, length = nrow(case.cont.dx)))
											
tmp.dx.index <- which(case.cont.dx[, 2] %in% phewas.code[, 'icd9'])
case.cont.dx <- case.cont.dx[tmp.dx.index, ]
store.phewas.code <- character(length=nrow(case.cont.dx))
store.exclude.range <- character(length=nrow(case.cont.dx))

icd9.codes <- data.frame(table(case.cont.dx[, 2]))
icd9.codes <- as.character(icd9.codes$Var1)
store.phewas.code <- character(length=nrow(case.cont.dx))
store.exclude.range <- character(length=nrow(case.cont.dx))
for(i in 1:length(icd9.codes)) {
  print(paste("Retrieving PheWAS code for ICD-9 code", i, "of", length(icd9.codes)))
  temp.phewas.code <- phewas.code[phewas.code[, 'icd9'] == icd9.codes[i], 'phewas_code']
  temp.exclude.range <- phewas.code[phewas.code[, 'icd9'] == icd9.codes[i], 'exclude_range']
  store.index <- which(case.cont.dx[, 2] == icd9.codes[i])
  store.phewas.code[store.index] <- temp.phewas.code
  store.exclude.range[store.index] <- temp.exclude.range
}

case.cont.dx[, 3] <- store.phewas.code
case.cont.dx[, 4] <- store.exclude.range

# index <- 1
# for(each.code in case.cont.dx[, 2]) {
  # print(paste("Retrieving PheWAS code for row", index, "of row", nrow(case.cont.dx), "..."))
  # store.phewas.code[index] <- phewas.code[phewas.code[, 'icd9'] == each.code, 'phewas_code']
  # store.exclude.range[index] <- phewas.code[phewas.code[, 'icd9'] == each.code, 'exclude_range']
  # index <- index + 1
# }

# case.cont.dx[, 3] <- store.phewas.code
# case.cont.dx[, 4] <- store.exclude.range

##### CREATE A DISEASE MATRIX BASED ON PheWAS GROUPINGS ##### 
##### IF ELEMENT i,j = 1, THEN CASE i HAS DISEASE j, 0 OTHERWISE #####

case.cont.dx <- case.cont.dx[complete.cases(case.cont.dx), ]
case.cont.dx <- case.cont.dx[case.cont.dx[, 4] != "555-564.99", ]
case.cont.dx[, 1] <- as.integer(case.cont.dx[, 1])
dis.names <- as.character(as.data.frame(table(case.cont.dx[, 3]))$Var1)
dis.matrix <- matrix(0, nrow = nrow(case.cont), ncol = length(dis.names))
colnames(dis.matrix) <- dis.names
dis.matrix <- cbind(patient_num = case.cont.ID, dis.matrix)

for(i in 1:nrow(case.cont.dx)) {
  print(paste("Operating on row", i, "of row", nrow(case.cont.dx), "..."))
  tmp.ID <- as.integer(case.cont.dx[i, 1])
  tmp.phewas <- case.cont.dx[i, 3]
  dis.matrix[dis.matrix[, 1] == tmp.ID, tmp.phewas] <- 1  
}

dis.matrix.case <- dis.matrix[1:nrow(case), ]
dis.matrix.cont <- dis.matrix[(nrow(case)+1):(nrow(case)+nrow(cont)), ]

dis.incid.case <- numeric()
dis.incid.cont <- numeric()

dis.matrix.col <- ncol(dis.matrix)-1 
for(each.dis in 1:dis.matrix.col) {
	dis.incid.case[each.dis] <- sum(dis.matrix.case[, each.dis+1])
	dis.incid.cont[each.dis] <- sum(dis.matrix.cont[, each.dis+1])
}

##### BEGIN COMPUTING ODDS RATIOS #####

# establish threshold for disease incidence
thres <- 5
reach.thres.case <- which(dis.incid.case >= thres) + 1
reach.thres.cont <- which(dis.incid.cont >= thres) + 1
reach.thres <- intersect(reach.thres.case, reach.thres.cont)
dis.matrix.thres <- dis.matrix[, reach.thres]

# drop Crohn's disease, if present
dis.matrix.thres <- dis.matrix.thres[, !(colnames(dis.matrix.thres) %in% "555.10")]

case.cont.data <- cbind(case.cont, dis.matrix.thres)

# compute logistic regression
dis.matrix.thres.col <- ncol(dis.matrix.thres)
coefs.nofacts <- numeric(dis.matrix.thres.col)
coefs.facts <- numeric(dis.matrix.thres.col)
pvalues.nofacts <- numeric(dis.matrix.thres.col)
pvalues.facts <- numeric(dis.matrix.thres.col)
zvalues.nofacts <- numeric(dis.matrix.thres.col)
zvalues.facts <- numeric(dis.matrix.thres.col)

dis.matrix.thres <- cbind(case.cont.data$status, case.cont.data$num_facts,
													dis.matrix.thres)
													
for(i in 1:dis.matrix.thres.col) {
	print(paste("Fitting logistic regression for disease", i, "..."))
	fit <- glm(dis.matrix.thres[, i+2] ~ dis.matrix.thres[, 1] + 
																			 dis.matrix.thres[, 2], 
						 family=binomial) # facts in model 
	fit2 <- glm(dis.matrix.thres[, i+2] ~ dis.matrix.thres[, 1],  
						  family=binomial) # no facts
	pvalues.facts[i] <- summary(fit)$coefficients[, 4][2]
	pvalues.nofacts[i] <- summary(fit2)$coefficients[, 4][2]
	coefs.facts[i] <- coef(fit)[2]
	coefs.nofacts[i] <- coef(fit2)[2]
	zvalues.facts[i] <- coef(summary(fit))[, "z value"][2]
	zvalues.nofacts[i] <- coef(summary(fit2))[, "z value"][2]
}
# compute odds ratios
odds.facts <- exp(coefs.facts)
odds.nofacts <- exp(coefs.nofacts)

# write out values to a table
temp.df <- data.frame(odds_facts = odds.facts, odds_nofacts = odds.nofacts, 
											coefs_facts = coefs.facts, coefs_nofacts = coefs.nofacts,
											pvalues_facts = pvalues.facts, pvalues_nofacts = pvalues.nofacts,
											zvalues_facts = zvalues.facts, zvalues_nofacts = zvalues.nofacts)
write.csv(temp.df, "values_matched.csv")

fdr <- locfdr(temp.df$zvalues_nofacts, bre = 100)$fdr
Efdr <- locfdr(temp.df$zvalues_nofacts, bre = 100)$Efdr
fdr.sig.index <- which(fdr < Efdr[3])
odds.sig.index <- which(temp.df$odds_nofacts > 1)
sig.index <- intersect(fdr.sig.index, odds.sig.index)
all.disease.names <- colnames(dis.matrix.thres)[3:ncol(dis.matrix.thres)]
index <- 1
comorbs <- character()
for(each.disease in all.disease.names[sig.index]) {
	comorbs[index] <- phewas.code[phewas.code[, 'phewas_code'] == each.disease, 'phewas_string'][1]
	
	index <- index+1
}
comorbs.df <- data.frame(Comorbidities = comorbs, Odds_Ratios = temp.df$odds_nofacts[sig.index],
												 PheWAS_Codes = all.disease.names[sig.index])
write.csv(comorbs.df, file="comorbidities_odds_ratios_long_nofacts_5.csv")

# conditional logistic regression
cond.coefs <- numeric(dis.matrix.thres.col)
cond.pvalues <- numeric(dis.matrix.thres.col)
for(i in 1:dis.matrix.thres.col) {
	print(paste("Fitting conditional logistic regression for disease", i, 
		"..."))
	cond.fit <- clogit(dis.matrix.thres[, i+2] ~ dis.matrix.thres[, 1] + dis.matrix.thres[, 2] + strata(case.cont.data$fact_bin))
	cond.pvalues[i] <- summary(cond.fit)$logtest[3]
	cond.coefs[i] <- coef(cond.fit)[1]
	cond.zvalues[i] <- 
}
# compute odds ratios
cond.odds <- exp(cond.coefs)

# plot histograms and densities
png.file.name <- paste("Odds_Ratios_", thres, ".png", sep = "")
png(png.file.name)
odds.density <- density(odds.facts)
hist(odds.facts, col=rgb(1,0,0,1/4), 
     prob=T, 
		 breaks=20, 
     main="Distribution of Odds Ratios", 
		 xlab="Odds Ratios", 
		 xaxt="n",
		 xlim=c(0,4),
		 ylim=c(0, max(odds.density$y)*1.5))
axis(side=1, at=seq(0, max(odds.facts)+1, by = 0.5))
lines(odds.density, col=rgb(1,0,0,1), lwd=3)
hist(odds.nofacts, col=rgb(0,0,1,1/4), prob=T, breaks=20, add=T)
lines(density(odds.nofacts), col=rgb(0,0,1,1), lwd=3)
hist(cond.odds, col=rgb(0,1,0,1/4), prob=T, breaks=20, add=T)
lines(density(cond.odds), col=rgb(0,1,0,1), lwd=3)
garbage <- dev.off()

##### IDENTIFY COMORBIDITIES OBTAINED FROM EACH MODEL #####

# identify comorbidities 
alpha <- 0.05 # set value for alpha

cond.dis.pvalues <- which(cond.pvalues < alpha)
dis.facts.pvalues <- which(pvalues.facts < alpha)
dis.nofacts.pvalues <- which(pvalues.nofacts < alpha)

cond.dis.greater.1 <- which(cond.odds > 1)
dis.facts.greater.1 <- which(odds.facts > 1)
dis.nofacts.greater.1 <- which(odds.nofacts > 1)

cond.dis <- intersect(cond.dis.pvalues, cond.dis.greater.1)
dis.facts <- intersect(dis.facts.pvalues, dis.facts.greater.1)
dis.nofacts <- intersect(dis.nofacts.pvalues, dis.nofacts.greater.1)

dis.included.names <- colnames(case.cont.data)
most.comorbs <- max(length(cond.dis), length(dis.facts), length(dis.nofacts))

# create table for storage of comorbidities
comorb.matrix <- matrix(NA, nrow = most.comorbs, ncol = 6) 
colnames(comorb.matrix) <- c("CLR_Facts", "CLR_OR", "LR_Facts", "LR_Facts_OR", 
														 "LR_Simple", "LR_Simple_OR")

# comorbidities using conditional logistic regression
index <- 1
for(each.dis in cond.dis) {
	tmp.dis.name <- dis.included.names[11+each.dis]
	tmp.comorb <- phewas.code[phewas.code[, 'phewas_code'] == tmp.dis.name, 
												'phewas_string'][1]
	comorb.matrix[index, 1] <- as.character(tmp.comorb)
	comorb.matrix[index, 2] <- round(cond.odds[each.dis], digits = 2)
	index <- index + 1											
}

# comorbidities using logistic regression with facts in model 
index <- 1
for(each.dis in dis.facts) {
	tmp.dis.name <- dis.included.names[11+each.dis]
	tmp.comorb <- phewas.code[phewas.code[, 'phewas_code'] == tmp.dis.name, 
												'phewas_string'][1]
	comorb.matrix[index, 3] <- as.character(tmp.comorb)
	comorb.matrix[index, 4] <- round(odds.facts[each.dis], digits = 2)
	index <- index + 1	
}

# comorbidities using simple logistic regression
index <- 1
for(each.dis in dis.nofacts) {
	tmp.dis.name <- dis.included.names[11+each.dis]
	tmp.comorb <- phewas.code[phewas.code[, 'phewas_code'] == tmp.dis.name, 
												'phewas_string'][1]
	comorb.matrix[index, 5] <- as.character(tmp.comorb)
	comorb.matrix[index, 6] <- round(odds.nofacts[each.dis], digits = 2)
	index <- index + 1	
}

# write out table
comorb.df <- as.data.frame(comorb.matrix)
write.csv(comorb.df, file = "comorbidities.csv")

##### RUNTIME #####
runtime <- proc.time() - ptm
print("Runtime:")
print(runtime)