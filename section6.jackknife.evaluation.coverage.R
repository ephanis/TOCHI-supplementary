
###############################################################################################
# Using Bailly et al's dataset to generate a dataset to evaluate the coverage of the jackknife method for Fleiss' kappa
#
#
# author: Theophanis Tsandilas
##################################################################

rm(list=ls()) # Clean up R's memory

source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")

includes <- function(val, low, high){
    if(val <= high && val >= low) TRUE
    else FALSE
}

include <- function(vals, lows, highs){
    as.numeric(vals>=lows & vals<=highs)
}


# n is the number of participants that I randomly draw from the overall population. R = number or repetitions
coverage.estimation <- function(signs.full, value.all, R = 100, n = 20, statistic=fleiss.kappa, conflev = .99){
    
    counter <- 0
    
    for(r in 1:R){
        signs.sample <- signs.full[,sample(1:ncol(signs.full), size = n)]
        
        ci <- jack.CI.random.raters(signs.sample, statistic, conflev)
        
        if(includes(value.all, ci[2], ci[3])) counter <- counter + 1
        
    }
    
    counter / R
}


# n is the number of participants that I randomly draw from the overall population. R = number or repetitions
coverage.estimation.referent <- function(signs.full, kappas.all, R = 100, n = 20, statistic=fleiss.kappa, conflev = .99){
    L <- length(kappas.all)
    coverage <- rep(0, L)
    
    for(r in 1:R){
        signs.sample <- signs.full[,sample(1:ncol(signs.full), size = n)]
        
        kappas <- sapply(1:L, function(i) jack.CI.random.raters.fleiss.kappa.for.item(signs.sample, i, conflev))
        
        coverage <- coverage + include(kappas.all, kappas[2,], kappas[3,])
    }
    
    coverage / R
}


# L is the number of referents and l is the number of random comparisons to make
coverage.estimation.referents.groups <- function(signs.full, kappas.all, l = 10, size = 10, R = 100, n = 20, statistic=fleiss.kappa, conflev = .99){
    L <- length(kappas.all)
    estimates <- list()
    
    for(i in 1:l){
        indices <- sample(1:L, size=size)
        kappa <- mean(kappas.all[indices])
        
        # cat("Kappa = ", kappa,  " ( n =", length(indices), ")\n")
        
        counter <- 0
        for(j in 1:R){
            signs.sample <- signs.full[,sample(1:ncol(signs.full), size = n)]
            ci <- jack.CI.random.raters.fleiss.kappa.for.item(signs.sample, indices, conflev)
            if(includes(kappa, ci[2], ci[3])) counter <- counter + 1
        }
        # cat("CI = ", counter/R, "%\n")
        
        estimates[[i]] <- counter / R
    }
    
    estimates
}


# L is the number of referents and l is the number of random comparisons to make
coverage.estimation.referents.diff <- function(signs.full, kappas.all, l = 10, R = 100, n = 20, statistic=fleiss.kappa, conflev = .99){
    L <- length(kappas.all)
    estimates <- list()
    
    for(i in 1:l){
        # Create two random groups, each of a random size from 1 to 8
        size1 <- sample(1:8,size=1)
        size2 <- min(sample(1:8,size=1), L - size1)
        
        combined <- sample(1:L, size=size1+size2)
        indices1 <- combined[1:size1]
        indices2 <- combined[(size1+1):(size1+size2)]
        
        diff <- mean(kappas.all[indices1]) -  mean(kappas.all[indices2])
        
        counter <- 0
        for(j in 1:R){
            signs.sample <- signs.full[,sample(1:ncol(signs.full), size = n)]
            
            ci <- jack.CI.diff.random.raters.fleiss.kappa(signs.sample, indices1, indices2, conflevel = conflev)
            
            if(includes(diff, ci[2], ci[3])) counter <- counter + 1
        }
        
        estimates[[i]] <- counter / R
    }
    
    estimates
}

report <- function(l){
    v <- unlist(l)
    
    cat(v, "\n")
    cat("Mean =", mean(v), ", Median =", median(v), ", SD =", sd(v), ", Min =", min(v), ", Max =", max(v), "\n\n")
}



## 1 - Load and prepare data

data <- read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)

nparticipants <- (ncol(data)-1)/5
nreferents <- nrow(data)

emptyframe <- function(nrows, ncols, colprefix) {
  as.data.frame(setNames(replicate(ncols,character(nrows), simplify = F), paste(colprefix, seq(1:ncols), sep="")))
}

signs.original <- emptyframe(nreferents, nparticipants, "P")

col <- function(data, header, suffix) {data[[paste(header, suffix, sep="")]]}
for (p in 1:nparticipants) {
    gestures <- col(data, "gesture", p)
    signs.original[p] <- gestures
}


# Preparation - To reduce calculation load, I replace all sign names by an numerical ID

signs.names <- names(table(t(signs.original)))
for(i in 1:length(signs.names)){
    signs.original[signs.original == signs.names[i]] = i
}


# I will now create a dataset of N participants by sampling with replacement from the 20 participants of the actual study
N = 6000
m <- matrix(, nrow=nreferents, ncol=N)
for(i in 1:nreferents){
    m[i,] <- sample(as.numeric(as.vector(as.matrix(signs.original[i,]))), N, replace=TRUE)
}



signs.simulated <- as.data.frame(m)

kappa <- fleiss.kappa.raw.noinference(signs.simulated)

R = 1600 # number of repetitions

cat("Fleiss Kappa in Full Population: ", kappa, "\n")

###############
kappas <- sapply(1:nreferents, function(i) fleiss.kappa.for.item(signs.simulated,i))
cat("Fleiss Kappas for each referent", kappas, "\n")
##################

################### CREATE GROUPS OF LOW, MEDIUM and HIGH AGREEMENT
mask.low = match(sort(kappas), kappas)[1:14]
mask.medium = match(sort(kappas), kappas)[15:28]
mask.high = match(sort(kappas), kappas)[29:42]

signs.low <- signs.simulated[mask.low,]
signs.medium <- signs.simulated[mask.medium,]
signs.high <- signs.simulated[mask.high,]

kappa.low <- fleiss.kappa.raw.noinference(signs.low)
kappa.medium <- fleiss.kappa.raw.noinference(signs.medium)
kappa.high <- fleiss.kappa.raw.noinference(signs.high)

cat("======================= LOW, MEDIUM, HIGH =============================\n")
cat("Fleiss Kappa in LOW Population: ", kappa.low, "\n")
cat("Fleiss Kappa in MEDIUM Population: ", kappa.medium, "\n")
cat("Fleiss Kappa in HIGH Population: ", kappa.high, "\n")

kappas.low <- sapply(1:14, function(i) fleiss.kappa.for.item(signs.low,i))
kappas.medium <- sapply(1:14, function(i) fleiss.kappa.for.item(signs.medium,i))
kappas.high <- sapply(1:14, function(i) fleiss.kappa.for.item(signs.high,i))

conflevel <- .95
cat("\n\n###################################################################################################################################\n")
cat("======================= CI = .95 =============================\n")
cat("###################################################################################################################################\n\n")

###################################### OVERALL ##############################################
cat("======================= OVERALL KAPPA =============================\n")
coverage.full <- coverage.estimation(signs.simulated, kappa, R, n=20, conflev = conflevel)
cat("FULL SET:", coverage.full, "\n")
##############################################################################################################################

######################################### SPLIT INTO THREE GROUPS ######################################################
cat("======================= OVERALL LOW, MEDIUM, HIGH =============================\n")
coverage.low <- coverage.estimation(signs.low, kappa.low, R, n=20, conflev = conflevel)
coverage.medium <- coverage.estimation(signs.medium, kappa.medium, R, n=20, conflev = conflevel)
coverage.high <- coverage.estimation(signs.high, kappa.high, R, n=20, conflev = conflevel)
cat("LOW:", coverage.low, "\n")
cat("============================================================================\n")
cat("MEDIUM:", coverage.medium, "\n")
cat("============================================================================\n")
cat("HIGH:", coverage.high, "\n")
cat("============================================================================\n")
##############################################################################################################################

###################################### RANDOM WITHIN-PARTICIPANT DIFFERENCES ##################################
cat("======================= WITHIN-PARTICIPANT DIFFERENCES - FULL SET =============================\n")
coverage.diff <- coverage.estimation.referents.diff (signs.simulated, kappas, l = 100, R, n=20, conflev = conflevel)
report(coverage.diff)

#cat("======================= WITHIN-PARTICIPANT DIFFERENCES - LOW =============================\n")
#coverage.diff.low <- coverage.estimation.referents.diff (signs.low, kappas.low, l = 100, R, n=20, conflev = conflevel)
#report(coverage.diff.low)

#cat("======================= WITHIN-PARTICIPANT DIFFERENCES - MEDIUM =============================\n")
#coverage.diff.medium <- coverage.estimation.referents.diff (signs.medium, kappas.medium, l = 100, R, n=20, conflev = conflevel)
#report(coverage.diff.medium)

#cat("======================= WITHIN-PARTICIPANT DIFFERENCES - HIGH =============================\n")
#coverage.diff.high <- coverage.estimation.referents.diff (signs.high, kappas.high, l = 100, R, n=20, conflev = conflevel)
#report(coverage.diff.high)


########################################## RANDOM GROUPS OF SIZE 5 or 10 #####################################################
cat("======================= RANDOM GROUPS =============================\n")
coverage.groups.5 <- coverage.estimation.referents.groups(signs.simulated, kappas, l = 100, size = 5, R, n=20, conflev = conflevel)
coverage.groups.10 <- coverage.estimation.referents.groups(signs.simulated, kappas, l = 100, size = 10, R, n=20, conflev = conflevel)
cat("============================ SIZE = 5 ====================================\n")
report(coverage.groups.5)
cat("============================ SIZE = 10 =======================================\n")
report(coverage.groups.10)
cat("============================================================================\n")
##############################################################################################################################


######################################### INDIVIDUAL REFERENTS ######################################################
cat("======================= INDIVIDUAL REFERENTS =============================\n")
coverage.individual <- coverage.estimation.referent(signs.simulated, kappas, R, n=20, conflev = conflevel)
report(coverage.individual)


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

conflevel <- .99

cat("\n\n###################################################################################################################################\n")
cat("======================= CI = .99 =============================\n")
cat("###################################################################################################################################\n\n")

###################################### OVERALL ##############################################
cat("======================= OVERALL KAPPA =============================\n")
coverage.full.2 <- coverage.estimation(signs.simulated, kappa, R, n=20, conflev = conflevel)
cat("FULL SET:", coverage.full.2, "\n")
##############################################################################################################################


######################################### SPLIT INTO THREE GROUPS ######################################################
cat("======================= OVERALL LOW, MEDIUM, HIGH =============================\n")
coverage.low.2 <- coverage.estimation(signs.low, kappa.low, R, n=20, conflev = conflevel)
coverage.medium.2 <- coverage.estimation(signs.medium, kappa.medium, R, n=20, conflev = conflevel)
coverage.high.2 <- coverage.estimation(signs.high, kappa.high, R, n=20, conflev = conflevel)
cat("LOW:", coverage.low.2, "\n")
cat("============================================================================\n")
cat("MEDIUM:", coverage.medium.2, "\n")
cat("============================================================================\n")
cat("HIGH:", coverage.high.2, "\n")
cat("============================================================================\n")
##############################################################################################################################


###################################### RANDOM WITHIN-PARTICIPANT DIFFERENCES ##################################
cat("======================= WITHIN-PARTICIPANT DIFFERENCES - FULL SET =============================\n")
coverage.diff.2 <- coverage.estimation.referents.diff (signs.simulated, kappas, l = 100, R, n=20, conflev = conflevel)
report(coverage.diff.2)

#cat("======================= WITHIN-PARTICIPANT DIFFERENCES - LOW =============================\n")
#coverage.diff.low.2 <- coverage.estimation.referents.diff (signs.low, kappas.low, l = 100, R, n=20, conflev = conflevel)
#report(coverage.diff.low.2)

#cat("======================= WITHIN-PARTICIPANT DIFFERENCES - MEDIUM =============================\n")
#coverage.diff.medium.2 <- coverage.estimation.referents.diff (signs.medium, kappas.medium, l = 100, R, n=20, conflev = conflevel)
#report(coverage.diff.medium.2)

#cat("======================= WITHIN-PARTICIPANT DIFFERENCES - HIGH =============================\n")
#coverage.diff.high.2 <- coverage.estimation.referents.diff (signs.high, kappas.high, l = 100, R, n=20, conflev = conflevel)
#report(coverage.diff.high.2)

########################################## RANDOM GROUPS OF SIZE 5 or 10 #####################################################
cat("======================= RANDOM GROUPS =============================\n")
coverage.groups.5.2 <- coverage.estimation.referents.groups(signs.simulated, kappas, l = 100, size = 5, R, n=20, conflev = conflevel)
coverage.groups.10.2 <- coverage.estimation.referents.groups(signs.simulated, kappas, l = 100, size = 10, R, n=20, conflev = conflevel)
cat("============================ SIZE = 5 ====================================\n")
report(coverage.groups.5.2)
cat("============================ SIZE = 10 =======================================\n")
report(coverage.groups.10.2)
cat("============================================================================\n")
##############################################################################################################################


######################################### INDIVIDUAL REFERENTS ######################################################
cat("======================= INDIVIDUAL REFERENTS =============================\n")
coverage.individual.2 <- coverage.estimation.referent(signs.simulated, kappas, R, n=20, conflev = conflevel)
report(coverage.individual.2)



