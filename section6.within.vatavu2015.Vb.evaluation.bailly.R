#===========================================================================================
# Evaluation of the Type I Error of the Vb test (Vatavu & Wobbrock, 2016)
# By generating populations from Bailly et al's dataset
# Author: Theophanis Tsandilas
#======================================================================================

#rm(list=ls()) # Clean up R's memory

source("section6.between.vatavu2016.Vb.R")
source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")

################################################################################################
################################################################################################
# Code for simulating the creation of pairs of samples
# (comparison between agreement rates for two only referents) with specific AR levels
################################################################################################
################################################################################################


# n is the number of participants that I randomly draw from the overall population. R = number or repetitions
typeI.estimation <- function(signs.full, R = 100, n = 20, alpha = .05){
    L <- length(signs.full[,1])
    errors <- rep(0, L)
    
    counter <- 0
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        samples1 <- signs.full[,sample(1:ncol(signs.full), size = n)]
        samples2 <- signs.full[,sample(1:ncol(signs.full), size = n)]

        for(i in 1:L){
            AR1 <- ar(t(samples1[i,]), n)
            AR2 <- ar(t(samples2[i,]), n)
            
            p <- Vb.pvalue(AR1, AR2, n, n)
            if(p < alpha) errors[i] <- errors[i] + 1
        }

        flush.console()
    }
    
    cat("\n")
    
    estimates <- list()
    
    for(i in 1:L){
        res <- binom.test(errors[i], R)
        estimates[[i]] <- c(res$estimate, res$conf.int[1], res$conf.int[2])
    }
    
    estimates
}


ar <- function(samples, N){
    t <- table(unlist(samples))
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
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


# I will now create a population of 1000 participants by sampling with replacement from the 20 participants of the actual study
N = 6000
m <- matrix(, nrow=nreferents, ncol=N)
for(i in 1:nreferents){
    m[i,] <- sample(as.numeric(as.vector(as.matrix(signs.original[i,]))), N, replace=TRUE)
}

signs.simulated <- as.data.frame(m)

agr <- percent.agreement.raw.noinference(signs.simulated)
cat("AR in Full Population: ", agr, "\n")
for(i in 1:length(signs.simulated[,1]))
cat("AR for referent ", i, " = ",  percent.agreement.raw.noinference(signs.simulated[i,]), "\n")

errors.05 <- typeI.estimation(signs.simulated, R=1600, alpha = .05)
print(errors.05)

errors.01 <- typeI.estimation(signs.simulated, R=1600, alpha = .01)
print(errors.01)

