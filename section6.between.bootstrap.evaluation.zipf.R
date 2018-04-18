

#===============================================================================================================
# Simulation experiment for evaluating the bootstrap method for comparing agreement rates of indpendent groups
# Implementation for two only groups
# Author: Theophanis Tsandilas
#===============================================================================================================

rm(list=ls()) # Clean up R's memory

#Requires the installation of Partitions package: https://cran.r-project.org/web/packages/partitions/partitions.pdf
library("partitions")

source("section6.between.bootstrap.R")

################################################
################################################
################################################
# Code for simulating the creation of populations with specific AR levels

######################
######################


library("Rmisc")
library("zipfR") # Check http://www.r-bloggers.com/the-zipf-and-zipf-mandelbrot-distributions/


crossesZero <- function(ci){
    if(ci[2] <= 0 && ci[3] >=0) TRUE
    else FALSE
}


typeI.estimation <- function(N, N1, N2, betas, R = 100, plevel = .05, B = 1000){
    L <- length(betas)
    errors <- rep(0, L)
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # Find populations of size N
        populations <- mapply(population.create, rep(N, L), betas)
        
        for(i in 1:L){
            population <- populations[, i]
            
            mask <- sample(length(population), N1)
            sample1 <- population[mask]
            sample2 <- sample(population[-mask], N2) # so that sample2 is independent of sample1
            
            ci <- bootstrap.diff.ci(sample1, sample2, B, plevel)
            #ci <- bootstrap.diff.ci.boot(sample1, sample2, B, plevel) # to test the alternative implementation (resuts should be the same)
            
            if(!crossesZero(ci)) errors[i] <- errors[i] + 1
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


# Create a population from a Zipf distribution
population.create <- function(N, beta){
    ZM <- lnre("zm", alpha = 1/2, B = beta)
    zmsample <- rlnre(ZM, n = N)
    
    #print(ar(zmsample, N))
    
    zmsample
}


ar <- function(samples, N){
    t <- table(unlist(samples))
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
}


########################################################################
########################################################################
###### Running this whole script can take very long...

# These are parameters for the Zip-Mandelbrot distribution.
# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
betas <- c(.306, .64, 1.05, 1.56, 2.25, 3.26, 5.0, 8.3, 18.4)

R <- 1600 # Number of trials for estimating a Type I error
B <- 3000 # bootstrap resamplig iterations (B = 10000 in our evaluations)

#################################################################
#################################################################
#################################################################
alpha <- .05
boot.errors <- typeI.estimation(100, 20, 20, betas, R, alpha, B)

cat("\n=========== Results for n = 20 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(boot.errors)
cat("=================================================\n\n")

#################################################################
#################################################################
#################################################################

