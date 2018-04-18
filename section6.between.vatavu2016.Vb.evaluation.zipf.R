
#===========================================================================================
# Evaluation of the Vb statistical test by Vatavu & Wobbrock 2015
# Author: Theophanis Tsandilas
#======================================================================================

rm(list=ls()) # Clean up R's memory

source("section6.between.vatavu2016.Vb.R")


# It requires the zipfR package for the Zip-Mandelbrot distribution
library("zipfR") # Check http://www.r-bloggers.com/the-zipf-and-zipf-mandelbrot-distributions/


typeI.estimation <- function(N, N1, N2, betas, R = 100, plevel = .05){
    L <- length(betas)
    errors <- rep(0, L)
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # Find populations of size N
        populations <- mapply(population.create, rep(N, L), betas)
        #samples1 <- mapply(population.create, rep(N1, L), betas)
        #samples2 <- mapply(population.create, rep(N2, L), betas)
        
        
        for(i in 1:L){
            #sample1 <- samples1[,i]
            #sample2 <- samples2[,i]
            population <- populations[, i]
            
            mask <- sample(length(population), N1)
            sample1 <- population[mask]
            sample2 <- sample(population[-mask], N2) # so that sample2 is independent of sample1
            
            AR1 <- ar(sample1, N1)
            AR2 <- ar(sample2, N2)
            
            p <- Vb.pvalue(AR1, AR2, N1, N2)
            if(p < plevel) errors[i] <- errors[i] + 1
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
    
    zmsample
}


zipf <- function(N, beta){
    ZM <- lnre("zm", alpha = 1/2, B = beta)
    
    zmsample <- rlnre(ZM, n = N)
    
    ar(zmsample, N)
}


ar <- function(samples, N){
    t <- table(unlist(samples))
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
}


########################################################################
########################################################################
###### Running this whole script can take long...

# These are parameters for the Zip-Mandelbrot distribution.
# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
betas <- c(.306, .64, 1.05, 1.56, 2.25, 3.26, 5.0, 8.3, 18.4)

R <- 1600 # Number of trials for estimating a Type I error

#################################################################
#################################################################
#################################################################
alpha <- .05
vb.errors <- typeI.estimation(100, 20, 20, betas, R, alpha)

cat("\n=========== Results for n = 20 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(vb.errors)
cat("=================================================\n\n")

