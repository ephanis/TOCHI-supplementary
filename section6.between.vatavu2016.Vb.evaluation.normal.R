
#===========================================================================================
# Evaluation of the Vb statistical test by Vatavu & Wobbrock 2015
# Author: Theophanis Tsandilas
#======================================================================================

rm(list=ls()) # Clean up R's memory

source("section6.between.vatavu2016.Vb.R")


typeI.estimation <- function(N, N1, N2, sds, R = 100, plevel = .05){
    L <- length(sds)
    errors <- rep(0, L)
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # Find populations of size N
        populations <- mapply(population.create, rep(N, L), sds)
        #samples1 <- mapply(population.create, rep(N1, L), sds)
        #samples2 <- mapply(population.create, rep(N2, L), sds)
        
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


# Create a population of size N from a normal frequency distribution with mean = 1 and std. dev = sd
# P is the size of the normal population - any large enough value will suffice
population.create <- function(N, sd, P = 1000){
    sample <- sample(c(1:P), N, replace = TRUE, prob=dnorm(mean=1,sd=sd,c(1:P)))
    
    sample
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

# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
sds <- c(5.42, 2.56, 1.63, 1.15, .88, .69, .576, .493, .416)

R <- 1600 # Number of trials for estimating a Type I error

#################################################################
#################################################################
#################################################################
alpha <- .05
vb.errors <- typeI.estimation(100, 20, 20, sds, R, alpha)

cat("\n=========== Results for n = 20 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(vb.errors)
cat("=================================================\n\n")

