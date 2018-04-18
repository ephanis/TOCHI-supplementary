#===========================================================================================
# Evaluation of the Type I Error of the Vrd test (Vatavu & Wobbrock, 2015)
# Half-Normal Distributions
# Author: Theophanis Tsandilas
#======================================================================================

rm(list=ls()) # Clean up R's memory

source("section6.within.vatavu2015.Vrd.R")

################################################################################################
################################################################################################
# Code for simulating the creation of pairs of samples
# (comparison between agreement rates for two only referents) with specific AR levels
################################################################################################
################################################################################################

typeI.estimation <- function(N, sds, R = 100, alpha = .05){
    L <- length(sds)
    errors <- rep(0, L)
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # 1st random sample of size N (1st referent)
        samples1 <- mapply(population.create, rep(N, L), sds)
        
        # 2nd random sample of size N (2nd referent)
        samples2 <- mapply(population.create, rep(N, L), sds)
        
        for(i in 1:L){
            mat <- cbind(samples1[, i], samples2[, i])
             pairs <- agreementPairs(mat)
             
             # Uncomment the following line (and comment the previous one) to apply Cochran's Q test on n - 1 independent pairs of participants as exaplained in the article
             #pairs <- agreementPairs.independent(mat)

            tryCatch({
                test <- vrd(pairs)
                if(pvalue(test) < alpha) errors[i] <- errors[i] + 1
            }, error = function(e) 1)
            
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
# P is the size of the population - any large enough value will suffice
population.create <- function(N, sd, P = 10000){
    sample <- sample(c(1:P), N, replace = TRUE, prob=dnorm(mean=1,sd=sd,c(1:P)))
    
    sample
}


ar <- function(samples){
	N <- length(samples)
	
    t <- table(unlist(samples))
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
}


########################################################################
########################################################################
###### Running this whole script can take long... You can reduce the number of trials for faster (but less accurate results)

# These are parameters for the discrete half normal distribution.
# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
sds <- c(5.42, 2.56, 1.63, 1.15, .88, .69, .576, .493, .416)

R <- 1600 # Number of trials for estimating a Type I error

#################################################################
#################################################################
alpha <- .05
errors <- typeI.estimation(20, sds, R, alpha)

cat("\n=========== Results for n = 20 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(errors)
cat("=================================================\n\n")

