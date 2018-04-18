#===========================================================================================
# Evaluation of the Type I Error of the Vrd test (Vatavu & Wobbrock, 2015)
# Zipfian distributions
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

# It requires the zipfR package for the Zip-Mandelbrot distribution
library("zipfR") # Check http://www.r-bloggers.com/the-zipf-and-zipf-mandelbrot-distributions/


typeI.estimation <- function(N, betas, R = 100, alpha = .05){
    L <- length(betas)
    errors <- rep(0, L)
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # 1st random sample of size N (1st referent)
        samples1 <- mapply(population.create, rep(N, L), betas)
        
        # 2nd random sample of size N (2nd referent)
        samples2 <- mapply(population.create, rep(N, L), betas)
        
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



# Create a population from a Zipf distribution
population.create <- function(N, beta){
    ZM <- lnre("zm", alpha = 1/2, B = beta)
    zmsample <- rlnre(ZM, n = N)
    
    zmsample
}

########################################################################
########################################################################
###### Running this whole script can take long... You can reduce the number of trials for faster (but less accurate results)

# These are parameters for the Zip-Mandelbrot distribution.
# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
betas <- c(.306, .64, 1.05, 1.56, 2.25, 3.26, 5.0, 8.3, 18.4)

R <- 1600 # Number of trials for estimating a Type I error

#################################################################
#################################################################
alpha <- .05
errors <- typeI.estimation(20, betas, R, alpha)

cat("\n=========== Results for n = 20 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(errors)
cat("=================================================\n\n")

