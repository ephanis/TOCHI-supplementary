

#===========================================================================================================================
# Simuation experiment for evaluating the Type I error of the bootstrap method for within-participants hypothesis testing
# Theophanis Tsandilas
#===========================================================================================================================

rm(list=ls()) # Clean up R's memory

source("coefficients/agreement.coefficients.R")
source("coefficients/agreement.CI.R")

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

typeI.estimation <- function(N, betas, R = 100, alpha = .05){
    L <- length(betas)
    errors <- rep(0, L)
    
    conflev = 1 - alpha
    
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # 1st random sample of size N (1st referent)
        samples1 <- mapply(population.create, rep(N, L), betas)
        
        # 2nd random sample of size N (2nd referent)
        samples2 <- mapply(population.create, rep(N, L), betas)
        
        for(i in 1:L){
            ci <- boot.CI.diff.random.raters(as.data.frame(t(samples1[, i])),
            as.data.frame(t(samples2[, i])), percent.agreement, confint = conflev, R=3000) #You can increase the number of iterations for higher accuracy (R=10000 in our experiments)
            
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
    
    zmsample
}


########################################################################
########################################################################
###### Running this whole script can take long...

# These are parameters for the Zip-Mandelbrot distribution.
# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
betas <- c(.306, .64, 1.05, 1.56, 2.25, 3.26, 5.0, 8.3, 18.4)

R <- 1600 # Number of trials for estimating the Type I error

#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
#################################################################
alpha <- .05
errors <- typeI.estimation(10, betas, R, alpha)

cat("\n=========== Results for n = 10 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(errors)
cat("=================================================\n\n")

#################################################################
#################################################################
#################################################################


