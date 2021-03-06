

#===========================================================================================
# Simuation experiment for evaluating the Type I error of the jackknife technique
# Theophanis Tsandilas
#======================================================================================

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

typeI.estimation <- function(N, sds, R = 100, alpha = .05){
    L <- length(sds)
    errors <- rep(0, L)
    
    conflev = 1 - alpha
    
    for(r in 1:R){
        cat("\r", r, " out of ", R, " : ", errors/r)
        
        # 1st random sample of size N (1st referent)
        samples1 <- mapply(population.create, rep(N, L), sds)
        
        # 2nd random sample of size N (2nd referent)
        samples2 <- mapply(population.create, rep(N, L), sds)
        
        for(i in 1:L){
            ci <- jack.CI.diff.random.raters(as.data.frame(t(samples1[, i])),
                as.data.frame(t(samples2[, i])), percent.agreement, confint = conflev)
                
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



# Create a population of size N from a normal frequency distribution with mean = 1 and std. dev = sd
# P is the size of the population - any large enough value will suffice
population.create <- function(N, sd, P = 10000){
    sample <- sample(c(1:P), N, replace = TRUE, prob=dnorm(mean=1,sd=sd,c(1:P)))
    
    sample
}


########################################################################
########################################################################
###### Running this whole script can take long...

# These are parameters for the normal distribution.
# They have been empirically approximated to produce distributions that correspond to different AR levels (AR = .1, .2, ..., .9)
sds <- c(5.42, 2.56, 1.63, 1.15, .88, .69, .576, .493, .416)

R <- 1600 # Number of trials for estimating the Type I error

#################################################################
#################################################################
#################################################################
alpha <- .05
errors <- typeI.estimation(20, sds, R, alpha)

cat("\n=========== Results for n = 20 and alpha = .05 (mean and 95% CIs) for AR = .1, .2, ..., .9  =================\n")
print(errors)
cat("=================================================\n\n")

#################################################################
#################################################################
#################################################################

