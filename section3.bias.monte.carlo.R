
#===========================================================================
#===========================================================================
# Monter Calro simulations showing the problem of chance agreement in the presence of bias for an infinite number of signs
#
# Theophanis Tsandilas, Inria 2017
#===========================================================================
#===========================================================================

rm(list=ls()) # Clean up R's memory

library("zipfR") # Check http://www.r-bloggers.com/the-zipf-and-zipf-mandelbrot-distributions/
source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")

emptyframe <- function(nrows, ncols, colprefix) {
    as.data.frame(setNames(replicate(ncols,character(nrows), simplify = F), paste(colprefix, seq(1:ncols), sep="")))
}

# Create a population from a Zipf distribution
zipf <- function(N, beta){
    ZM <- lnre("zm", alpha = 1/2, B = beta)
    zmsample <- rlnre(ZM, n = N)
    
    zmsample
}


# Create a population of size N from a half-normal frequency distribution with mean = 1 and std. dev = sd
# P is the size of the original population - any large enough value will suffice
normal <- function(N, sd, P = 10000){
    sample <- sample(c(1:P), N, replace = TRUE, prob=dnorm(mean=1,sd=sd,c(1:P)))
    
    sample
}


#Someone can take the linear combination of the two distributions to generate new distribution functions
mix <- function(N, beta, sd, a.zipf, a.normal){
    sample.zipf <- as.numeric(as.character(zipf(N, beta)))
    sample.normal <- as.numeric(as.character(normal(N, sd)))
    
    sample(c(sample.zipf, sample.normal), size=N, prob=c(rep(a.zipf,N), rep(a.normal,N)))
}


monte.carlo.simulation.zipf <- function(nparticipants, nreferents, rep = 100, B = .65){
    ar <- rep(0, rep)
    chance <- rep(0, rep)
    kappa <- rep(0, rep)
    alpha <- rep(0, rep)

    signs <- rep(0, rep)
    
    for(r in 1:rep){
        data <- emptyframe(nreferents, nparticipants, "P")
        
        for(i in 1:nparticipants){
            data[i] <- zipf(nreferents, B)
        }
        
        ar[r] <- percent.agreement(data)
        chance[r] <- fleiss.chance.agreement.raw.noinference(data)
        kappa[r] <- fleiss.kappa.raw.noinference(data)
        alpha[r] <- krippen.alpha.raw.noinference(data)
        
        signs[r] <- length(table(unlist(data)))
    }
    
    list(ar, chance, kappa, alpha, signs)
}


monte.carlo.simulation.normal <- function(nparticipants, nreferents, rep = 100, sd = 2){
    ar <- rep(0, rep)
    chance <- rep(0, rep)
    kappa <- rep(0, rep)
    alpha <- rep(0, rep)
    
    signs <- rep(0, rep)
    
    for(r in 1:rep){
        data <- emptyframe(nreferents, nparticipants, "P")
        
        for(i in 1:nparticipants){
            data[i] <- normal(nreferents, sd)
        }
        
        ar[r] <- percent.agreement(data)
        chance[r] <- fleiss.chance.agreement.raw.noinference(data)
        kappa[r] <- fleiss.kappa.raw.noinference(data)
        alpha[r] <- krippen.alpha.raw.noinference(data)
        
        signs[r] <- length(table(unlist(data)))
    }
    
    list(ar, chance, kappa, alpha, signs)
}


monte.carlo.simulation.mix <- function(nparticipants, nreferents, rep = 100, B = .65, sd = 2, a.zipf = .5, a.normal = .5){
    ar <- rep(0, rep)
    chance <- rep(0, rep)
    kappa <- rep(0, rep)
    alpha <- rep(0, rep)
    
    signs <- rep(0, rep)
    
    for(r in 1:rep){
        data <- emptyframe(nreferents, nparticipants, "P")
        
        for(i in 1:nparticipants){
            data[i] <- mix(nreferents, B, sd, a.zipf, a.normal)
        }
        
        ar[r] <- percent.agreement(data)
        chance[r] <- fleiss.chance.agreement.raw.noinference(data)
        kappa[r] <- fleiss.kappa.raw.noinference(data)
        alpha[r] <- krippen.alpha.raw.noinference(data)
        
        signs[r] <- length(table(unlist(data)))
    }
    
    list(ar, chance, kappa, alpha, signs)
}

getResult <- function(res, index, statistic){
    statistic(unlist(res[index]))
}

# Examples
res <- monte.carlo.simulation.zipf(20, 40, rep, B)
cat("################ Zipf Distribution #######################", "\n")
cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")

res <- monte.carlo.simulation.normal(20, 40, rep, sd)
cat("################ Normal Distribution #######################", "\n")
cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")


B = .306
sd = 5.42
# Examples
res <- monte.carlo.simulation.zipf(20, 40, rep, B)
cat("################ Zipf Distribution #######################", "\n")
cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")

res <- monte.carlo.simulation.normal(20, 40, rep, sd)
cat("################ Normal Distribution #######################", "\n")
cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")


B = .151
sd = 11.02
# Examples
res <- monte.carlo.simulation.zipf(20, 40, rep, B)
cat("################ Zipf Distribution #######################", "\n")
cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")

res <- monte.carlo.simulation.normal(20, 40, rep, sd)
cat("################ Normal Distribution #######################", "\n")
cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")


# This is for the same experiment but taking the linear combination of zipf and normal
#B = .646
#sd = 2.58
#res <- monte.carlo.simulation.mix(20, 40, rep, B, sd)
#cat("################ Mixed Distribution #######################", "\n")
#cat("Mean number of observed signs =", getResult(res, 5, mean), "\n")
#cat("Percent Agreement (AR) =", getResult(res, 1, mean), "\n")
#cat("Fleiss' Chance Agreement =", getResult(res, 2, mean), "\n")
#cat("Fleiss' Kappa =", getResult(res, 3, mean), "\n")
#cat("Krippendorf's alpha =", getResult(res, 4, mean), "\n")

