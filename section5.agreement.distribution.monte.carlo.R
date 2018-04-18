
#===========================================================================
#===========================================================================
# Script used to derive the null distribution for AR (or Fleiss' Kappa) values.
# We assume that bias distributions are Zipfian
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


monte.carlo.simulation.ar <- function(nparticipants, nreferents, rep = 100, B = .65){
    ar <- rep(0, rep)
    
    for(r in 1:rep){
        data <- emptyframe(nreferents, nparticipants, "P")
        
        for(i in 1:nparticipants){
            data[i] <- zipf(nreferents, B)
        }
        
        ar[r] <- percent.agreement(data)
    }
    
    ar
}

monte.carlo.simulation.kappa <- function(nparticipants, nreferents, rep = 100, B = .65){
    kappa <- rep(0, rep)
    
    for(r in 1:rep){
        data <- emptyframe(nreferents, nparticipants, "P")
        
        for(i in 1:nparticipants){
            data[i] <- zipf(nreferents, B)
        }
        kappa[r] <- fleiss.kappa.raw.noinference(data)
    }
    
    kappa
}


plot.distribution <- function(agr.samples, N = 40, min = 0, max = 1.0, axis = TRUE){
    binsize = 1 / N
    
    samples <- table(seq(from = min, to = max, by = binsize)) * 0
    
    agr.samples <- round(N * agr.samples)/N
    freq <- table(unlist(agr.samples))
    # names(freq) <- as.numeric(names(freq)) / N
    
    print(freq)
    
    samples[names(freq)] <- freq[names(freq)]
    
    samples <- t(samples / sum(samples))
    
    points <- barplot(samples, col = rgb(.8,.8,.95), ylim=range(c(0, 1)), axisnames = FALSE, axes = FALSE)
    
    #####
    size <- length(points)
    filtered <- points[seq(1, size, size/(10*(max-min)))]

    filtered <- append(filtered, points[size])

    axis(1, at = filtered, labels = seq(min, max, .1))
    
    filtered <- filtered + (filtered[2] - filtered[1]) /2
    axis(1, at = filtered, labels = NA, tck = -0.02)
    
    if(axis){
        axis(2, at=pretty(seq(0, 1, 0.5)), lab=pretty(seq(.1, .9, by = .1))*100, las=1)
        axis(2, at=seq(.1, .9, 0.1), labels = NA, tck = -0.025)
    }
    #####
    
    samples
}


rep = 100 #For better results, please increase the number of Monte Carlo iterations
N = 80
min = 0
max = .5

par(mfrow=c(2,3), mar = c(2,3,2,0.1))

B = .646
#B = .151
res <- monte.carlo.simulation.ar(20, 40, rep, B)
#Uncomment the following for Kappa
#res <- monte.carlo.simulation.kappa(20, 40, rep, B)
plot.distribution(res, N, min, max)

B = .306
res <- monte.carlo.simulation.ar(20, 40, rep, B)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)

B = .151
res <- monte.carlo.simulation.ar(20, 40, rep, B)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)


B = .646
#B = .151
res <- monte.carlo.simulation.ar(10, 40, rep, B)
plot.distribution(res, N, min, max)

B = .306
res <- monte.carlo.simulation.ar(10, 40, rep, B)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)

B = .151
res <- monte.carlo.simulation.ar(10, 40, rep, B)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)
