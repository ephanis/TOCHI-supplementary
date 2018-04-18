
#===========================================================================
#===========================================================================
# Script used to derive the distributions of agreement differences in Figure 6.
#
# Theophanis Tsandilas, Inria 2017
#===========================================================================
#===========================================================================

rm(list=ls()) # Clean up R's memory

library("zipfR") # Check http://www.r-bloggers.com/the-zipf-and-zipf-mandelbrot-distributions/

ar <- function(samples){
    N <- length(samples)
    t <- table(samples)
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
}

# Create a population from a Zipf distribution
zipf <- function(N, beta){
    ZM <- lnre("zm", alpha = 1/2, B = beta)
    zmsample <- rlnre(ZM, n = N)
    
    zmsample
}

normal <- function(N, sd, P = 10000){
  sample <- sample(c(1:P), N, replace = TRUE, prob=dnorm(mean=1,sd=sd,c(1:P)))
  
  sample
}

monte.carlo.simulation.ar <- function(nparticipants, rep = 100, distr, param){
    ars <- rep(0, rep)
    
    for(r in 1:rep){
        data1 <- distr(nparticipants, param)
        data2 <- distr(nparticipants, param)
        
        ars[r] <- ar(data1) - ar(data2)
    }
    
    ars
}


plot.distribution <- function(agr.samples, N = 40, min = 0, max = 1.0, axis = TRUE){
    binsize = 1 / N
    
    samples <- table(seq(from = min, to = max, by = binsize)) * 0
    
    agr.samples <- round(N * agr.samples)/N
    freq <- table(unlist(agr.samples))
    freq <- freq/sum(freq)

    samples[names(samples)] <- freq[names(samples)]
    samples <- t(samples)
    
    points <- barplot(samples, col = rgb(.8,.8,.95), ylim=range(c(0, .8)), axisnames = FALSE, axes = FALSE)
    
    #####
    size <- length(points)
    filtered <- points[seq(1, size, size/(10*(max-min)))]

    filtered <- append(filtered, points[size])

    axis(1, at = filtered, labels = seq(min, max, .1))
    
    filtered <- filtered + (filtered[2] - filtered[1]) /2
 #   axis(1, at = filtered, labels = NA, tck = -0.02)
    
    if(axis){
        axis(2, at=pretty(seq(0, 1, 0.5)), lab=pretty(seq(.1, .9, by = .1))*100, las=1)
        #axis(2, at=seq(.1, .9, 0.1), labels = NA, tck = -0.025)
    }
    #####
    
    samples
}


# Example - We assume here that q = 10 (pe = .125), n = 10, and h = 1/20 = .05 (bin size)
#res1 <- monte.carlo(10, 10, 20)

rep = 5000
N = 20
min = -.5
max = .5

par(mfrow=c(2,3), mar = c(2,3,2,0.1))

B = .306
res <- monte.carlo.simulation.ar(20, rep, zipf, B)
plot.distribution(res, N, min, max)

B = 2.25
res <- monte.carlo.simulation.ar(20, rep, zipf, B)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)

B = 18.4
res <- monte.carlo.simulation.ar(20, rep, zipf, B)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)


sd = 5.42
res <- monte.carlo.simulation.ar(40, rep, normal, sd)
plot.distribution(res, N, min, max)

sd = .88
res <- monte.carlo.simulation.ar(40, rep, normal, sd)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)

sd = .416
res <- monte.carlo.simulation.ar(40, rep, normal, sd)
plot.distribution(res, N, min, max, axis = FALSE)
#plot.distribution(res, N, min, max)
