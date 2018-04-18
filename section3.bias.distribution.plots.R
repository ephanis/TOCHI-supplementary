
#===========================================================================
#===========================================================================
# Script used to plot theoretical bias distribution functions
#
# Theophanis Tsandilas
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

# Create a population of size N from a normal frequency distribution with mean = 1 and std. dev = sd
# P is the size of the original population - any large enough value will suffice
normal <- function(N, sd, P = 10000){
    sample <- sample(c(1:P), N, replace = TRUE, prob=dnorm(mean=1,sd=sd,c(1:P)))
    
    sample
}


zipf.plot <- function(beta, nsamples, N, ymax = .8, color = rgb(0.3, 0.3, 0.8)){
    sample <- zipf(nsamples, beta)
    
    freq <- table(sample)
    total <- sum(freq)
    
    freq <- head(freq, n = N) / total
    
    plot(freq, xlim = c(1, N), ylim=c(0, ymax), col = color, cex=1,  lwd=1, xaxt='n', yaxt='n', ann=FALSE, bty="n",  type = "o")
    #lines(freq, col = color)
    
    freq
}


normal.plot <- function(sd, nsamples, N, ymax = .8, color = rgb(0.3, 0.3, 0.8)){
    sample <- normal(nsamples, sd)
    
    freq <- table(sample)
    total <- sum(freq)
    
    freq <- head(freq, n = N) / total
    
    plot(freq, xlim = c(1, N), ylim=c(0, ymax), col = color, cex=1,  lwd=1, xaxt='n', yaxt='n', ann=FALSE, bty="n", type = "o")
    #lines(freq, col = color)
    
    freq
}


#Plot the distributions
par(mfrow=c(2,3))

B = .646
zipf.plot(B, 200000, 15, ymax = 0.4)
axis(2, at = pretty(seq(0, 0.4, by = .1)), lab=pretty(seq(0, 0.4, by = .1))*100, las=1)
#axis(1, at = seq(1,15, by=1))


B = .306
zipf.plot(B, 200000, 15, ymax = 0.4)
axis(2, at = pretty(seq(0, 0.4, by = .1)), lab=pretty(seq(0, 0.4, by = .1))*100, las=1)

B = .151
zipf.plot(B, 200000, 15, ymax = 0.4)
axis(2, at = pretty(seq(0, 0.4, by = .1)), lab=pretty(seq(0, 0.4, by = .1))*100, las=1)


sd = 2.58
normal.plot(sd, 20000000, 15, ymax = 0.4)
axis(2, at = pretty(seq(0, 0.4, by = .1)), lab=pretty(seq(0, 0.4, by = .1))*100, las=1)

sd = 5.42
normal.plot(sd, 200000, 15, ymax = 0.4)
axis(2, at = pretty(seq(0, 0.4, by = .1)), lab=pretty(seq(0, 0.4, by = .1))*100, las=1)

sd = 11.02
normal.plot(sd, 200000, 15, ymax = 0.4)
axis(2, at = pretty(seq(0, 0.4, by = .1)), lab=pretty(seq(0, 0.4, by = .1))*100, las=1)



