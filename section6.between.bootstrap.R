

#===========================================================================================
# Statistical tests for between-participant AR comparisons with bootstrapping techniques
# Implementation for two only groups
# Theophanis Tsandilas
#======================================================================================

library(boot)


# Basic calculation of AR measure
boot.ar <- function(samples, N){
    t <- table(unlist(samples))
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
}

boot.ar.simple <- function(samples){
    N <- length(samples)
    
    t <- table(unlist(samples))
    agreements <- t*(t-1)
    ar <- sum(agreements) / (N * (N-1))
    
    ar
}

# Basic calculation of difference with AR measures
boot.ar.diff <- function(samples1, N1, samples2, N2){
    diff <- boot.ar(samples1, N1) - boot.ar(samples2, N2)
    
    diff
}

boot.ar.diff.simple <- function(samples1, samples2){
    N1 <- length(samples1)
    N2 <- length(samples2)
    
    diff <- boot.ar(samples1, N1) - boot.ar(samples2, N2)
    
    diff
}

bootstrap.diff <- function(group1, N1, group2, N2, R){
    diffs <- replicate(R, boot.ar.diff(sample(group1, replace = TRUE), N1, sample(group2, replace = TRUE), N2))
    
    diffs
}

bootstrap.diff.ci <- function(group1, group2, R = 1000, plevel = .05){
    N1 <- length(group1)
    N2 <- length(group2)
    
    diff <- sort(bootstrap.diff(group1, N1, group2, N2, R))
  
   # uncomment this to check the distributions of AR differences that we get...
   #hist(diff)
  
    perc <- R*plevel/2
    
    low <- floor(perc)
    upper <- ceiling(R - perc)
    
    c(boot.ar.diff(group1, N1, group2, N2), diff[low], diff[upper])
}



# Alternative implementation by using R's boot package to calulcate confidence intervals
# The percentile method seems to work best. All the other methods yield large Type I errors.
# It's probably because difference distributions are symmetric
bootstrap.diff.ci.boot <- function(group1, group2, R = 1000, plevel = .05){
     stat <- function(data, data2, indices) {
            d1 <- data[indices]
         
            d2 <- sample(data2, replace = TRUE)
         
            return(boot.ar.diff.simple(d1, d2))
         }
     
     bootobj <- boot(data = group1, statistic = stat, R = R, data2 = group2)
     ci <- tryCatch(boot.ci(bootobj, type='perc', conf = 1 - plevel), error=function(...) list(fail=TRUE))
     
     if(length(ci$fail) && ci$fail) {
         c(0, 0, 0)
     } else {
         c(boot.ar.diff.simple(group1, group2), ci$perc[4], ci$perc[5])
         #c(ci$t0, ci$bca[4], ci$bca[5])
     }
}







