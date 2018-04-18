

#===========================================================================================
# Statistical tests for between-participant AR comparisons as proposed by Vatavu & Wobbrock 2016
# Implementation for two only groups
# Author: Theophanis Tsandilas
#======================================================================================

# Requires the installation of Partitions package: https://cran.r-project.org/web/packages/partitions/partitions.pdf
library("partitions")

Vb.pvalue <- function(ar1, ar2, N1, N2){
    n1 <- N1 * (N1 - 1) /2
    n2 <- N2 * (N2 - 1) /2
    
    Pr(ar1 * n1, ar2 * n2, N1, N2)
}


# Calculation of probability of observing a1, a2 or more extreme proporions for groups with N1 and N2 participants
Pr <- function(a1, a2, N1, N2){
    n1 <- N1 * (N1 - 1) /2
    n2 <- N2 * (N2 - 1) /2
    
    pr1 <- probabilities(N1)
    pr2 <- probabilities(N2)
    
    v_ <- Vb(a1, a2, n1, n2)
    
    p12 <- 0.5 * pr1[as.character(a1)]*pr2[as.character(a2)]
    if(is.na(p12)) p12 = 0
    
    for(e1 in as.integer(names(pr1))){
        for(e2 in as.integer(names(pr2))){
            v <- Vb(e1, e2, n1, n2)
            
            if(v > v_){
                p <- pr1[as.character(e1)]*pr2[as.character(e2)]
                p12 <- p12 + p
            }

        }
    }
  
    p12
}


#Calucation of Vb for groups with pairs numbers n1 and n2 and agreement configurations e1 and e2
Vb <- function(e1, e2, n1, n2){
    
    e = e1 + e2
    a = e / (n1 + n2)
    expected1 = n1*a
    expected2 = n2*a
    
    (e1 - expected1)^2 + (e2 - expected2)^2
}


# Calculate the probability of each possible number of agreeing pairs for N participants
probabilities <- function(N){
    partitions <- as.matrix(parts(N))
    agreements <- partitions.agreements(partitions)
    
    freq <- partitions.frequencies(partitions)
    
    frame <- data.frame(agreements=agreements, f=freq)
    
    with(frame, tapply(f, agreements, FUN = sum)) / sum(freq)
}


# Calculate the number of agreeing pairs for all partitions
partitions.agreements <- function(partitions){
    m <- partitions*(partitions - 1) / 2
    colSums(m)
}


# Calculate the frequency of observing partitions
partitions.frequencies <- function(partitions){
    npartitions <- length(partitions[1,])
    nparticipants <- length(partitions[,1])
    
    f <- rep(1, npartitions)
    
    #print(partitions)
    
    for(i in 1:npartitions){
        partition <- partitions[,i]
    
        n <- nparticipants
        for(j in 1:nparticipants){
            f[i] <- f[i]*choose(n, partition[j])
            
            n <- n - partition[j]
        }
        
        t <- table(partition)
        
        t <- t[names(t) != '0'] # Do not conider repititions of 0
        t <- t[t > 1]
       
        if(length(t) == 0) denom <- 1
        else denom <- prod(factorial(t))
       
        f[i] <- f[i] / denom
    }
    
    f
}


# Calculate the mean AR for a given number of N participants
AR.mean <- function(N){
    partitions <- as.matrix(parts(N))
    agreements <- partitions.agreements(partitions)
    
    # print(agreements)
    
    freq <- partitions.frequencies(partitions)
    
    # print(freq)
    
    mean <- agreement.mean(N, agreements, freq)
    
    print(sprintf("N = %2.0f, AR = %.3f", N, mean))
    
}


# Calculate the mean AR for a given number of N participants, where agreements and frequencies vectors have already been computed
agreement.mean <- function(N, agreements, frequencies){
    pairs.max <- N * (N - 1) /2
    
    sum(agreements*frequencies / sum(frequencies) / pairs.max)
}

