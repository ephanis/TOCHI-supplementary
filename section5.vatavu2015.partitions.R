

#============================================================================================================
# Producing the partitions for the probabilistic reasoning of Vatavu & Wobbrock 2015 (magnitude of agreement)
# Theophanis Tsandilas
#============================================================================================================

rm(list=ls()) # Clean up R's memory

#Requires the installation of Partitions package: https://cran.r-project.org/web/packages/partitions/partitions.pdf
library("partitions")


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

# Calculate the mean AR for a given number of N participants, where agreements and frequencies vectors have already been computed
agreement.mean <- function(N, agreements, frequencies){
    pairs.max <- N * (N - 1) /2
    
    sum(agreements*frequencies / sum(frequencies) / pairs.max)
}


# Calculate the mean AR for a given number of N participants when partition frequencies are considered
AR.mean <- function(N){
    partitions <- as.matrix(parts(N))
    agreements <- partitions.agreements(partitions)
    
    freq <- partitions.frequencies(partitions)
    
    mean <- agreement.mean(N, agreements, freq)
    
    cat("Partition frequencies: ", freq, "\n")
    print(sprintf("N = %2.0f, AR = %.3f", N, mean))
    
}


# Calculate the mean AR for a given number of N participants when partition frequencies are not considered
AR.simple.mean <- function(N){
    partitions <- as.matrix(parts(N))
    agreements <- partitions.agreements(partitions)
    
    freq <- rep(1, length(partitions))
    
    mean <- agreement.mean(N, agreements, freq)
    
    print(sprintf("N = %2.0f, AR = %.3f", N, mean))
}



# Plots the distribution of AR values accroding to the solution (incorrect) by Vatavu & Wobbrock 2015
# n: number of participants, h: bin size
AR.distribution.plot <- function(N, h){
    total <- N*(N-1)/2
    
    partitions <- as.matrix(parts(N))
    
    agr.samples <- partitions.agreements(partitions) / total
    
    n <- 1/h
    agr.samples <- round(agr.samples / h)
    samples <- table(seq(from = 0, to = 1.0, by = h)) * 0
    freq <- table(unlist(agr.samples))
    names(freq) <- as.numeric(names(freq)) * h
    samples[names(freq)] <- freq[names(freq)]
    samples <- t(samples / sum(samples))
    
    points <- barplot(samples, col = rgb(.8,.8,.95), ylim=range(c(0, .3)), axisnames = FALSE, axes = FALSE)
    
    ####
    size <- length(points)
    filtered <- points[seq(1, size, size/10)]
    filtered <- append(filtered, points[size])
    axis(1, at = filtered, labels = seq(0, 1, .1))
    filtered <- filtered + (filtered[2] - filtered[1]) /2
    axis(1, at = filtered, labels = NA, tck = -0.02)
    
    axis(2, at=seq(0, .3, 0.1))
    
    samples
}



# Examples
cat("Matrix with all possible partitions for 6 participants:\n")
print(as.matrix(parts(6)))

cat("\n", "Mean AR for 6 participants if considering that all partitions appear with the same frequency\n")
AR.simple.mean(6)

cat("\n\nMean AR for 6 participants if taking into account the frequency of partitions\n")
AR.mean(6)

# Plot the distribution of AR values (Vatavu & Wobbrock 2015) for 20 participants and bin size h = .05
AR.distribution.plot(20, .05)
