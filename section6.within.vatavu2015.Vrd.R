
#===========================================================================================
# Statistical test (Vrd) for within-participant AR comparisons as proposed by Vatavu & Wobbrock 2015
# Author: Theophanis Tsandilas
#======================================================================================

#Requires the installation of the coin package, which provides an implementation of Cohran's Q test:
# https://cran.r-project.org/web/packages/coin/coin.pdf
# -> install.packages("coin")


rm(list=ls()) # Clean up R's memory
library("coin")

# data is a matrix containing the signs proposed for each referent (columns) by each participant (rows)
# return a binary matrix with agreement values (1 = agreement, 0 = disagreement) for all the pairs of participants (rows) for each referent (columns)
agreementPairs <- function(data){
    nrows <- nrow(data)
    ncols <- ncol(data)
    
    mat <- matrix(, nrow = nrows * (nrows - 1) / 2, ncol = ncols)
    index <- 1
    
    for(i1 in 1:(nrows - 1)){
        for(i2 in (i1 + 1):nrows){
            for(r in 1:ncols){
                mat[index, r] <- as.integer(data[i1, r] == data[i2, r])
            }
            
            index <- index + 1
        }
    }
    
    mat
}


# coin implementation : https://cran.r-project.org/web/packages/coin/coin.pdf
vrd <- function(pairs){
    nrefs <- ncol(pairs)
    npairs <- nrow(pairs)
    
    #vectorize the pairs
    response <- c(t(pairs))
    fact <- gl(nrefs, 1, nrefs*npairs)
    block <- gl(npairs, nrefs)
    
    symmetry_test(response~fact|block,teststat = "quad")
}


############## Correct generation of pairs for Cohran's Q test by only considering the n - 1 independent pairs
######################################################################################
agreementPairs.independent <- function(data){
    nrows <- nrow(data)
    ncols <- ncol(data)
    
    mat <- matrix(, nrow = nrows - 1, ncol = ncols)
    
    for(i in 1:(nrows - 1)){
        for(r in 1:ncols){
            mat[i, r] <- as.integer(data[i, r] == data[i + 1, r])
        }
    }
    
    mat
}
##############################################################################################

######################################################################################
######################################################################################
# Vatavu's implemetation for 2 referents (k = 2)
# Note that their equation 9 (k > 2) seems to be incorrect
# Vatavu et al' equations are unnecessary. There are numerous statistical packages that calculate Cohran's Q test (see above)

# agreement
ar <- function(v){
    sum(v)/length(v)
}

# co-agreement
cr <- function(v1, v2){
    sum(v1*v2)/length(v1)
}


vrd.vatavu <- function(pairs){
    n <- nrow(pairs)
    
    ar.vector <- apply(pairs, 2, ar)
    cr.total <- cr(pairs[,1], pairs[,2])
    
    n*(ar.vector[1] - ar.vector[2])^2 / (ar.vector[1] + ar.vector[2] - 2*cr.total)
}
######################################################################################
######################################################################################





