

#===========================================================================================
# Statistical tests for between-participant Kappa comparisons with bootstrapping techniques
# Implementation for two only groups
# Theophanis Tsandilas
#======================================================================================
source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")

################## Assuming here that each group can have its own chance agreement

bootstrap.diff <- function(r, group1, group2, R){
    diffs <- replicate(R, fleiss.kappa.for.item(sample(group1, replace = TRUE), r) -  fleiss.kappa.for.item(sample(group2, replace = TRUE), r))
    
    diffs
}


# r is the indexes of the referents of interest
bootstrap.diff.ci <- function(r, group1, group2, R = 1000, plevel = .05){
    diff <- sort(bootstrap.diff(r, group1, group2, R))
    
    # uncomment this to check the distributions of Kappa differences that we get...
    #hist(diff)
    
    perc <- R*plevel/2
    
    low <- floor(perc)
    upper <- ceiling(R - perc)
    
    c(diff[low], diff[upper])
}


