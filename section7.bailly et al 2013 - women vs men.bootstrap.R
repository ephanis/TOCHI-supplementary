
##################################################################
# Analysis of Bailly et al's dataset to test Vatavu & Wobbrock's conclusions
# about differences between women and men.
#
# Using the bootstrap method for independent groups
#
# author: Theophanis Tsandilas
##################################################################

rm(list=ls()) # Clean up R's memory

source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")
source("section6.between.bootstrap.R") # Implementation of Bootstrap methods for Agreement Rates

crossesZero <- function(ci){
    if(ci[2] <= 0 && ci[3] >=0) TRUE
    else FALSE
}

bootstrap.compare <- function(group1, group2){
    nref <- nrow(group1)
    N1 <- ncol(group1)
    N2 <- ncol(group2)
    
    counter <- 0
    for(r in 1:nref){
        signs1 <- group1[r,]
        signs2 <- group2[r,]
        
        ci <- bootstrap.diff.ci(signs1, signs2, 3000)
        # Using 3000 bootstrap iterations here, but we need more for a better accuracy

        if(!crossesZero(ci)){
            cat("x")
            counter <- counter + 1
        }
        else cat("-")
     
    }
    
    cat("\n> ", counter, " differences found to be statistically significant\n")
    
    counter
}


bootstrap.compare.random <- function(all, N1){
    N <- ncol(all)
    
    mask <- sample(1:N, N1)
    group1 <- all[mask]
    group2 <- all[-mask]
    
    bootstrap.compare(group1, group2)
}


## 1 - Load and prepare data

data <- read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)

# This mask identifies the ID of male participants
males_mask <- c(1, 5, 6, 8, 9, 13, 15, 18, 20)

signs_gestures <- data[, seq(3, ncol(data), by=5)] # These are participants' proposals of key gestures

signs_males <- signs_gestures[males_mask]
signs_females <- signs_gestures[-males_mask]

percent <- jack.CI.random.raters(signs_males, percent.agreement)
printCI("AR for Males: ", percent)

percent <- jack.CI.random.raters(signs_females, percent.agreement)
printCI("AR for Females: ", percent)

cat("\nThe following computations can take overnight!!!\n")
count <- bootstrap.compare(signs_males, signs_females)
cat("Number of statistically significant differences for Males/Females (alpha = .05) = ", count, "\n")

repetitions <- 30 # Increase for a better but slower estimate ( = 1000 in the evaluation presented in the article)
cat("\nRejections of null hypothesis for ", repetitions, " random partitions of 9 and 11 participants (this can take long...)\n", sep="")

result <- replicate(repetitions, bootstrap.compare.random(signs_gestures, 9))

#cat("Results: ", result, "\n")
cat("\nMean number of statistically significant differences (alpha = .05) in random partitions: ", mean(result), "\n", sep="")



