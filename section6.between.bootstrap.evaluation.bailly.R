
##################################################################
# Analysis of Bailly et al's dataset to test Vatavu & Wobbrock's conclusions
# about differences between women and men.
#
# Using the bootstrap method for independent groups
#
# author: Theophanis Tsandilas
##################################################################

rm(list=ls()) # Clean up R's memory

source("section6.between.bootstrap.kappa.R")


crossesZero <- function(ci){
    if(ci[1] <= 0 && ci[2] >=0) TRUE
    else FALSE
}

bootstrap.compare <- function(all, n1, R.Boot = 1000, conflev = .95){
    n <- ncol(all)
    
    mask <- sample(1:n, n1)
    group1 <- all[mask]
    group2 <- all[-mask]
    
    nref <- nrow(group1)
    
    counter <- 0
    for(r in 1:nref){
        signs1 <- group1[r,]
        signs2 <- group2[r,]
        
        ci <- bootstrap.diff.ci(r, group1, group2, R = R.Boot, 1 - conflev)
        
        # Using 3000 bootstrap iterations here, but we need more for a low significance level

        if(!crossesZero(ci)){
            cat("x")
            counter <- counter + 1
        }
        else cat("-")
    }
    
    cat("\n> ", counter, " differences found to be statistically significant\n")
    
    counter
}



## 1 - Load and prepare data

data <- read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)

nparticipants <- (ncol(data)-1)/5
nreferents <- nrow(data)

emptyframe <- function(nrows, ncols, colprefix) {
  as.data.frame(setNames(replicate(ncols,character(nrows), simplify = F), paste(colprefix, seq(1:ncols), sep="")))
}

signs_gestures <- emptyframe(nreferents, nparticipants, "P")


col <- function(data, header, suffix) {data[[paste(header, suffix, sep="")]]}
for (p in 1:nparticipants) {
    gestures <- col(data, "gesture", p)
    signs_gestures[p] <- gestures
}


repetitions <- 100 # Increase for a better but slower estimate ( = 1000 in the evaluation presented in the article)
cat("\nRejections of null hypothesis for ", repetitions, " random partitions of 10 and 10 participants (this can take long...)\n", sep="")

result <- replicate(repetitions, bootstrap.compare(signs_gestures, 10))

#cat("Results: ", result, "\n")
ci <- CI(result)
cat("\nMean number of statistically significant differences (alpha = .05) in random partitions: ", ci[2], ", 95% CI = [", ci[3], ", ", ci[1], "]\n", sep="")



