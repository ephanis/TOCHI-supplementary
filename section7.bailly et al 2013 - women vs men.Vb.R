
##################################################################
# Analysis of Bailly et al's dataset to test Vatavu & Wobbrock's (2016) conclusions
# about differences between women and men.
#
# Using Vatavu & Wobbrock's Vb significance test
#
# author: Theophanis Tsandilas
##################################################################

rm(list=ls()) # Clean up R's memory

source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")
source("section6.between.vatavu2016.Vb.R")


vb.compare <- function(group1, group2){
    nref <- nrow(group1)
    N1 <- ncol(group1)
    N2 <- ncol(group2)
    
    counter <- 0
    for(r in 1:nref){
        signs1 <- group1[r,]
        signs2 <- group2[r,]
        
        ar1 <- percent.agreement.raw.noinference(signs1)
        ar2 <- percent.agreement.raw.noinference(signs2)
    
        pval <- Vb.pvalue(ar1, ar2, N1, N2)
        
        if(pval <.05){
            cat("x")
            counter <- counter + 1
        }
        else cat("-")
    }
    
    cat("\n> ", counter, " differences found to be statistically significant\n")
        
    counter
}

vb.compare.random <- function(all, N1){
    N <- ncol(all)
    
    mask <- sample(1:N, N1)
    group1 <- all[mask]
    group2 <- all[-mask]
    
    vb.compare(group1, group2)
}


## 1 - Load and prepare data

data <- read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)

# This mask identifies the ID of male participants
males_mask <- c(1, 5, 6, 8, 9, 13, 15, 18, 20)

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

signs_males <- signs_gestures[males_mask]
signs_females <- signs_gestures[-males_mask]

percent <- jack.CI.random.raters(signs_males, percent.agreement)
printCI("AR for males: ", percent)

percent <- jack.CI.random.raters(signs_females, percent.agreement)
printCI("AR for females: ", percent)

count <- vb.compare(signs_males, signs_females)
cat("\nNumber of statistically significant differences for Males/Females (alpha = .05) = ", count, "\n")

repetitions <- 30 # Increase for a better but slower estimate ( = 1000 in the evaluation presented in the article)
cat("\nRejections of null hypothesis for ", repetitions, " random partitions of 9 and 11 participants (this can take long...)\n", sep="")

result <- replicate(repetitions, vb.compare.random(signs_gestures, 9))

#cat("Results: ", result, "\n")
cat("\nMean number of statistically significant differences (alpha = .05) in random partitions: ", mean(result), "\n", sep="")

