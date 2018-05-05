
##################################################################
# Analysis of Bailly et al's dataset (from original raw data)
#
# authors: Pierre Dragicevic, Theophanis Tsandilas
##################################################################

# This code uses Wobbrock's terminology:
#   participant = rater
#   referent    = subject (item)
#   sign        = category
#
# There are 20 participants and 42 referents

rm(list=ls()) # Clean up R's memory

source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")

## 1 - Load and prepare data

data = read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)

# For each participant, there are five columns, where the first captures the key and the second captures the gesture
signs_keys <- data[, seq(2, ncol(data), by=5)] # These are participants' proposals of keys
signs_gestures <- data[, seq(3, ncol(data), by=5)] # These are participants' proposals of key gestures

# Replace the column names by the participant IDs
names(signs_keys) <- paste0("P", 1:ncol(signs_keys))
names(signs_gestures) <- paste0("P", 1:ncol(signs_gestures))

## 2 - Compute agreements with CIs

analyze <- function(ratings) {
  printRatingsInfo(ratings)

  fleiss.pe <- fleiss.chance.agreement.raw.noinference(ratings)
  cat("Chance Agreement: ", fleiss.pe, "\n")
  

  percent.with.replacement <- jack.CI.random.raters(ratings, percent.agreement.with.replacement)
  printCI("Percent agreement with replacement", percent.with.replacement)
  
  percent <- jack.CI.random.raters(ratings, percent.agreement)
  printCI("                 Percent agreement", percent)

  Kappa <- jack.CI.random.raters(ratings, fleiss.kappa)
  printCI("                     Fleiss' Kappa", Kappa)

  alpha <- jack.CI.random.raters(ratings, krippen.alpha)
  printCI("               Krippendorf's alpha", alpha)  
  
  coef <- jack.CI.random.raters(ratings, bp.coeff)
  printCI("    Brennan-Prediger's coefficient", coef)  
  
}


## 2 bis - Compute agreement differences with CIs

analyze.diff <- function(all.ratings, ratings1, ratings2) {
    printRatingsInfo(ratings1)
    printRatingsInfo(ratings2)
    
    diff.percent.with.replacement <- jack.CI.diff.random.raters(ratings1, ratings2, percent.agreement.with.replacement)
    printCI("  Difference in percent agr. with replacement", diff.percent.with.replacement)
    
    diff.percent <- jack.CI.diff.random.raters(ratings1, ratings2, percent.agreement)
    printCI("  Difference in percent agreement", diff.percent)
    
    # This assumes that chance agreement is the same for the two sets and equal to the overall chance agreement
    fleiss.pe <- fleiss.chance.agreement.raw.noinference(all.ratings)
    fleiss.pooled <- function(ratings) {generic.kappa(ratings, fleiss.pe)}
    diff.Kappa <- jack.CI.diff.random.raters(ratings1, ratings2, fleiss.pooled)
    printCI("  Difference in Fleiss' Kappa (pooled pe)", diff.Kappa)
    
    # This is a different approach (different pe for each set of referents) -- Normally, it is not what we want.
    diff.fleiss <- jack.CI.diff.random.raters(ratings1, ratings2, fleiss.kappa)
    printCI("  Difference in Fleiss' Kappa (different pe for each set)", diff.fleiss)
}


cat("\n## Keys only ## \n\n")
analyze(signs_keys)

cat("\n## Gestures only ## \n\n")
analyze(signs_gestures)


## 3 - Estimate the difference in agreement between directional and other gestures
cat("\n## Directional minus non-directional gestures ##\n")

directional <- c("top", "bottom", "left", "right", "previous", "next", "Previous", "Next")
allcmds <- data[,1]
matched.rows <- grep(paste(directional,collapse="|"), allcmds)
directional.ratings <- signs_gestures[matched.rows,]
other.ratings <- signs_gestures[-matched.rows,]
cat("Selected", nrow(directional.ratings), "and", nrow(other.ratings), "rows.\n")

analyze.diff(signs_gestures, directional.ratings, other.ratings)



## 4 - Calculate and plot individual agreement coefficients (and their CIs) for individual gesture signs

refs <- list()
k <- list()
l <- list()
u <- list()

cat("\n## Fleiss for Individual Gesture Signs ## \n\n")
for(index in 1:nrow(signs_gestures)){
    ci <- jack.CI.random.raters.fleiss.kappa.for.item(signs_gestures, index)
    
    print( sprintf("%s, Kappa = %.3f, 95%% CI = [%.3f, %.3f]", data[index, 1], ci[1], ci[2], ci[3]))
    
    refs[index] <- data[index, 1]
    k[index] <- ci[1]
    l[index] <- ci[2]
    u[index] <- ci[3]
}

fleiss.df <- data.frame(Referent = as.character(refs), Mean = as.double(k), Low = as.double(l), Upper = as.double(u))
fleiss.df <- fleiss.df[with(fleiss.df, order(-Mean)),]

referents <- fleiss.df$Referent

library(ggplot2)

plot <- ggplot(fleiss.df, aes(x = Referent, y = Mean)) +
geom_point(size = 1.8, color = "#0000aa", shape = 1) +
scale_x_discrete(name ="Referent", limits = referents) +
scale_y_continuous(breaks=c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1))+
theme_bw() +
geom_errorbar(aes(ymax = Upper, ymin = Low), width=0.3, color = "#0000aa") +
theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 53, vjust = 1, hjust = 1, size = 12), axis.text.y = element_text(size = 12)) +
theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#000000", size=0.03))

print(plot)

