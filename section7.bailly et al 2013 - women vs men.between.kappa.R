
##################################################################
# Analysis of Bailly et al's dataset to test agreement between women and men
#
# Using the bootstrap method for independent groups
#
# author: Theophanis Tsandilas
##################################################################

rm(list=ls()) # Clean up R's memory

source("coefficients/agreement.CI.R")
source("coefficients/agreement.coefficients.R")

data <- read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)

# This mask identifies the ID of male participants
males_mask <- c(1, 5, 6, 8, 9, 13, 15, 18, 20)

gestures <- data[, seq(3, ncol(data), by=5)] # These are participants' proposals of key gestures

# This mask identifies the ID of male participants
men_mask <- c(1, 5, 6, 8, 9, 13, 15, 18, 20)
gestures_men <- gestures[men_mask]
gestures_women <- gestures[-men_mask]

kappa.men <- jack.CI.random.raters(gestures_men, fleiss.kappa)
printCI("Fleiss' Kappa for Men: ", kappa.men)

kappa.women <- jack.CI.random.raters(gestures_women, fleiss.kappa)
printCI("Fleiss' Kappa for Women: ", kappa.women)

kappa.delta <- fleiss.kappa.bootstrap.diff.ci(1:nrow(gestures), gestures_women, gestures_men, R = 3000)
# R: number of bootstrap samples
printCI("Difference in Fleiss' Kappa between Women and Men: ", kappa.delta)
