source("section6.within.vatavu2015.Vrd.R")

############################################################################
# Testing the implementaiton of te Vrd test on the dataset of Bailly et. al.
#
# Author: Theophanis Tsandilas
############################################################################


data = read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)
nparticipants <- (ncol(data)-1)/5
nreferents <- nrow(data)

emptyframe <- function(nrows, ncols, colprefix) {
    as.data.frame(setNames(replicate(ncols,character(nrows), simplify = F), paste(colprefix, seq(1:ncols), sep="")))
}

gest.signs <- emptyframe(nreferents, nparticipants, "P")

col <- function(data, header, suffix) {data[[paste(header, suffix, sep="")]]}
for (p in 1:nparticipants) {
    gestures <- col(data, "gesture", p)
    gest.signs[p] <- gestures
}

mat <- t(as.matrix(gest.signs))

cat("==================================================================================", "\n")
cat("Using Vrd test to evaluate the overall effect of referent type on agreement rate: ", "\n")
cat("==================================================================================", "\n")
pairs <- agreementPairs(mat)
cochran.qtest <- vrd(pairs)
print(cochran.qtest)

#cat("=======================================================================================", "\n")
#cat("Using alternative use of Cochran's Q test to account for dependencies:", "\n")
#cat("=======================================================================================", "\n")
#pairs <- agreementPairs.independent(mat)
#cochran.qtest <- vrd(pairs)
#print(cochran.qtest)
