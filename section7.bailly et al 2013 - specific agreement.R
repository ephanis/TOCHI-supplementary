
###################################################################################
# Analysis of Bailly et al's dataset using specific agreement (Only for gesture signs)
#
# author: Theophanis Tsandilas
###################################################################################

rm(list=ls()) # Clean up R's memory
source("coefficients/specific.agreement.CI.R")

## 1 - Load and prepare data

data = read.csv("data/bailly et al 2013 - dataset.csv", stringsAsFactors=F)
nparticipants <- (ncol(data)-1)/5
nreferents <- nrow(data)

emptyframe <- function(nrows, ncols, colprefix) {
  as.data.frame(setNames(replicate(ncols,character(nrows), simplify = F), paste(colprefix, seq(1:ncols), sep="")))
}
signs <- emptyframe(nreferents, nparticipants, "P")

col <- function(data, header, suffix) {data[[paste(header, suffix, sep="")]]}
for (p in 1:nparticipants) {
  gestures <- col(data, "gesture", p)
  signs[p] <- gestures
}


# 2 - Calculate Specific Agreement for gesture signs

cat("##### ", "Sign Frequency and Agreement Specific to Individual Signs" , " #####\n\n")


ci.df <- specific.agreement.CI(signs)

# To correct for chance agreement specific to signs, use instead the following line
# ci.df <- chance.corrected.specific.agreement.CI(signs)

ci.df <- ci.df[with(ci.df, order(-Freq)),]

for(q in 1:nrow(ci.df)){
    cat(toStringFormat(ci.df[q,]))
    cat("\t")
    cat(as.character(ci.df[q,1]))
    cat("\n")
}



# 3 - Plot Specific Agreement

###########################################################
###########################################################
categories <- ci.df$Category

library(ggplot2)

x_angle <- 40

plot <- ggplot(ci.df, aes(x = Category, y = Specific)) +
geom_point(size = 1.8, color = "#0000aa", shape = 1) +
# ylim(-0.3, 1.0) +
scale_x_discrete(name ="Sign", limits = categories) +
scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1))+
theme_bw() +
geom_errorbar(aes(ymax = Upper, ymin = Lower), width=0.1, color = "#0000aa") +
theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = x_angle, vjust = 1, hjust = 1, size = 10), axis.text.y = element_text(size = 12)) +
theme(panel.border = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "#000000", size=0.03))

print(plot)

