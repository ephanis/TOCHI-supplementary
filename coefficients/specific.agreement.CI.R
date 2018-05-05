# author: Theophanis Tsandilas

# Code for caclulating category-specific agreement and its jackknife CIs
# The implementation of the jackknife CIs is based on the code in agreement.CI.R but has been vectorized for faster calcuating agreement for all sign categories

# REFERENCES
#
# http://www.john-uebersax.com/stat/raw.htm
# https://github.com/jmgirard/mReliability/wiki/Specific-agreement-coefficient


# ==============================================================
# This is an r function for trimming leading and trealing blanks
# ==============================================================
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) }


categories <- function(ratings){
  ratings.mat <- as.matrix(ratings) 

  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
 #   categ <- as.vector(na.omit(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
  }    

  categ
}


specific.agreement <- function(ratings, categ){
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  
  q <- length(categ)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  agree.mat <- matrix(0, nrow=n, ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else
      agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  
  sa <- agree.mat[1, ]
  
  for(k in 1:q){
    actual <- 0
    possible <- 0 
    
    for(i in 1:n){
      actual <- actual + agree.mat[i,k]*(agree.mat[i,k] - 1)
      possible <- possible + agree.mat[i,k] * (r - 1)
    }
    
    sa[k] <- actual / possible
    if(is.nan(sa[k])) sa[k] <- 0
  }
  
  names(sa) <- categ
  
  sa
}

#
# Generates the jackknife samples for rating data, assuming a fixed set of items
# and an infinite population of raters.
#
jackknife.specific.samples.random.raters <- function(ratings) {
    
    # sample size for the jackknife = number of raters
    n <- ncol(ratings)
    
    # generate jacknife samples
    jackknife_samples <- list()
    for (r in 1:n) {
        # remove all data from rater r
        ratings_subset <- ratings[,-r]
        # add the dataset to the list of samples
        jackknife_samples[[r]] <- ratings_subset
    }
    
    jackknife_samples
}


specific.agreement.CI <- function(ratings, conflevel = 0.95){
  bounds <- c(0, 1)
  
  categ <- categories(ratings)
  
  jackknife.samples <- jackknife.specific.samples.random.raters(ratings)
  jackknife.specific.CI(ratings, jackknife.samples, specific.agreement, categ, bounds, 0.95)
}


# Spitzer and Fleiss [1974] argue that specific agreement should be corrected for chance
# Here, we used the simplest chance correction described by Uebersax [1982], when chance areement specific to
# categories is common accross participants
# Frequencies of categories = chance agreement specific to categories
#
frequencies <- function(ratings, categ) {
    ratings.mat <- as.matrix(ratings)
    n <- nrow(ratings.mat) # number of subjects
    r <- ncol(ratings.mat) # number of raters
    
    q <- length(categ)
    
    counts <- as.table(rep(0, length(categ)))
    names(counts) <- categ
    observed <- table(unlist(ratings))
    
    counts[names(observed)] <- observed[names(observed)]
    
    counts / sum(counts)
}

# Chance corrected agreement specific to categories
chance.corrected.specific.agreement <- function(ratings, categ){
    sa <- specific.agreement(ratings, categ)
    pe <- frequencies(ratings, categ)
    
    (sa - pe)/(1 - pe)
}

chance.corrected.specific.agreement.CI <- function(ratings, conflevel = 0.95){
    bounds <- c(-1, 1)
    
    categ <- categories(ratings)
    
    jackknife.samples <- jackknife.specific.samples.random.raters(ratings)
    jackknife.specific.CI(ratings, jackknife.samples, chance.corrected.specific.agreement, categ, bounds, 0.95)
}


#
# Computes CIs for a statistic given a set of jackknife samples, by first computing a
# variance estimate using Equation 5.3.29 from Gwet (2014, p.153), then turning
# this variance into a CI using the t distribution.
#
# NOTE: When used with fixed raters, this function returns CIs that are slightly different and sometimes
# noticeably larger than CIs from Gwet's R scripts. His function does not seem to use his
# Equation 5.3.29 but another method, and they yield different variances. The variances produced by our function
# however agree with the variances he reports in his book on Table 5.6 p.154.
#
jackknife.specific.CI <- function(full.sample, jackknife.samples, statistic, categories, bounds, conflevel = 0.95) {
  
  # 1 - Compute the mean estimate using the full sample
  
  mean.estimate <- statistic(full.sample, categories)
  
  # 2 - Compute the Jacknife variance estimate from jackknife samples
  #     using Equation 5.3.29 from Gwet (2014, p.153).
  
  # Compute the statistic for each sample
  #jackknife.values <- sapply(jackknife.samples, statistic)
  jackknife.values <- sapply(jackknife.samples, function(x, y) statistic(x, categories))
  
  # Jacknife sample size = rater sample size
  n <- length(jackknife.values[1,])

  # Seems irrelevant to  
  # Remove NaN values (TODO: explain)
#  if (sum(is.nan(jackknife.values)) != 0 && fixNaN) {
#    jackknife.values <- ifelse(is.nan(jackknife.values), 1, jackknife.values)
#    jackknife.values <- jackknife.values[!is.na(jackknife.values)]    
#  }
  jn <- length(jackknife.values[1,])
  q <- length(categories)

  # Compute the jackknife mean estimate  
  jackknife.mean.estimate <- array(1:length(jackknife.values[,1]))
  for(k in 1:q){
    jackknife.mean.estimate[k] <- mean(jackknife.values[k,])
  }
  
  # The following assumes an infinite population, so the (1 - gr) term disappears
  variance <-  array(1:length(jackknife.values[,1]))
  for(k in 1:q){
    variance[k] <- (jn - 1) / jn * sum((jackknife.values[k,] - jackknife.mean.estimate[k]) ^ 2)  
  }
  
  # 1 - Turn the mean and variance into a CI using the t distribution.
  #     This is not explained in Gwet's book but it is the method he uses in his
  #     R scripts. Also see Equation 7 from Abdi and Williams (2010, p.3).
  stderr <- sqrt(variance)
  
  # lower confidence bound
  lower.bound <- mean.estimate - stderr * qt(1 - (1 - conflevel) / 2, n - 1)
  
  # upper confidence bound
  upper.bound <- mean.estimate + stderr * qt(1 - (1 - conflevel) / 2, n - 1)
  
  # Limit to range (Gwet does limits upper bound to 1 in his R scripts)
  lower.bound <- cap(lower.bound, bounds[1], bounds[2])
  upper.bound <- cap(upper.bound, bounds[1], bounds[2])  
  
  # Return the mean and confidence limits
  # Create a 3-row matrix to store: (1) the mean estimate, (2) the lower bound, (3) the upper bound 
  #cis <- matrix(ncol = q, nrow = 3)
  #cis[1,] <- mean.estimate
  #cis[2,] <- lower.bound
  #cis[3,] <- upper.bound
  
  tab <- table(unlist(full.sample))
  target <- c(names(tab))
  #print(tab)
  
  cis.df <- data.frame(
    Category = categories, 
    Specific = c(mean.estimate),
    Lower = c(lower.bound),
    Upper = c(upper.bound),
    Count = 0:0,
    Freq = 0:0
    )

  cis.df <- cis.df[match(target,cis.df$Category),]
  cis.df$Count <- as.numeric(tab)
  cis.df$Freq <- as.numeric(100 * tab / sum(tab))
  
  cis.df
}

#
# Cap value
#
cap <- function(x, xmin, xmax) {
  n <- length(x)
  for(i in 1:n) x[i] <- min(xmax, max(xmin, x[i]))
  
  x
}

toStringFormat <- function(ci){
    sprintf("%2.1f%%, %.2f, 95%% CI [%.2f, %.2f]", ci[6], ci[2], ci[3], ci[4])
} 

