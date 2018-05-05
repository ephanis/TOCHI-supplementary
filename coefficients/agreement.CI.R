# Authors: Pierre Dragicevic and Theophanis Tsandilas


# REFERENCES
#
# Gwet (2004)
# Abdi, H. & Williams, L.J. (2010). Jackknife. In N.J. Salkind, D.M., Dougherty, & B. Frey (Eds.): Encyclopedia of Research Design. Thousand Oaks (CA): Sage. pp. 655-660.
# Mike Meredith (2010) Bias and confidence limits for indices [online] http://www.wcsmalaysia.org/analysis/Biod_indices_bstrap.htm
# Mark Gardener () Community Ecology: Analytical Methods Using R and Excel
# Dawson (2011) Measurement of work group diversity. Aston University PhD thesis.

# This assumes that coefficients are already loaded, e.g.:
# source("agreement.coefficients.R")

library(boot) # Library for boostrap CIs

#
# Computes a 95% jackknife CI for a given statistic on a rating dataset, assuming a fixed set of
# items and random drawing from an infinite population of raters. Uses Gwet's
# method of "inference conditionally upon the item sample" (Gwet, 2014, p.152-154).
jack.CI.random.raters <- function(ratings, statistic, confint = 0.95, bounds = "auto", fixNaN = FALSE, debug= FALSE) {
  
  if (bounds == "auto") {bounds <- statistic.bounds(statistic)}
  
  jackknife.samples <- jackknife.samples.random.raters(ratings)
  jackknife.CI(ratings, jackknife.samples, statistic, bounds, confint, fixNaN, debug)
}



#
# Uses pooled pe from the full population (Fleiss)
#
jack.CI.diff.random.raters.fleiss.kappa <- function(all.ratings, item.indexes.1, item.indexes.2, conflevel = 0.95) {
    
    statistic <- function(ratings) {
        # Compute Fleiss' chance agreement
        pe <- fleiss.chance.agreement.raw.noinference(ratings)
        
        # Compute Fleiss kappa for the first group of items
        pa1 <- percent.agreement(ratings[item.indexes.1,])
        
        # Compute Fleiss kappa for the second group of items
        pa2 <- percent.agreement(ratings[item.indexes.2,])
        
        (pa1 - pa2) / (1 - pe)
    }
    
    jackknife.samples <- jackknife.samples.random.raters(all.ratings)
    jackknife.CI(all.ratings, jackknife.samples, statistic, c(-2, 2), conflevel, fixNaN=F, debug=F)
}



# An alternative implementation for computing the CI on the difference. For kappa coefficients, it is not completely accurate, as it does not calculate pe over the full set of referents
# The previous function should be preferered!!!
# Computes a CI on the difference between two datasets with the exact same raters: statistic(ratings1) - statistic(ratings2).
# ratings1 and ratings2 should include the same raters, in the same order.
#
jack.CI.diff.random.raters <- function(ratings1, ratings2, statistic, confint = 0.95, bounds = "auto", fixNaN = FALSE, debug= FALSE) {
  
  if (bounds == "auto") {
    bounds <- statistic.bounds(statistic)
    bounds <- c(bounds[1] - bounds[2], bounds[2] - bounds[1])
  }
    
  # Resamples should be in the same order in the two rating datasets. This wouldn't work with bootstrapping!
  jackknife.samples1 <- jackknife.samples.random.raters(ratings1)
  jackknife.samples2 <- jackknife.samples.random.raters(ratings2)
  
  jackknife.CI.diff(ratings1, ratings2, jackknife.samples1, jackknife.samples2, statistic, bounds, confint, fixNaN, debug)
}



# Computes a 95% BCa bootstrap CI for a given statistic on a rating dataset, assuming a fixed set of
# items and random drawing from an infinite population of raters.
#
# The results are often not very satisfactory. The jackknife seems to result in better Type I errors
#
boot.CI.random.raters <- function(ratings, statistic, confint = 0.95, bounds = "auto", R = 3000) {

  if (bounds == "auto") {bounds <- statistic.bounds(statistic)}
  
  tratings <- t(ratings)
  
  stat <- function(data, indices) {
    d <- data[indices,]
    return(statistic(t(d)))
  }
  
  bootobject <- boot(data = tratings, statistic = stat, R=R)
  #The percentile method seems to produce better results than other methods (e.g., the bca method).
  
  b <- tryCatch(boot.ci(bootobject, type='perc'),
                    error=function(...) list(fail=TRUE))
                    
  if(length(b$fail) && b$fail) {
      ci <- c(NA,NA)
      warning('could not obtain bootstrap confidence interval')
    } else {
      ci <- c(b$perc[4], b$perc[5])
  }

  # Limit to theoretical bounds
  ci[1] <- cap(ci[1], bounds[1], bounds[2])
  ci[2] <- cap(ci[2], bounds[1], bounds[2])  

  c(statistic(ratings), ci[1], ci[2])
}


##########
##########
# Uses pooled q
#
jack.CI.random.raters.bp.coeff.for.item <- function(all.ratings, item.index, conflevel = 0.95) {

  statistic <- function(ratings) {
    
    # Compute total number of categories
    q <- length(get.categories(ratings))
    
    # Extract ratings for specific item
    item.ratings <- ratings[item.index,]
        
    # Compute bp coef for the item
    bp.coeff.q(item.ratings, q)
  }
  
  jackknife.samples <- jackknife.samples.random.raters(all.ratings)
  jackknife.CI(all.ratings, jackknife.samples, statistic, c(-1, 1), conflevel, fixNaN=F, debug=F)
}




#
# Uses pooled pe from Fleiss
#
jack.CI.random.raters.fleiss.kappa.for.item <- function(all.ratings, item.index, conflevel = 0.95) {
  
  statistic <- function(ratings) {
    
    # Compute Fleiss' chance agreement
    pe <- fleiss.chance.agreement.raw.noinference(ratings)
    
    # Compute Fleiss kappa for the specific item
    pa <- percent.agreement(ratings[item.index,])
    
    (pa - pe) / (1 - pe)
  }
  
  jackknife.samples <- jackknife.samples.random.raters(all.ratings)
  jackknife.CI(all.ratings, jackknife.samples, statistic, c(-1, 1), conflevel, fixNaN=F, debug=F)
}


#
# Bootstrap Percent Agreement for Item
boot.CI.random.raters.fleiss.kappa.for.item <- function(all.ratings, item.index, R = 1000) {
    
    tratings <- t(all.ratings)
    
    stat <- function(data, pos, indices) {
        d <- data[indices,]

        # Compute Fleiss' chance agreement
        pe <- fleiss.chance.agreement.raw.noinference(t(d))
        #pe <- 0
        
        # Compute Fleiss kappa for the specific item
        pa <- percent.agreement(t(t(d)[pos,]))
        
        return((pa - pe) / (1 - pe))
    }
    
    bootobject <- boot(data = tratings, statistic = stat, R = R, pos = item.index)
    
    b <- tryCatch(boot.ci(bootobject, type='perc'),
    error=function(...) list(fail=TRUE))
    
    if(length(b$fail) && b$fail) {
        ci <- c(NA,NA)
        warning('could not obtain the bootstrap confidence interval')
    } else {
        ci <- c(b$perc[4], b$perc[5])
    }
    
    # Limit to theoretical bounds
    ci[1] <- cap(ci[1], -1, 1)
    ci[2] <- cap(ci[2], -1, 1)
    
    c(stat(tratings, item.index), ci[1], ci[2])
}



##########################################################################################
#
# Returns a statistic's theoretical bounds
#
statistic.bounds <- function(statistic) {
  if (identical(statistic, percent.agreement)) {
    return (c(0, 1))
  } else if (identical(statistic, percent.agreement.with.replacement)) {
    return (c(0, 1))  
  } else {
    return (c(-1, 1)) # assume it's a corrected agreement score
  }
  
  # return c(-inf, inf)
}

#
# Computes the percent agreement
#
percent.agreement <- function(ratings) {
  percent.agreement.raw.noinference(ratings)
}

#
# Compute the percent agreement with replacement, or Hill's N2 index, or Wobbrock's initial agreement rate.
#
percent.agreement.with.replacement <- function(ratings) {
  
  # Compute this from regular percent agreement using the linear relationship as explained by
  # Dawson (2011) p.70, equation at the middle of the page. Gwet optimized this already so we're saving time.
  
  R <- ncol(ratings) # number of raters
  pa <- percent.agreement.raw.noinference(ratings)
  # Note: pa = (R*par - 1) / (R - 1)
  par <- (pa * (R - 1) + 1) / R
  par
}

#
# Extracts the agreement coefficient from
# Gwet's implementation of Fleiss Kappa function
#
fleiss.kappa <- function(ratings) {
  fleiss.kappa.raw.noinference(ratings)
}

#
# Extracts the agreement coefficient from
# Gwet's implementation of Fleiss Kappa function
#
chance.agreement <- function(ratings) {
    fleiss.chance.agreement.raw.noinference(ratings)
}

#
# Extracts the agreement coefficient from
# Gwet's implementation of Krippendorf's alpha function
#
krippen.alpha <- function(ratings) {
  
  # This method throws an error if there's only one item:
  #    Error in weights.mat * (pi.vec %*% t(pi.vec)) : non-conformable arrays 
  # So try to avoid this and return NA instead
  if (nrow(ratings) == 1)
    return (NA)
  else
    krippen.alpha.raw.noinference(ratings)
}


#
# Extracts the agreement coefficient from
# Gwet's implementation of Brennan-Prediger coefficient function
#
bp.coeff <- function(ratings) {
  bp.coeff.raw.noinference(ratings)
}


#
# Computes the Brennan-Prediger coefficient with the number of categories q provided as an argument.
#
# If q is equal to the number of categories in the ratings, returns the same result as bp.coeff.
# Used for computing item-specific agreements.
#
# TODO: replace with generic kappa
#
bp.coeff.q <- function(ratings, q) {
  pa <- percent.agreement.raw.noinference(ratings)
  (pa - 1/q) / (1 - 1/q)
}

generic.kappa <- function(ratings, pe) {
  pa <- percent.agreement.raw.noinference(ratings)
  (pa - pe) / (1 - pe)  
}

###################################################################################################

#
# Generates the jackknife samples for rating data, assuming a fixed set of items
# and an infinite population of raters.
#
jackknife.samples.random.raters <- function(ratings) {

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


# Computes CIs for a statistic given a set of jackknife samples, by first computing a
# variance estimate using Equation 5.3.29 from Gwet (2014, p.153), then turning
# this variance into a CI using the t distribution.
#
# NOTE: When used with fixed raters, this function returns CIs that are slightly different and sometimes
# noticeably larger than CIs from Gwet's R scripts. His function does not seem to use his
# Equation 5.3.29 but another method, and they yield different variances. The variances produced by our function
# however agree with the variances he reports in his book on Table 5.6 p.154.
#
jackknife.CI <- function(full.sample, jackknife.samples, statistic, bounds, conflevel = 0.95, fixNaN, debug) {
  
  # 1 - Compute the mean estimate using the full sample
  
  mean.estimate <- statistic(full.sample)
  
  # 2 - Compute the Jacknife variance estimate from jackknife samples
  #     using Equation 5.3.29 from Gwet (2014, p.153).
  
  # Compute the statistic for each sample
  jackknife.values <- sapply(jackknife.samples, statistic)
  
  # Jacknife sample size = rater or item sample size
  n <- length(jackknife.values)
  
  # Remove NaN values (TODO: explain)
  if (sum(is.nan(jackknife.values)) != 0 && fixNaN) {
    jackknife.values <- ifelse(is.nan(jackknife.values), 1, jackknife.values)
    #jackknife.values <- jackknife.values[!is.na(jackknife.values)]    
  }
  jn <- length(jackknife.values)

  # Compute the jackknife mean estimate
  jackknife.mean.estimate <- mean(jackknife.values)
  
  # The following assumes an infinite population, so the (1 - gr) term disappears.
  variance <- (jn - 1) / jn * sum((jackknife.values - jackknife.mean.estimate) ^ 2)
  
  # 1 - Turn the mean and variance into a CI using the t distribution.
  #     This is not explained in Gwet's book but it is the method he uses in his
  #     R scripts. Also see Equation 7 from Abdi and Williams (2010, p.3).

  stderr <- sqrt(variance)
  
  if (debug) {
    cat("variance: ", variance, "\n")
    cat("stderr: ", stderr, "\n")
  }
  
  # lower confidence bound
  lower.bound <- mean.estimate - stderr * qt(1 - (1 - conflevel) / 2, n - 1)
  
  # upper confidence bound
  upper.bound <- mean.estimate + stderr * qt(1 - (1 - conflevel) / 2, n - 1)
  
  # Limit to range (Gwet does limits upper bound to 1 in his R scripts)
  lower.bound <- cap(lower.bound, bounds[1], bounds[2])
  upper.bound <- cap(upper.bound, bounds[1], bounds[2])  
  
  # Return the mean and confidence limits
  if(mean.estimate == 1){
      if(conflevel == 0.95) c(1, 1-3/n, 1) # This is the rule of 3: https://en.wikipedia.org/wiki/Rule_of_three_(statistics)
      else c(1, 1-4.61/n, 1) # Assuming a 99% CI
  }
  
  else c(mean.estimate, lower.bound, upper.bound)
}

#
# TODO: factorize with the above function
#
jackknife.CI.diff <- function(full.sample1, full.sample2, jackknife.samples1, jackknife.samples2, statistic, bounds, conflevel = 0.95, fixNaN, debug) {
  
  # 1 - Compute the mean estimate using the full sample
  
  mean.estimate <- statistic(full.sample1) - statistic(full.sample2)
    
  # 2 - Compute the Jacknife variance estimate from jackknife samples
  #     using Equation 5.3.29 from Gwet (2014, p.153).
  
  
  # Compute the statistic for each sample
  jackknife.values <- sapply(jackknife.samples1, statistic) - sapply(jackknife.samples2, statistic)
  
  ############ starting from here, the function is exactly the same
  
  # Jacknife sample size = rater or item sample size
  n <- length(jackknife.values)
  
  # Remove NaN values (TODO: explain)
  if (sum(is.nan(jackknife.values)) != 0 && fixNaN) {
    jackknife.values <- ifelse(is.nan(jackknife.values), 1, jackknife.values)
    #jackknife.values <- jackknife.values[!is.na(jackknife.values)]    
  }
  jn <- length(jackknife.values)
  
  # Compute the jackknife mean estimate
  jackknife.mean.estimate <- mean(jackknife.values)
  
  # The following assumes an infinite population, so the (1 - gr) term disappears.
  variance <- (jn - 1) / jn * sum((jackknife.values - jackknife.mean.estimate) ^ 2)
  
  # 1 - Turn the mean and variance into a CI using the t distribution.
  #     This is not explained in Gwet's book but it is the method he uses in his
  #     R scripts. Also see Equation 7 from Abdi and Williams (2010, p.3).
  
  stderr <- sqrt(variance)
  
  if (debug) {
    cat("variance: ", variance, "\n")
    cat("stderr: ", stderr, "\n")
  }
  
  # lower confidence bound
  lower.bound <- mean.estimate - stderr * qt(1 - (1 - conflevel) / 2, n - 1)
  
  # upper confidence bound
  upper.bound <- mean.estimate + stderr * qt(1 - (1 - conflevel) / 2, n - 1)
  
  # Limit to range (Gwet does limits upper bound to 1 in his R scripts)
  lower.bound <- cap(lower.bound, bounds[1], bounds[2])
  upper.bound <- cap(upper.bound, bounds[1], bounds[2])  
  
  # Return the mean and confidence limits
  c(mean.estimate, lower.bound, upper.bound)
}



#
# Lists all categories in a rating table
#
get.categories <- function(ratings) {
  unique(as.vector(as.matrix(ratings)))
}

#
# Cap value
#
cap <- function(x, xmin, xmax) {
  min(xmax, max(xmin, x))
}

#
# Display info on rating table
#
printRatingsInfo <- function(ratings) {
  r <- ncol(ratings)
  s <- nrow(ratings)
  k <- length(get.categories(ratings))
  cat("Rating data has ", s, " items, ", r, " raters, and ", k, " categories.\n", sep="")
}

#
# Display mean and CI in a human-readable fashion
#
rd <- function(x) {sprintf("%.3f", round(x,3))}
printCI <- function(name, result) {
  cat(name, " = ", rd(result[1]), ", 95% CI [", rd(result[2]), ", ", rd(result[3]), "]\n", sep="")
}




#
# Computes a CI on the difference between two datasets with the exact same raters: statistic(ratings1) - statistic(ratings2).
# ratings1 and ratings2 should include the same raters, in the same order.
#
boot.CI.diff.random.raters <- function(ratings1, ratings2, statistic, confint = 0.95, R = 3000, bounds = "auto") {
    
    if (bounds == "auto") {
        bounds <- statistic.bounds(statistic)
        bounds <- c(bounds[1] - bounds[2], bounds[2] - bounds[1])
    }
    
    boot.CI.diff(ratings1, ratings2, statistic, bounds, confint, R)
}



######
# Computes a 95% percentile bootstrap CI for a given statistic on a rating dataset, assuming a fixed set of
# items and random drawing from an infinite population of raters.
#
boot.CI.diff <- function(ratings1, ratings2, statistic, bounds, confinf = 0.95, R) {
    
    #if (bounds == "auto") {bounds <- statistic.bounds(statistic)}
    
    tratings1 <- t(ratings1)
    tratings2 <- t(ratings2)

    stat <- function(data, data2, indices) {
        d1 <- data[indices,]
        d2 <- data2[indices,]
        
        return(statistic(t(d1)) - statistic(t(d2)))
    }

    bootobject <- boot(data = tratings1, data2 = tratings2, statistic = stat, R=R)


    # Use the percentile method, as it seems to work better than others (e.g. bca)
    b <- tryCatch(boot.ci(bootobject, type='perc', confinf), error=function(...) list(fail=TRUE))
    
    if(length(b$fail) && b$fail) {
        ci <- c(0, 0, 0)
        #ci <- c(NA,NA)
        #warning('could not obtain BCa bootstrap confidence interval')
    } else {
        ci <- c(b$perc[4], b$perc[5])
    }
    
    lower.bound <- cap(ci[1], bounds[1], bounds[2])
    upper.bound <- cap(ci[2], bounds[1], bounds[2])
    
    # Return the mean and confidence limits
    c(statistic(ratings1) - statistic(ratings2), lower.bound, upper.bound)
    
    
    # Limit to theoretical bounds
    #ci[1] <- cap(ci[1], bounds[1], bounds[2])
    #ci[2] <- cap(ci[2], bounds[1], bounds[2])
    
    #  c(statistic(ratings1) - statistic(ratings2), ci[1], ci[2])
}


#===========================================================================================
# Statistical tests for between-participant Kappa comparisons with bootstrapping techniques
# Implementation for two only independent groups
# Theophanis Tsandilas
#======================================================================================
################## Assuming here that each group can have its own chance agreement

bootstrap.distribution.diff <- function(r, group1, group2, R){
    diffs <- replicate(R, fleiss.kappa.for.item(sample(group1, replace = TRUE), r) -  fleiss.kappa.for.item(sample(group2, replace = TRUE), r))
    
    diffs
}


# r is the indexes of the referents of interest
# group1 and group2 are two independent groups
fleiss.kappa.bootstrap.diff.ci <- function(r, group1, group2, R = 1000, plevel = .05){
    diff <- sort(bootstrap.distribution.diff(r, group1, group2, R))
    # uncomment this to check the distribution of bootstrap samples that we get...
    #hist(diff)
    
    perc <- R*plevel/2
    
    low <- floor(perc)
    upper <- ceiling(R - perc)
    
    point.estimate <- fleiss.kappa.for.item(group1, r) - fleiss.kappa.for.item(group2, r)
    
    c(point.estimate, diff[low], diff[upper])
}

