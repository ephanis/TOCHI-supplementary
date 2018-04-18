#  							AGREE.COEFF3.RAW.R
#						 		 (November 8, 2014)
#Description: This script file contains a series of R functions for computing various agreement coefficients
#		  for multiple raters (2 or more) when the input data file is in the form of nxr matrix or data frame showing 
#             the actual ratings each rater (column) assigned to each subject (in row). That is n = number of subjects, and r = number of raters.
#             A typical table entry (i,g) represents the rating associated with subject i and rater g. 
#Author: Kilem L. Gwet, Ph.D.
#
#

# ==============================================================
# This is an r function for trimming leading and trealing blanks
# ==============================================================
trim <- function( x ) {
  gsub("(^[[:space:]]+|[[:space:]]+$)", "", x) }


#
# Gwet's original implementation accounts for weighted items. Weigths are not relevenat in gesture elicitation, and so always use an identity matrix for weights.
#
identity.weights<-function(categ){
    weights<-diag(length(categ))
    return (weights)
}


#
# Computes only percent agreement
#
percent.agreement.raw.noinference <- function(ratings, conflev=0.95, N=Inf) {
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)
  
  weights.mat<-identity.weights(categ)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else
      agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pa
}


#
# Computes Fleiss' chance agreement
#
fleiss.chance.agreement.raw.noinference <- function(ratings, conflev=0.95, N=Inf) {
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
    categ <- sort(as.vector(na.omit(categ.init)))
  }else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)
  
  weights.mat<-identity.weights(categ)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else
      agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
    
  # calculating Fleiss's chance agreement
  ri.vec <- agree.mat%*%rep(1,q)  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  
  pe
}


#=====================================================================================
#fleiss.kappa.raw: This function computes Fleiss' generalized kappa coefficient (see Fleiss(1971)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Fleiss, J. L. (1981). Statistical Methods for Rates and Proportions. John Wiley & Sons.
#======================================================================================
fleiss.kappa.raw.noinference <- function(ratings, conflev=0.95, N=Inf){
  pa <- percent.agreement.raw.noinference(ratings, conflev, N)
  pe <- fleiss.chance.agreement.raw.noinference(ratings, conflev, N)

  fleiss.kappa <- (pa-pe)/(1-pe)
  
  fleiss.kappa
}


#=====================================================================================
#krippen.alpha.raw: This function computes Krippendorff's alpha coefficient (see Krippendorff(1970, 1980)) and 
#		   its standard error for 3 raters or more when input dataset is a nxr matrix of alphanumeric 
#		   ratings from n subjects and r raters.
#-------------
#The algorithm used to compute krippendorff's alpha is very different from anything that was published on this topic. Instead,
#it follows the equations presented by K. Gwet (2010)
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Gwet, K. (2012). Handbook of Inter-Rater Reliability: The Definitive Guide to Measuring the Extent of Agreement Among 
#	Multiple Raters, 3rd Edition.  Advanced Analytics, LLC; 3rd edition (March 2, 2012)
#Krippendorff (1970). "Bivariate agreement coefficients for reliability of data." Sociological Methodology,2,139-150
#Krippendorff (1980). Content analysis: An introduction to its methodology (2nd ed.), New-bury Park, CA: Sage.
#======================================================================================
krippen.alpha.raw.noinference <- function(ratings, conflev=0.95, N=Inf){
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init)){
    categ <- sort(as.vector(na.omit(categ.init)))
  }else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)
  
  weights.mat<-identity.weights(categ)
  
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else
      agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  # calculating krippendorff's alpha coefficient
  
  ri.vec <- agree.mat%*%rep(1,q)
  agree.mat<-agree.mat[(ri.vec>=2),]
  agree.mat.w <- agree.mat.w[(ri.vec>=2),]
  ri.vec <- ri.vec[(ri.vec>=2)]
  ri.mean <- mean(ri.vec)
  n <- nrow(agree.mat)
  epsi <- 1/sum(ri.vec)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  pa <- (1-epsi)* sum(sum.q/(ri.mean*(ri.vec-1)))/n + epsi
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/ri.mean))
  pe <- sum(weights.mat * (pi.vec%*%t(pi.vec)))
  krippen.alpha <- (pa-pe)/(1-pe)
  
  krippen.alpha
}


#===========================================================================================
#bp.coeff.raw: Brennan-Prediger coefficient (see Brennan & Prediger(1981)) and its standard error for multiple raters when input 
#		   dataset is a nxr matrix of alphanumeric ratings from n subjects and r raters 
#-------------
#The input data "ratings" is a nxr data frame of raw alphanumeric ratings
#from n subjects and r raters. Exclude all subjects that are not rated by any rater.
#Bibliography:
#Brennan, R.L., and Prediger, D. J. (1981). ``Coefficient Kappa: some uses, misuses, and alternatives."
#           Educational and Psychological Measurement, 41, 687-699.
#======================================================================================
bp.coeff.raw.noinference <- function(ratings, conflev=0.95, N=Inf){
  ratings.mat <- as.matrix(ratings) 
  n <- nrow(ratings.mat) # number of subjects
  r <- ncol(ratings.mat) # number of raters
  f <- n/N # final population correction 
  
  # creating a vector containing all categories used by the raters
  
  categ.init <- unique(as.vector(ratings.mat))
  if (is.numeric(categ.init))
    categ <- sort(as.vector(na.omit(categ.init)))
  else{
    categ.init <- trim(categ.init) #trim vector elements to remove leading and trailing blanks
    categ <- categ.init[nchar(categ.init)>0]
  }
  q <- length(categ)
  
  weights.mat<-identity.weights(categ)
  
 
  # creating the nxq agreement matrix representing the distribution of raters by subjects and category
  agree.mat <- matrix(0,nrow=n,ncol=q)
  for(k in 1:q){
    if (is.numeric(ratings.mat)){
      k.mis <-(ratings.mat==categ[k])
      in.categ.k <- replace(k.mis,is.na(k.mis),FALSE)
      agree.mat[,k] <- in.categ.k%*%rep(1,r) 
    }else
      agree.mat[,k] <- (trim(ratings.mat)==categ[k])%*%rep(1,r)
  }
  
  agree.mat.w <- t(weights.mat%*%t(agree.mat))
  
  
  ri.vec <- agree.mat%*%rep(1,q)
  sum.q <- (agree.mat*(agree.mat.w-1))%*%rep(1,q)
  n2more <- sum(ri.vec>=2)
  pa <- sum(sum.q[ri.vec>=2]/((ri.vec*(ri.vec-1))[ri.vec>=2]))/n2more
  
  pi.vec <- t(t(rep(1/n,n))%*%(agree.mat/(ri.vec%*%t(rep(1,q)))))
  pe <- sum(weights.mat) / (q^2)
  bp.coeff <- (pa-pe)/(1-pe)
  
  bp.coeff
}



#
# Computes Fleiss' chance agreement for indivindual items
fleiss.kappa.for.item <- function(ratings, item.index) {
    # Compute Fleiss' chance agreement
    pe <- fleiss.chance.agreement.raw.noinference(ratings)
    
    # Compute Fleiss kappa for the specific item
    pa <- percent.agreement.raw.noinference(ratings[item.index,], )
    
    (pa - pe) / (1 - pe)
}

