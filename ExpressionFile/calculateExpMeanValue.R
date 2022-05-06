{
  # function
  # average
  fReplicateMean <- function(x,names)
  {
    for (name in names)
    {
      itemName <- paste("average_",name, sep = "")
      #x[itemName] <- rowMeans(x[,regexpr(name,colnames(x))>0],na.rm=TRUE, dim=1)
      x[itemName] <- apply(x[,regexpr(name,colnames(x))>0], 1, trimMean)
    }
    averageNames <- paste("average_",names,sep="")
    x <- x[,c(averageNames)]
    return(x)
  }
  
  # trimMean
  trimMean <- function(vec,p1=0.2, p2=0.9) {
    #######################################################################
    #
    # This function computes the trimmed mean for a vector.  Note that this
    # implementation follows the Affymetrix code, which gives different
    # results than the standard R function mean().
    #
    # Arguments:
    #       vec - vector of values for computing trimmed mean
    #       p1,
    #       p2 - lower and upper percentage for trimming, expressed as a
    #            decimal fraction (not whole number)
    #
    # Value:
    #	a numeric value representing the trimmed mean for the given
    #       data vector
    #
    #######################################################################
    
    whole <- vec
    total <- length(vec)
    if (total==0) return(0)
    
    whole <- sort(whole)
    
    dG1 <- total * p1 + 1
    dG2 <- total * (1 - p2) + 1
    g1 <- floor(dG1)
    g2 <- floor(dG2)
    r1 <- dG1 - g1
    r2 <- dG2 - g2
    last <- total - g2 + 1
    if (last <= 0) last <- 0
    
    sum <- (1 - r1) * whole[g1] + (1 - r2) * whole[last]
    sum <- sum + sum(whole[(g1+1):(last-1)])
    
    subtotal <- last - g1 - 1
    if (subtotal <= 0) subtotal <- 0
    subtotal <- subtotal + 2 - r1 - r2
    return(sum / subtotal)
  }
}