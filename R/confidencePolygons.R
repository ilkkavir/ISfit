confidencePolygons <- function(d,h,limits){
# 
# I. Virtanen 2010
# 

  # sort the limits, because the largest polygon needs to be plotted first
  limits <- sort(unique(limits),decreasing=T)
  
  n <- length(limits)
  np <- length(d[1,])

  # a list for the confidence polygons
  confPols <- list()

  for (k in seq(n)){

    # confidence intervals at all ranges
    confInt <- matrix(ncol=2,nrow=dim(d)[1])
    for(r in seq(dim(d)[1])){
      confInt[r,] <- sort(d[r,])[c( (1-limits[k]) , (1+limits[k]) )*np/2]
    }
    confPols[[k]]     <- list()
    confPols[[k]]$x   <- c(confInt[,1],rev(confInt[,2]),confInt[1,1])
    confPols[[k]]$y   <- c(h,rev(h),h[1])
    confPols[[k]]$col <- rgb((1-k/(n+1)),(1-k/(n+1)),(1-k/(n+1)))

    # take care of NA values
    for(r in seq(2,length(confPols[[k]]$x))){
      if(is.na(confPols[[k]]$x[r])){
        confPols[[k]]$x[r] <- confPols[[k]]$x[r-1]
        confPols[[k]]$y[r] <- confPols[[k]]$y[r-1]
      }
    }

  }

  return(confPols)

} # confidencePolygons