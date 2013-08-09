scaleCovar <- function(covar,scale,inverse=F,...){
#
# Simple scaling to adjust step sizes in iteration
#
#
# I. Virtanen 2013
#

  if(inverse){
    return(t(t(covar*scale)*scale))
  }else{
    return(t(t(covar/scale)/scale))
  }
}
