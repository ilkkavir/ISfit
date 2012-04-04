scaleParams <- function(param,scale,inverse=F,...){
#
# Simple scaling to adjust step sizes in iteration
#
#
# I. Virtanen 2012
#

  if(inverse){
    return(param*scale)
  }else{
    return(param/scale)
  }
}
