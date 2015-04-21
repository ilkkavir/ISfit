frequencyAmbiguity <- function( lag , freq , flen=0){
#
# Frequecy ambiguity function
# 
# INPUT:
#  lag     time lag
#  freq    frequency axis points
#  flen    filter impulse response length
#
# OUTPUT:
#  freqAmb frequency ambiguity function
#
#

  nf <- length(freq)
  
  # differences between frequency points
  fdiff <- diff(freq)

  # scale factor according ot frequency differencees
  s     <- c(fdiff[1],(fdiff[1:(nf-2)]+fdiff[2:(nf-1)])/2,fdiff[nf-1])

  # we assume a boxcar impulse response, so frequency response will be a sinc function
  sincarg <- freq * flen
  sinc <- sin(sincarg)/sincarg
  sinc[sincarg==0] <- 1

  return(exp(1i*2*pi*freq*lag)*sinc*s)

} # frequencyAmbiguity

