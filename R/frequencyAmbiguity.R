frequencyAmbiguity <- function( lag , freq ){
#
# Frequecy ambiguity function
# 
# INPUT:
#  lag     time lag
#  freq    frequency axis points
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

  return(exp(1i*2*pi*freq*lag)*s)

} # frequencyAmbiguity

