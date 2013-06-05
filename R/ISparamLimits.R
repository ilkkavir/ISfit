ISparamLimits <- function(nIon,nSite){
#
# Default lower and upper limits for parameters, in physical units
#
#
#
#
#
#
#

  nPar                            <- (9+nIon+nSite)

  paramLimits                     <- matrix(0,nrow=2,ncol=nPar)

  # velocities
  paramLimits[1,7:9] <- -1e4
  paramLimits[2,7:9] <-  1e4


  # electron density
  paramLimits[1,1]                <- 1e8
  paramLimits[2,1]                <- 1e13

  # ion temperatures
  paramLimits[1,2:3] <- 10
  paramLimits[2,2:3] <- 1e4

  # electron temperatures
  paramLimits[1,4:5] <- 10
  paramLimits[2,4:5] <- 5e4

  # collision frequency
  paramLimits[2,6] <- 1e20 # is this correct?

  # ion abundaces
  paramLimits[2,10:(9+nIon)] <- 1

  # acf scales
  paramLimits[1,(10+nIon):nPar] <- 1e-2
  paramLimits[2,(10+nIon):nPar] <- 100

  return(paramLimits)

}
