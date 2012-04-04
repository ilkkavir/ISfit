ISparamLimits.default <- function(nIon,nSite){
#
# Default lower and upper limits for parameters, in physical units
#
#
#
#
#
#
#

  nPar                            <- (7+8*nIon+nSite)

  paramLimits                     <- matrix(0,nrow=2,ncol=nPar)

  # velocities
  paramLimits[1,seq(5,nPar,by=8)] <- -1e4
  paramLimits[1,seq(6,nPar,by=8)] <- -1e4
  paramLimits[1,seq(7,nPar,by=8)] <- -1e4
  paramLimits[2,seq(5,nPar,by=8)] <- 1e4
  paramLimits[2,seq(6,nPar,by=8)] <- 1e4
  paramLimits[2,seq(7,nPar,by=8)] <- 1e4


  # electron and ion densities
  paramLimits[1,1]                <- 1e4
  paramLimits[2,1]                <- 1e13
  paramLimits[2,seq(9,nPar,by=8)] <- 1e13

  # temperatures
  paramLimits[1,seq(2,nPar,by=8)] <- 10
  paramLimits[1,seq(3,nPar,by=8)] <- 10
  paramLimits[2,seq(2,nPar,by=8)] <- 1e4
  paramLimits[2,seq(3,nPar,by=8)] <- 1e4

  # collision frequencies
  paramLimits[2,seq(4,nPar,by=8)] <- 1e20 # is this correct?

  # ion masses
  paramLimits[2,seq(8,nPar,by=8)] <- 200

  # acf scales
  paramLimits[1,seq((8*(nIon+1)),nPar)] <- 1e-2
  paramLimits[2,seq((8*(nIon+1)),nPar)] <- 100

  return(paramLimits)
g  
} # ISparamLimits.default
