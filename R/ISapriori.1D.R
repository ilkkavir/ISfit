ISapriori.1D <- function( aprioriParam ,  nIon , ... ){
#
# default apriori theory matrix, "measurements", and covariance matrix for 1D IS parameter fits
#
# INPUT:
#  aprioriParam    apriori parameter values
#  ...             arbitrary parameters to be passed forward to other functions, mainly for compatability reasons
#
# OUTPUT:
#  aprioriTheory      apriori theory matrix
#  aprioriMeas        apriori "measurements"
#  invAprioriCovar    inverse of apriori covariance matrix
#
#  I. Virtanen 2012, 2013
#
  # length of the parameter vector
  nPar                         <- length(aprioriParam)

  # number of imaginary apriori "measurements"
  nApriori                     <- ( nPar + 5 )
  
  # apriori theory matrix
  aprioriTheory                <- matrix( 0 , nrow=nApriori , ncol=nPar )

  # apriori measurement vector
  aprioriMeas                  <- vector(mode='numeric',length=nApriori)

  # the apriori covariance matrix will be diagonal, so we begin with
  # a vector of standard deviations, which is easier.
  aprioriStd                   <- vector(mode='numeric',length=nApriori)
  
  # apriori parameter values
  aprioriTheory[1:nPar,1:nPar] <- diag(rep(1,nPar))
  
  aprioriMeas[1:nPar]          <- aprioriParam
  
  aprioriStd[1]                <- 1e5                # electron density
  aprioriStd[2]                <- 1                  # paralell ion temperature
  aprioriStd[3]                <- 1                  # perpendicular ion temperature
  aprioriStd[4]                <- 1                  # parallel electron temperature
  aprioriStd[5]                <- 1                  # perpendicular electron temperature
  aprioriStd[6]                <- 1e-3               # ion-neutral collision frequency
  aprioriStd[7]                <- 1e4                # ion velocity, x-component
  aprioriStd[8]                <- 1e4                # ion velocity, y-component
  aprioriStd[9]                <- 1e4                # ion velocity, z-component
  aprioriStd[10:(9+nIon)]       <- 1e-3              # ion abundances
  aprioriStd[(10+nIon):nPar]    <- 1e-3               # fix the ACF scaling (there should be only one site)
  
  # force certain parameter differences close to zero
  curRow                         <- nPar + 1
  
  # assume isotropic electron temperature
  aprioriTheory[curRow,c(4,5)]   <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1

  # assume isotropic ion temperature
  aprioriTheory[curRow,c(2,3)] <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1

  # assume that all velocity components have equal amplitude
  # this is not true, but allows proper estimation of beam-aligned velocity component
  aprioriTheory[curRow,c(7,8)]   <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1
  aprioriTheory[curRow,c(7,9)] <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1

  # Sum of ion abundances must be unity
  aprioriTheory[curRow,10:(nIon+9)] <- 1
  aprioriMeas[curRow] <- 1
  aprioriStd[curRow] <- 1e-3
  curRow <- curRow+1

  aprioriTheory <- aprioriTheory
  aprioriStd    <- aprioriStd
  aprioriMeas   <- aprioriMeas
  return(list(aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas))
  
}
