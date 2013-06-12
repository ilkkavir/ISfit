ISapriori <- function( aprioriParam ,  nIon , absCalib=FALSE , TiIsotropic=FALSE , ... ){
#
# default apriori theory matrix, "measurements", and covariance matrix for 3D IS parameter fits
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
  nApriori                     <- ( nPar + 3 )
  
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
  aprioriStd[2]                <- .1                 # paralell ion temperature
  aprioriStd[3]                <- 1                  # perpendicular ion temperature
  aprioriStd[4]                <- .1                 # parallel electron temperature
  aprioriStd[5]                <- 1                  # perpendicular electron temperature
  aprioriStd[6]                <- 1e-3               # ion-neutral collision frequency
  aprioriStd[7]                <- 1                  # ion velocity, x-component
  aprioriStd[8]                <- 1                  # ion velocity, y-component
  aprioriStd[9]                <- 1                  # ion velocity, z-component
  aprioriStd[10:(9+nIon)]      <- 1e-3               # ion abundances

  aprioriStd[nIon+10]           <- 1e-3                  # the first site is a reference, do not scale
  if(nPar>(nIon+10)){
    if(absCalib){
      aprioriStd[(nIon+11):length(aprioriParam)] <- 1e-3 # fix all sites to the same ACF scale
    }else{
      aprioriStd[(nIon+11):length(aprioriParam)] <- .1   # allow scaling for other sites
    }
  }

  # force certain parameter differences close to zero
  curRow                         <- nPar + 1
  
  # assume isotropic electron temperature
  aprioriTheory[curRow,c(4,5)]   <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1
  
  # ion temperature anisotropy
  aprioriTheory[curRow,c(2,3)]   <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  if(TiIsotropic){
    aprioriStd[curRow]             <- 1e-3
  }else{
    aprioriStd[curRow]             <- .1
  }
  curRow                         <- curRow + 1

  # Sum of ion abundances must be unity
  aprioriTheory[curRow,10:(nIon+9)] <- 1
  aprioriMeas[curRow] <- 1
  aprioriStd[curRow] <- 1e-3
  curRow <- curRow+1


  return(list(aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas))
  
}
