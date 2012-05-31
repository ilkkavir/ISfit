ISapriori.default.3D <- function( aprioriParam , nIon , absCalib=FALSE , ... ){
#
# default apriori theory matrix, "measurements", and covariance matrix for 3D IS parameter fits
#
# INPUT:
#  aprioriParam    apriori parameter values
#  nIon            number of ions in the parameter vector
#  absCalib        if TRUE, the different sites are assumed to be absolutely calibrated and ACF scales are fixed to unity for all sites
#                  if FALSE, only the first site is fixed and others are allowed to scale in the parameter fit
#  ...             arbitrary parameters to be passed forward to other functions, mainly for compatability reasons
#
# OUTPUT:
#  aprioriTheory      apriori theory matrix
#  aprioriMeas        apriori "measurements"
#  invAprioriCovar    inverse of apriori covariance matrix
#
#  I. Virtanen 2012
#
  # length of the parameter vector
  nPar                         <- length(aprioriParam)

  # number of imaginary apriori "measurements"
  nApriori                     <- ( nPar + (nIon-1)*5 + 5 )
  
  # apriori theory matrix
  aprioriTheory                <- matrix( 0 , nrow=nApriori , ncol=nPar )

  # apriori measurement vector
  aprioriMeas                  <- vector(mode='numeric',length=nApriori)

  # the apriori covariance matrix will be diagonal, so we begin with a vector of standard deviations, which is easier.
  aprioriStd                   <- vector(mode='numeric',length=nApriori)

  
  # apriori parameter values
  aprioriTheory[1:nPar,1:nPar] <- diag(rep(1,nPar))
  
  aprioriMeas[1:nPar]          <- aprioriParam
  
  aprioriStd[1]                <- .5                # electron density
  aprioriStd[2]                <- .5                # parallel electron temperature
  aprioriStd[3]                <- .5                # perpendicular electron temperature
  aprioriStd[4]                <- 1e-3              # electron-neutral collision frequency
  aprioriStd[5]                <- 1e4               # electron velocity, x-component
  aprioriStd[6]                <- 1e4               # electron velocity, y-component
  aprioriStd[7]                <- 1e4               # electron velocity, z-component
  aprioriStd[8]                <- 1e-3              # mass of first ion species
  aprioriStd[9]                <- 1e-3              # abundance of first ion species
  aprioriStd[10]               <- .5                # parallel temperature of first ion species
  aprioriStd[11]               <- .5                # perpendicular temperature of first ion species
  aprioriStd[12]               <- 1e-3              # ion-neutral collision frequency for the first species
  aprioriStd[13]               <- 1e3               # velocity, x-component
  aprioriStd[14]               <- 1e3               # velocity, y-component
  aprioriStd[15]               <- 1e3               # velocity, z-component

  if(nIon>1){
    aprioriStd[16:(nIon*8+7)]  <- rep(
                                         c(
                                           1e-3,       # ion mass
                                           1e-3,       # ion abundance
                                           .5 ,        # parallel ion temperature
                                           .5 ,        # perpendicular ion temperature
                                           1e-3  ,     # ion-neutral collision frequency
                                           1e3 ,       # velocity, x-component
                                           1e3 ,       # velocity, y-component
                                           1e3         # velocity, z-component
                                           ),
                                      (nIon-1)
                                         )
  }
  aprioriStd[(nIon+1)*8]         <- 1e-3                  # the first site is a reference, do not scale
  if(nPar>(nIon+1)*8){
    if(absCalib){
      aprioriStd[((nIon+1)*8+1):length(aprioriParam)] <- 1e-3 # fix all sites to the same ACF scale
    }else{
      aprioriStd[((nIon+1)*8+1):length(aprioriParam)] <- .5   # allow scaling for other sites
    }
  }
  
  # apriori information of certain parameter differences
  curRow                         <- nPar + 1
  
  # electron temperature is assumed to be isotropic
  aprioriTheory[curRow,c(2,3)]   <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1

#  # ion temperature anisotropy
#  aprioriTheory[curRow,c(10,11)] <- c(1,-1)
#  aprioriMeas[curRow]            <- 0
#  aprioriStd[curRow]             <- 1000
#  curRow                         <- curRow + 1

  # sum of all ion abundances must be 1
#  aprioriTheory[curRow,1]        <- 1
  aprioriTheory[curRow,seq(9,(8*nIon+7),by=8)] <- 1
  aprioriMeas[curRow]            <- 1
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1

  # all ion species are in same temperature
  if(nIon>1){
    for(k in seq(nIon-1)){
      aprioriTheory[curRow,c(10,(10+k*8))] <- c(1,-1)
      aprioriMeas[curRow]                  <- 0
      aprioriStd[curRow]                   <- 1e-3
      curRow                               <- curRow + 1
      aprioriTheory[curRow,c(11,(11+k*8))] <- c(1,-1)
      aprioriMeas[curRow]                  <- 0
      aprioriStd[curRow]                   <- 1e-3
      curRow                               <- curRow + 1
    }
  }

  # electrons and ions usually have the same velocity
  aprioriTheory[curRow,c(5,13)]  <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1
  aprioriTheory[curRow,c(6,14)]  <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1
  aprioriTheory[curRow,c(7,15)]  <- c(1,-1)
  aprioriMeas[curRow]            <- 0
  aprioriStd[curRow]             <- 1e-3
  curRow                         <- curRow + 1

  # all ions have have the same velocity
  if(nIon>1){
    for(k in seq(nIon-1)){
      aprioriTheory[curRow,c(13,(13+k*8))] <- c(1,-1)
      aprioriMeas[curRow]                  <- 0
      aprioriStd[curRow]                   <- 1e-3
      curRow                               <- curRow + 1      
      aprioriTheory[curRow,c(14,(14+k*8))] <- c(1,-1)
      aprioriMeas[curRow]                  <- 0
      aprioriStd[curRow]                   <- 1e-3
      curRow                               <- curRow + 1      
      aprioriTheory[curRow,c(15,(15+k*8))] <- c(1,-1)
      aprioriMeas[curRow]                  <- 0
      aprioriStd[curRow]                   <- 1e-3
      curRow                               <- curRow + 1      
    }
  }

  aprioriTheory <<- aprioriTheory
  aprioriStd    <<- aprioriStd
  aprioriMeas   <<- aprioriMeas
  return(list(aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas))
  
} # ISapriori.default.3D
