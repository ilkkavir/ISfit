ISapriori.default.guisdap <- function( aprioriParam ,  ... ){
#
# default apriori theory matrix, "measurements", and covariance matrix for guisdap-style IS parameter fits
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
#  I. Virtanen 2012
#
  # length of the parameter vector
  nPar                         <- length(aprioriParam)

  # apriori theory matrix
  aprioriTheory                <- diag(rep(1),nPar)

  # apriori measurement vector
  aprioriMeas                  <- vector(mode='numeric',length=nPar)

  # the apriori covariance matrix will be diagonal, so we begin with a vector of standard deviations, which is easier.
  aprioriStd                   <- vector(mode='numeric',length=nPar)

  
  # apriori parameter values
  aprioriMeas[1:nPar]          <- aprioriParam
  
  aprioriStd[1]                <- 1e5                # electron density
  aprioriStd[2]                <- 1                  # ion temperature
  aprioriStd[3]                <- 1                  # temperature ratio
  aprioriStd[4]                <- 1e-3               # ion-neutral collision frequency
  aprioriStd[5]                <- 1e4                # ion velocity
  aprioriStd[6:nPar]           <- 1e-3

  aprioriTheory <<- aprioriTheory
  aprioriStd    <<- aprioriStd
  aprioriMeas   <<- aprioriMeas
  return(list(aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas))
  
} # ISapriori.default.3D
