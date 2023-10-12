dtGradient <- function( dirtheData , directTheory , param , ... ){
#
# Numerical measurement gradient
#
# INPUT:
#  dirtheData    direct theory data with the current parameter values
#  directTheory  direct theory function
#  param         current parameter values
#  ...           other parameters needed by directTheory
#
# OUTPUT:
#  length(param) x length(dirtheData) matrix of numerical derivatives
#

  nd <- length(dirtheData)
  np <- length(param)
  dtdiff <- matrix(0,nrow=np,ncol=nd)

  for(k in seq(length(param))){

    # a small difference to the current parameters
    param2     <- param
    param2[k]  <- param[k] + .0001

    # derivative
    dtdiff[k,] <- ( directTheory( param2 , ... ) - dirtheData ) / .0001
  }

  return(dtdiff)

}
