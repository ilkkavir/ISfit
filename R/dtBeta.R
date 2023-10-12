dtBeta <- function( measData , measVar , dirtheGrad , dirtheData , param , aprioriTheory , aprioriMeas , invAprioriCovar){
#
# Gradient of the cost function, only real part to keep the output real
#
# INPUT:
#  measData        measured data
#  measVar         measurement variances
#  dirtheGrad      matrix of numerical differences
#  dirtheData      current direct theory values
#  param           parameter values
#  aprioriTheory   theory matrix of the linear apriori information
#  aprioriMeas     "apriori measurements"
#  invAprioriCovar inverse of apriori measurement covariance matrix
#
# I. Virtanen 2012
#
  return( (Re(dirtheGrad) %*% (Re(measData-dirtheData)/measVar)) + (Im(dirtheGrad) %*% (Im(measData-dirtheData)/measVar)) +
         t(aprioriTheory)%*%invAprioriCovar%*%aprioriMeas - (t(aprioriTheory)%*%invAprioriCovar%*%aprioriTheory)%*%param )
  
}
