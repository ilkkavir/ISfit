chiSquare <- function( measData , measVar , dirtheData , param , aprioriTheory , aprioriMeas , invAprioriCovar , ... ){
#
# sum of residuals 
#
#
# INPUT:
#  measData         a vector of measured data values
#  measVar          variances of the measurements
#  dirtheData       data values from direct theory
#  param            current parameter values
#  aprioriTheory    theory matrix of linear apriori information
#  aprioriMeas      imaginary "apriori measurements"
#  invAprioriCovar  inverse of apriori measurement covariance matrix (i.e. apriori Fisher information matrix)
#
# OUTPUT:
#  chi-squared with the current parameters
#  
# I. Virtanen 2012
#
  return( c( sum( abs(measData - dirtheData)**2 / measVar )  +
         ( t( aprioriTheory%*%param - aprioriMeas ) %*% invAprioriCovar %*% ( aprioriTheory%*%param - aprioriMeas ) ) ) )

}
