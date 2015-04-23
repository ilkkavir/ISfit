ISparamScales <- function(param,nIon){
#
#
#
#
#
#
#
#
#
#
#
#

  parScales <- param

  # do not allow electron density step sizes smaller than 1e11
  parScales[1] <- max(parScales[1],1e11)

  # do not allow temperature steps smaller than 300 K, this
  # is mainly to keep reasonable width for the prior distributions
#  parScales[2:5] <- pmax(parScales[2:5],300)
  parScales[2:5] <- max(parScales[2:5],300)

  # make sure that the collision frequency steps are non-zero
  parScales[6] <- max(parScales[6],10)

  # set ion abundance and acf scale steps to unity
  parScales[seq(10,length(parScales))] <- 1

  # set ion velocity steps to 1000
  parScales[7:9] <- 1000

  return(parScales)
  
}
