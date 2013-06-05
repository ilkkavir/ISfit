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

  # do not allow electron density step sizes smaller than 1e10
  parScales[1] <- max(parScales[1],1e10)

  # make sure that the collision frequency steps are non-zero
  parScales[6] <- max(parScales[6],10)

  # set ion abundance and acf scale steps to unity
  parScales[seq(10,length(parScales))] <- 1

  # set ion velocity steps to 100
  parScales[7:9] <- 100

  return(parScales)
  
}
