ISparamScales.general <- function(param,nIon){
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
  parScales[seq(4,length(parScales),by=8)] <- pmax(parScales[seq(4,length(parScales),by=8)],10)

#  # ion density steps should match with electron density steps
#  parScales[seq(9,length(parScales),by=8)] <- parScales[1]
  # ion abundance step is unity
  parScales[seq(9,length(parScales),by=8)] <- 1

  # set ion velocity steps to 100
  parScales[seq(5,length(parScales),by=8)] <- 100
  parScales[seq(6,length(parScales),by=8)] <- 100
  parScales[seq(7,length(parScales),by=8)] <- 100

  # set the acf scale step sizes to unity
  parScales[seq((8*(nIon+1)),length(parScales))] <- 1

  return(parScales)
  
} # ISparamScales.general
