heightTicks <- function(hlim){
# time ticks and corresponding labels that fit well for the time axis given by times [s]
  htick <- vector()
  hstr <- character()
  result <- list()

  htick <- pretty(hlim,n=5)
  htick <- htick[htick>hlim[1]]
  htick <- htick[htick<.95*hlim[2]]

  hstr <- as.character(htick)

  result[['tick']] <- htick
  result[['string']] <- hstr

  return(result)

} #timeTicks


