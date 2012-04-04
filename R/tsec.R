tsec <- function(sec){
  h <- sec%%(3600*24)
  m <- sec%%3600
  s <- sec%%60
  d <- (sec-h)/(3600*24)
  h <- (h-m)/3600
  m <- (m-s)/60
  return(c(d,h,m,s))
} # tsec

