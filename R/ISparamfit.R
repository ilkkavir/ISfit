ISparamfit <- function( acf , var , lags , iSite , fSite , aSite , kSite , B , initParam ,
                       aprioriTheory , aprioriMeas , invAprioriCovar , 
                       nIon , paramLimits , fitFun=leastSquare.lvmrq , scaleFun=scaleParams , directTheory=ISdirectTheory , ... ){
#
# Incoherent scatter plasma parameter fit
# 
# INPUT:
#  acf            a vector of complex autocorrelation fucntion values
#  var            variances of the acf values
#  lags           time-lags of the acf points
#  iSite          receiver site indices for all data points
#  fSite          transmitter frequencies for all receiver sites
#  aSite          angle between incident and scattered wave vector for each receiver site, in degrees
#  kSite          a list of scattering wave vectors at each site
#  B              magnetic field direction at target position
#  nIon           number of ion masses
#  fitFun         function used for the leas-squares optimisation
#
#  initParam      initial values for the fitting routine, the actual format
#                 depends on requirements of directTheory
#  aprioriTheory  theory matrix for the linear apriori theory, must match with directTheory
#  aprioriMeas    the imaginary apriori "measurements", must match with directTheory
#  invAprioriCovar  inverse of the  covariance matrix of the apriori "measurements"
#  paramLimits    2 x length(initParam) matrix of smallest (paramLimits[1,]) and
#                 largest (paramLimits[2,]) parameter values. Must match with directTheory
#  scaleFun       the function used for scaling step sizes in Levenberg-Marquardt iteration. Must
#                 match with directTheory
#  directTheory   direct theory function for the fitting routine
#  ...            additional input needed by scaleFun
#
# OUTPUT:
#  
#   Output from fitFun, contents vary depending on fitFun and directTheory
#
# I. Virtanen 2012, 2013
#

  # number of data points
  nData <- length(acf)
  
  # check that acf and var are of the same length
  if(nData != length(var)) stop('different lengths of acf and variance vectors.')

  # use the same lag-values repeatedly if the vector is too short
  if(length(lags) < nData) lags <- rep(lags,length.out=nData)

  # use the same site indeces repeatedly if the vector is too short
  if(length(iSite) < nData) iSite <- rep(iSite,length.out=nData)
  
  # count the receiver sites
  ns <- length(unique(iSite))

  # repeat the carrier frequencies and scattering angles if necessary
  if(length(fSite) < ns) fSite <- rep(fSite,length.out=ns)
  if(length(aSite) < ns) aSite <- rep(aSite,length.out=ns)

  # create frequency axes for all sites
  freq <- list()
  for(k in seq(ns)){
      fmax <- 4/min(diff(unique(sort(lags[which(iSite==k)]))))
      fstep <- min(1/max(lags[which(iSite==k)])/4,fSite[k]/1e6)
      if(is.infinite(fmax)) fmax <- fSite[k]/1e4
      if(is.infinite(fstep)) fstep <- fSite[k]/1e6
      if(fmax<1e3) fmax <- fSite[k]/1e4
      if(fstep>=fmax) fstep <- fmax / 100
      freq[[k]] <- seq(-fmax,fmax,by=fstep)
  }

  # frequency ambiguity functions
  nf <- max(sapply(freq,length))
  nl <- length(lags)
  fAmb <- matrix(nrow=nf,ncol=nl)
  for(k in seq(nl)) fAmb[,k] <- frequencyAmbiguity( lags[k] , freq[[iSite[k]]] )

  # the actual fitting
  nlsParam <- fitFun( measData=acf , measVar=var , initParam=initParam , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , paramLimits=paramLimits , fAmb=fAmb , scaleFun=scaleFun , nIon=nIon , nSite=ns , iSite=iSite , fSite=fSite , aSite=aSite , kSite=kSite , B=B , xSite=freq , directTheory=directTheory , ... )

  return(nlsParam)

  
}

