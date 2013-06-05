ISparamfit.guisdap <- function( acf , var , lags , fSite , aSite ,  initParam , mIon=c(30.5,16,1) ,
                       aprioriTheory , aprioriMeas , invAprioriCovar ,
                       paramLimits , fitFun=leastSquare.lvmrq , scaleFun=scaleParams ,
                               directTheory=ISdirectTheory.guisdap , ... ){
#
# Incoherent scatter plasma parameter
# 
# INPUT:
#  acf            a vector of complex autocorrelation fucntion values
#  var            variances of the acf values
#  lags           time-lags of the acf points
#  fSite          transmitter frequency
#  aSite          angle between incident and scattered wave vector, in degrees
#  initParam      start values for the fitting routine
#                   c( Ne , Ti , Te/Ti , nu_in , vi , pm2 , pm3 , ... )
#  aprioriTheory  theory matrix for the linear apriori theory
#  aprioriMeas    the imaginary apriori "measurements"
#  invAprioriCovar  inverse of the  covariance matrix of the apriori "measurements"
#  mIon           ion masses
#  paramLimits    2 x length(initParam) matrix of smallest (paramLimits[1,]) and largest (paramLimits[2,]) parameter values
#  fitFun         function used for the leas-squares optimisation
#  scaleFun       the function used for scaling step sizes in Levenberg-Marquardt iteration.
#  ...            additional input needed by scaleFun
#
# OUTPUT:
#  
#
#
# I. Virtanen 2012
#

  # number of data points
  nData <- length(acf)
  
  # check that acf and var are of the same length
  if(nData != length(var)) stop('different lengths of acf and variance vectors.')

  # use the same lag-values repeatedly if the vector is too short
  if(length(lags) < nData) lags <- rep(lags,length.out=nData)

  # create a frequency axis
  fmax <- 1/min(diff(unique(sort(lags))))
  fstep <- min(1/max(lags)/2,fSite/1e6)
  if(is.infinite(fmax)) fmax <- fSite/1e4
  if(is.infinite(fstep)) fstep <- fSite/1e6
  freq <- seq(-fmax,fmax,by=fstep)

  # frequency ambiguity functions
  nf <- length(freq)
  nl <- length(lags)
  fAmb <- matrix(nrow=nf,ncol=nl)
  for(k in seq(nl)) fAmb[,k] <- frequencyAmbiguity( lags[k] , freq )

  # the actual fitting
  nlsParam <- fitFun( measData=acf , measVar=var , initParam=initParam , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , paramLimits=paramLimits , fAmb=fAmb , directTheory=directTheory , scaleFun=scaleFun , mIon=mIon , fSite=fSite , aSite=aSite , xSite=freq , ... )

  return(nlsParam)

  
} # ISparamfit.guisdap
