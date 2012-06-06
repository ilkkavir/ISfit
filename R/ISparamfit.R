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
#  B              magnetic field direction
#  initParam      start values for the fitting routine
#                   c( Ne , Te_par , Te_perp , nu_en , ve_x , ve_y , ve_z ,
#                 m1 , N1 , T1_par , T1_perp , nu_1n , v1_x , v1_y , v1_z ,
#                 m2 , N2 , T2_par , ...
#                 ...
#                 mnIon , ...                                   , vnIon_z ,
#                 s_1 , s_2 , ... , s_Nsite)
#                 s_1 ... s_Nsite are acf scaling factor for different sites
#  aprioriTheory  theory matrix for the linear apriori theory
#  aprioriMeas    the imaginary apriori "measurements"
#  invAprioriCovar  inverse of the  covariance matrix of the apriori "measurements"
#  nIon           number of ion masses
#  paramLimits    2 x length(initParam) matrix of smallest (paramLimits[1,]) and largest (paramLimits[2,]) parameter values
#  fitFun         function used for the leas-squares optimisation
#  scaleFun       the function used for scaling step sizes in Levenberg-Marquardt iteration.
#  direcTheory    direct theory function for the fitting routine
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

  # use the same site indeces repeatedly if the vector is too short
  if(length(iSite) < nData) iSite <- rep(iSite,length.out=nData)
  
  # count the receiver sites
  ns <- length(unique(iSite))

  # repeat the carrier frequencies and scattering angles if necessary
  if(length(fSite) < ns) fSite <- rep(fSite,length.out=ns)
  if(length(aSite) < ns) aSite <- rep(aSite,length.out=ns)

  # create frequency axes for all sites
  freq <- list()
  for(k in seq(ns)) freq[[k]] <- seq(-100000,100000,by=1000)*fSite[k]/1e9

  # frequency ambiguity functions
  nf <- length(freq[[1]])
  nl <- length(lags)
  fAmb <- matrix(nrow=nf,ncol=nl)
  for(k in seq(nl)) fAmb[,k] <- frequencyAmbiguity( lags[k] , freq[[iSite[k]]] )

  # the actual fitting
  nlsParam <- fitFun( measData=acf , measVar=var , initParam=initParam ,
                      aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar ,
                      paramLimits=paramLimits , fAmb=fAmb ,
                      scaleFun=scaleFun , nIon=nIon , nSite=ns , iSite=iSite , fSite=fSite , aSite=aSite ,
                      kSite=kSite , B=B , xSite=freq , directTheory=directTheory , ... )

  return(nlsParam)

  
} # ISparamfit
