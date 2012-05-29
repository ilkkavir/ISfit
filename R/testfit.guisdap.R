testfit.guisdap <- function(refPoint  = KIR,
                       hTarg     = 200,
                       vion      = 0,
                       fradar    = 235e6,
                       aSite     = 180,
                       nacf      = 50,
                       lags      = seq(50)*5e-5,
                       dataStd   = 1e-19,

                       plotTest  = F,
                       absLimit  = 5,
                       diffLimit = 1e-2,
                       maxLambda = 1e30,
                       maxIter   = 100,

                       time      = c(2009,7,1,11,0,0),
                       heights   = seq(1,1000),
                       initErr   = c(1e10,100,.3,0,0,0,0)
                       ){
#
# Test the guisdap-style plasma parameter fit with simulated ACF data
#
# INPUT:
#   refPoint   c(lat,lon) reference point, above which iri-parameters are calculated
#   hTarg      target height above the reference point
#   vion       ion velocity 
#   fradar     transmitter frequency [Hz]
#   aSite      scattering angle
#   nacf       number of autocorrelation functions (integration periods) 
#   lags       time-lags  [us]
#   dataStd    lag profile standard deviation
#
#   plotTest   logical, if T, the measurements and direct theory values are plotted at each iteration step
#   absLimit   upper limit for fit residual, when the residual is reduced below absLimit, the solver starts to look for end condition
#              based on diffLimit
#   diffLimit  upper limit for the fractional change of residual in one iteration step. The iteration is stopped when residual
#              is below absLimit and the fractional change of residual is below diffLimit
#   maxLambda  Maximum value for the lambda parameter in levenberg-marquardt iteration
#   maxIter    Maximum number of iterations (Hess matrix calculations) in the levenberg-marquardt iteration
#
#   time       c(year,month,day,hour,min,sec) time for which plasma parameters are calculated with the IRI model
#   heibeg, heiend, heistp  parameters for the IRI-model call. The model parameters are originally computed at heights
#              seq(heibeg,heiend,by=heistp). Most probably the user does not need to change these parameters
#   initErr    Error added to initial parameters before solving
#
#
# I. Virtanen 2012  
#

  mIon <- c(30.5,16,1)

  aSite <- 180

  # parameters from iri model
  ptmp           <- iriParams(time=time,latitude=refPoint[1],longitude=refPoint[2],heights=heights)

  # height point closest to the user input value
  h              <- which(abs(heights-hTarg)==min(abs(heights-hTarg)))[1]

  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
  # the densities in outfmsis are in cm^-3
  ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,h])['NO+',] )
#  ioncoll        <- (2.44e-16*ptmp$outfmsis[2,h] + 4.34e-16*ptmp$outfmsis[3,h] + 4.28e-16*ptmp$outfmsis[4,h])*1e6
  
  # parameters
#  par            <- c( ptmp$outf[c(1,3),h] , ptmp$outf[4,h]/ptmp$outf[3,h] , ioncoll , vion , ptmp$outf[5,h]/100 , ptmp$outf[6,h]/100)
  par            <- c( ptmp['e-',h] , ptmp['Ti',h] , ptmp['Te',h]/ptmp['Ti',h] , ioncoll , vion , ptmp['O+',h]/ptmp['e-',h] , ptmp['H+',h]/ptmp['e-',h] )


  # initial (and apriori) parameter values
  initParam      <- par + initErr
  initParam[5]   <- 0
  initParam[3]   <- 1

  # parameter scaling factors
  parScales      <- par
  parScales[3]   <- 1
  parScales[5:7]  <- 1

  # scale the initial parameter values
  initParam      <- scaleParams( initParam , parScales , inverse=F)

  # parameter value limits
  parLimits      <- matrix( c( 1e4 , 10 , 1e-2 , 0 , -1e4 , 0 , 0    ,     1e13 , 1e4 , 100 , 1e20 , 1e4 , 1 , 1) , nrow=2 , byrow=T )

  # scale the parameter limits
  limitParam     <- parLimits
  limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
  limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
  
  # apriori information
  apriori        <- ISapriori.default.guisdap( initParam , nIon )

  #
  # generate the simulated ACF data and other
  #
  nData    <- length(lags) * nacf
  freq     <- seq(-100000,100000,by=1000)*fradar/1e9
  ele <- c(par[1],par[2]*par[3],par[2]*par[3],par[4]*.35714,par[5],0,0)
  ion <-  list(
                            c(mIon[1],(1-sum(par[6:7]))*par[1],par[2],par[2],par[4],par[5],0,0),
                            c(mIon[2],par[6]*par[1],par[2],par[2],par[4],par[5],0,0),
                            c(mIon[3],par[7]*par[1],par[2],par[2],par[4],par[5],0,0)
                            )
  simuData <- rep(
                  simuACF(
                          ele = ele,
                          ion = ion,
                          kdir = c(1,0,0),
                          fradar = fradar,
                          scattAngle = aSite,
                          freq = freq,
                          lags = lags,
                          Bdir = c(1,0,0)
                          ),
                  nacf
                  ) + (rnorm(nData) + 1i*rnorm(nData))*dataStd/sqrt(2)
  simuVar   <- rep(dataStd**2,nData)


  plot(Re(simuACF(
                          ele = ele,
                          ion = ion,
                          kdir = c(1,0,0),
                          fradar = fradar,
                          scattAngle = aSite,
                          freq = freq,
                          lags = lags,
                          Bdir = c(1,0,0)
                          )))
  
  print(
        system.time(
                    fitpar   <- ISparamfit.guisdap(
                                           acf             = simuData,
                                           var             = simuVar,
                                           nData           = nData,
                                           fSite           = fradar,
                                           aSite           = aSite,
                                           initParam       = initParam,
                                           invAprioriCovar = apriori$invAprioriCovar,
                                           aprioriTheory   = apriori$aprioriTheory,
                                           aprioriMeas     = apriori$aprioriMeas,
                                           mIon            = c(30.5,16,1),
                                           paramLimits     = limitParam,
                                           directTheory    = ISdirectTheory.guisdap,
                                           absLimit        = absLimit,
                                           diffLimit       = diffLimit,
                                           scaleFun        = scaleParams,
                                           scale           = parScales,
                                           lags            = rep(lags,nacf),
                                           plotTest        = plotTest,
                                           maxLambda       = maxLambda,
                                           maxIter         = maxIter
                                           )
                    )
        )




  parfit         <- scaleParams(fitpar$param,scale=parScales,inverse=T)
  fitstd         <- scaleParams(sqrt(diag(fitpar$covar)),scale=parScales,inverse=T)

  for (k in seq(length(par))){
    cat(sprintf("%25.10f %25.10f %25.10f %25.10f \n",par[k],parfit[k],fitstd[k],parfit[k]/par[k]))
  }
  fitpar$parScales <- parScales
  invisible(fitpar)
} # testfit.guisdap
