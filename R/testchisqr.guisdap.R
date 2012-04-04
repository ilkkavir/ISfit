testchisqr.guisdap <- function(refPoint  = KIR,
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
                       heibeg    = 1,
                       heiend    = 1000,
                       heistp    = 1,
                       testParams= c(2),
                       testValues= list(seq(100,5000)),
                       initParam = NULL
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
#   testParams  indeces of the parameter to test for (1=Ne,2=Te,3=TeTi,4=coll,5=V,6=pm2,7=pm3)
#   testValues values to test for testParams
#   initParam  optional input plasma parameter vector, if NULL, parameters are taken from IRI-model at the given location and time
#
# I. Virtanen 2012  
#

  mIon <- c(30.5,16,1)

  aSite <- 180

  # parameters from iri model
  ptmp           <- iri(time=time,latitude=refPoint[1],longitude=refPoint[2],heibeg=heibeg,heiend=heiend,heistp=heistp)

  # height point closest to the user input value
  h              <- which(abs(seq(heibeg,heiend,by=heistp)-hTarg)==min(abs(seq(heibeg,heiend,by=heistp)-hTarg)))[1]

  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
  # the densities in outfmsis are in cm^-3
  ioncoll        <- (2.44e-16*ptmp$outfmsis[2,h] + 4.34e-16*ptmp$outfmsis[3,h] + 4.28e-16*ptmp$outfmsis[4,h])*1e6
  
  # parameters
  par            <- c( ptmp$outf[c(1,3),h] , ptmp$outf[4,h]/ptmp$outf[3,h] , ioncoll , vion , ptmp$outf[5,h]/100 , ptmp$outf[6,h]/100)

  # print iri parameter
  print(par)
  
  # initial (and apriori) parameter values
  if(is.null(initParam)) initParam      <- par

  # parameter scaling factors
  parScales      <- par
  parScales[3]   <- 1
  parScales[5:7] <- 1

  # scale the test values
  if(is.list(testValues)){
    testValuesScale <- testValues
  }else{
    testValuesScale <- list(testValues)
  }
  nTest           <- length(testParams)
  for(k in seq(nTest)){
    for(m in seq(length(testValuesScale[[k]]))){
      parTmp                   <- initParam
      parTmp[testParams[[k]]]  <- testValues[[k]][m]
      testValuesScale[[k]][m]  <- scaleParams( parTmp , parScales , inverse=F )[testParams[k]]
    }
  }

  # scale the initial parameter values
  initParam      <- scaleParams( initParam , parScales , inverse=F)


  # apriori information
  apriori        <- ISapriori.default.guisdap( scaleParams(par,parScales,inverse=F) , nIon )


  #
  # generate the simulated ACF data and other
  #
  nData    <- length(lags) * nacf
  freq     <- seq(-100000,100000,by=1000)*fradar/1e9
  simuData <- rep(
                  simuACF(
                          ele = c(par[1],par[2]*par[3],par[2]*par[3],par[4]*.35714,par[5],0,0),
                          ion = list(
                            c(mIon[1],(1-sum(par[6:7]))*par[1],par[2],par[2],par[4],par[5],0,0),
                            c(mIon[2],par[6]*par[1],par[2],par[2],par[4],par[5],0,0),
                            c(mIon[3],par[7]*par[1],par[2],par[2],par[4],par[5],0,0)
                            ),
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

  # frequency ambiguity function
  ntot <- nacf*length(lags)
  lags <- rep(lags,length.out=ntot)
  nf <- length(freq)
  nl <- length(lags)
  fAmb <- matrix(nrow=nf,ncol=nl)
  for(k in seq(nl)) fAmb[,k] <- frequencyAmbiguity( lags[k] , freq )

  chisqrs   <- chisqrTest(
                         measData        = simuData,
                         measVar         = simuVar,
                         nData           = length(measData),
                         fSite           = fradar,
                         aSite           = aSite,
                         initParam       = initParam,
                         invAprioriCovar = apriori$invAprioriCovar,
                         aprioriTheory   = apriori$aprioriTheory,
                         aprioriMeas     = apriori$aprioriMeas,
                         mIon            = c(30.5,16,1),
                         directTheory    = ISdirectTheory.guisdap,
                         scaleFun        = scaleParams,
                         scale           = parScales,
                         lags            = lags,
                         fAmb            = fAmb,
                         xSite           = freq,
                         testParams      = testParams,
                         testValues      = testValuesScale
                         )

  return(chisqrs)

} # testchisqr.guisdap
