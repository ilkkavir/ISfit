testfit.3D <- function(refPoint  = KIR,
                       locTrans  = list(c(0,0)),
                       locRec    = list( c(    0 ,    0 ),
                                         c(  100 ,    0 ),
                                         c( -100 ,    0 ),
                                         c(    0 ,  100 ),
                                         c(    0 , -100 )
                         ),
                       locTarg   = c(0,0,200),
                       locxy     = T,

                       TperpTpar = 1,
                       vion      = c(0,0,0),
                       vele      = vion,
                       
                       cSite     = rep(1,length(locRec)*length(locTrans)),
                       fradar    = rep(235e6,length(locTrans)),
                       nacf      = rep(50,length(locRec)*length(locTrans)),
                       lags      = rep(list(seq(50)*5e-5),length(locRec)*length(locTrans)),
                       dataStd   = rep(1e-19,length(locRec)*length(locTrans)),

                       plotTest  = F,
                       absLimit  = 5,
                       diffLimit = 1e-2,
                       maxLambda = 1e30,
                       maxIter   = 100,

                       time      = c(2009,7,1,11,0,0),
                       heibeg    = 1,
                       heiend    = 1000,
                       heistp    = 1
                       ){
#
# Test the 3D plasma parameter fit with simulated ACF data
#
# INPUT:
#   refPoint   c(lat,lon) reference point, from which the distance x and y are measured.
#   locTrans   list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the transmitter antennas.
#              Angles in degrees.
#              OR list(c(x,y)) position of the transmitter, IF locxy=T
#   locRec     list(c(lat,lon)) or list(c(lat,lon,height)) latitudes, longitudes (and heights) of the receiver antennas.
#              Angles in degrees.
#              OR list(c(x,y)) IF locxy=T
#   locTarg    c(x,y,h) if locxy=T, or c(lat,lon,h) if locxy=F, target plasma volume location
#   locxy      logical, if T, the positions of the transmitter and receiver are given in cartesian coordinates
#
#   TperpTpar  ratio of perpendicular and parallel ion temperatures (with respect to local magnetic field direction)
#   vion       ion velocity vector in cartesian coordinates (z-axis upwards)
#   vele       electron velocity
#
#   cSite      ACF scaling factors for each bistatic combination (scaling with distance and beam shapes are not included in this test)
#   fradar     transmitter frequencies [Hz]
#   nacf       number of autocorrelation functions (integration periods) from each bistatic combination
#   lags       time-lags measured with each bistatic combination [us]
#   dataStd    standard deviations of measurement noise for each bistatic combination
#
#              - the lists for bistatic combinations are given with receiver index running faster. e.g. with N transmitters and
#                M receivers the order would be T1R1, T1R2, ... , T1RM, T2R1, ... , T2RM, T3R1, ... , TNRM
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
#
#
#
# I. Virtanen 2012  
#

  # if only one location is given as a vector, convert it into a list
  if(!is.list(locTrans)) locTrans <- list(locTrans)
  if(!is.list(locRec))   locRec   <- list(locRec)

  # list lengths
  nTrans <- length(locTrans)
  nRec   <- length(locRec)
  nComb  <- nTrans*nRec

  # convert all coordinates to a cartesian system with origin at the centre of Earth, x-axis pointing to zero-meridian
  # at equator, y-axis to 90-degrees east at equator, and z-axis towards north pole

  # transmitters
  xyzTrans <- vector(mode='list',length=nTrans)
  for(k in seq(nTrans)){
    if(locxy){
      latlonTrans <- planarToSpherical.geographic(x=locTrans[[k]][1],y=locTrans[[k]][2],refPoint=refPoint)
    }else{
      latlonTrans <- locTrans[[k]]
    }
    xyzTrans[[k]] <- sphericalToCartesian(c(latlonTrans$lat,latlonTrans$lon))
  }
  # receivers
  xyzRec   <- vector(mode='list',length=nRec)
  for(k in seq(nRec)){
    if(locxy){
      latlonRec  <- planarToSpherical.geographic(x=locRec[[k]][1],y=locRec[[k]][2],refPoint=refPoint)
    }else{
      latlonRec  <- locRec[[k]]
    }
    xyzRec[[k]]  <- sphericalToCartesian(c(latlonRec$lat,latlonRec$lon))
  }
  # target (here we need also latitude and longitude, because they are used by the IRI-model)
  if(locxy){
    latlonTarg   <- planarToSpherical.geographic(x=locTarg[1],y=locTarg[2],refPoint=refPoint)
  }else{
    latlonTarg   <- locTarg[1:2]
  }
  xyzTarg        <- sphericalToCartesian(c(latlonTarg$lat,latlonTarg$lon,EarthRadius()+locTarg[3]))
  
  # lists for site parameters
  kSite <- vector(mode='list',length=nComb)
  aSite <- vector(mode='numeric',length=nComb)

  for(t in seq(nTrans)){
    for(r in seq(nRec)){
      # scattering wave vector, in a coordinate system with horizontal x- and y-axes
      kSite[[(t-1)*nRec+r]] <- rotateHorizontal.vector.cartesian(
                                 scatterPlaneNormal.cartesian(xyzTrans[[t]],xyzRec[[r]],xyzTarg),
                                 xyzTarg
                               )                                                                 
      # scattering angle (angle between incident and scattered waves)
      aSite[(t-1)*nRec+r]   <- vectorAngle.cartesian((xyzTarg-xyzTrans[[t]]),(xyzRec[[r]]-xyzTarg),degrees=T)
    }
  }

  # parameters from iri model
  ptmp           <- iri(time=time,latitude=latlonTarg[['lat']],longitude=latlonTarg[['lon']],heibeg=heibeg,heiend=heiend,heistp=heistp)

  # height point closest to the user input value
  h              <- which(abs(seq(heibeg,heiend,by=heistp)-locTarg[3])==min(abs(seq(heibeg,heiend,by=heistp)-locTarg[3])))[1]

  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
  # the densities in outfmsis are in cm^-3
  ioncoll        <- (2.44e-16*ptmp$outfmsis[2,h] + 4.34e-16*ptmp$outfmsis[3,h] + 4.28e-16*ptmp$outfmsis[4,h])*1e6
  # an approximation for electron-neutral collision frequency
  elecoll        <- ioncoll*.35714

  # asymmentric ion temperatures
  tion           <- ptmp$outf[3,h]*c(1,TperpTpar)

  # electron parameters
  ele            <- c(ptmp$outf[c(1,4,4),h],elecoll,vele)

  # ion parameters
  ion             <- list(
                         c(30.5,(ptmp$outf[8,h]+ptmp$outf[9,h])*ele[1]/100,tion,ioncoll,vion),
                         c(16.0,ptmp$outf[5,h]*ele[1]/100,tion,ioncoll,vion),
                         c(1.0,ptmp$outf[6,h]*ele[1]/100 ,tion,ioncoll,vion)
                         )

  nIon           <- length(ion)


  # parameters in the form used by ISparamfit
  par            <- ISparamList2Vec(ele,ion,cSite)

  # initial (and apriori) parameter values
  # set initial velocity to zero
  initParam                         <- par
  initParam[seq(5,(8*nIon+7),by=8)] <- 0
  initParam[seq(6,(8*nIon+7),by=8)] <- 0
  initParam[seq(7,(8*nIon+7),by=8)] <- 0
  # initially isotropic temperatures and Te=Ti
  initParam[c(10,11)]                     <- initParam[10]
  for(k in seq(0,nIon)) initParam[c(10,11)+((k-1)*8)] <- initParam[c(10,11)]

  # parameter scaling factors
  parScales      <- ISparamScales.default(initParam,nIon)

  # scale the initial parameter values
  initParam      <- scaleParams( initParam , parScales , inverse=F)

  # parameter value limits
  parLimits      <- ISparamLimits.default(nIon,nComb)

  # scale the parameter limits
  limitParam     <- parLimits
  limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
  limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
  
  # apriori information
  apriori        <- ISapriori.default.3D( initParam , nIon )


  # magnetic field direction
  Btmp           <- igrf(date=time[1:3],lat=latlonTarg[['lat']],lon=latlonTarg[['lon']],height=locTarg[3],isv=0,itype=1)
  B              <- c(Btmp$x,-Btmp$y,-Btmp$z) # the model has y-axis to east and z-axis downwards
  
  #
  # generate the simulated ACF data and other
  #
  nData    <- sum( sapply(lags,length) * nacf )
  simuData <- vector(mode='complex',length=nData)
  simuVar  <- vector(mode='numeric',length=nData)
  iSite    <- vector(mode='integer',length=nData)
  n        <- 0
  for(t in seq(nTrans)){
    for(r in seq(nRec)){
      k                      <- (t-1)*nRec+r
      nd                     <- length(lags[[k]]) * nacf[k]
      freq                   <- seq(-100000,100000,by=1000)*fradar[t]/1e9
      simuData[(n+1):(n+nd)] <- cSite[k] *
                                  rep(
                                      simuACF(
                                              ele        = ele,
                                              ion        = ion,
                                              kdir       = kSite[[k]],
                                              fradar     = fradar[t],
                                              scattAngle = aSite[k],
                                              freq       = freq,
                                              lags       = lags[[k]],
                                              Bdir       = B
                                              ),
                                      nacf[k]
                                      ) + (rnorm(nd) + 1i*rnorm(nd))*dataStd[k]/sqrt(2)
      simuVar[(n+1):(n+nd)]  <- dataStd[k]**2
      iSite[(n+1):(n+nd)]    <- k
      n                      <- n + nd
    }
  }
  print(
        system.time(
                    fitpar   <- ISparamfit(
                                           acf             = simuData,
                                           var             = simuVar,
                                           nData           = nData,
                                           iSite           = iSite,
                                           fSite           = fradar,
                                           aSite           = aSite,
                                           initParam       = initParam,
                                           invAprioriCovar = apriori$invAprioriCovar,
                                           aprioriTheory   = apriori$aprioriTheory,
                                           aprioriMeas     = apriori$aprioriMeas,
                                           nIon            = 3,
                                           paramLimits     = limitParam,
                                           directTheory    = ISdirectTheory,
                                           absLimit        = absLimit,
                                           diffLimit       = diffLimit,
                                           scaleFun        = scaleParams,
                                           scale           = parScales,
                                           lags            = unlist(lags),
                                           plotTest        = plotTest,
                                           B               = B,
                                           kSite           = kSite,
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
} # testfit.3D
