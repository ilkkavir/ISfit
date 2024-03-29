testfit <- function(
    general=FALSE,
    refPoint  = KIR,
    locTrans  = list(c(0,0)),
    locRec    = list( c(    0 ,    0 ),
        c(  100 ,    0 ),
        c( -100 ,    0 ),
        c(    0 ,  100 ),
        c(    0 , -100 )
        ),
    locTarg   = c(0,0,200),
    locxy     = T,
    refSite = 1,

    fwhmTrans = c(1),
    fwhmRec   = c(1),
    fwhmRange = c(1),
    fwhmIonSlab = 100,
    resNS,
    resEW,
    resH,
    Pt        = 1e6,
    Tnoise    = 200,
    fradar    = 235e6,
    phArrTrans= FALSE,
    phArrRec  = FALSE,
    absCalib  = FALSE,
    TiIsotropic=FALSE,
    TeIsotropic=TRUE,

    TperpTpar = 1,
    vion      = c(0,0,0),
    vele      = vion,

    nlags     = 30,
    integrationTime = 10,
    dutyCycle = .1,

    plotFit  = T,
    plotTest  = F,
    absLimit  = 5,
    diffLimit = 1e-2,
    maxLambda = 1e30,
    maxIter   = 100,

    time      = c(2009,7,1,11,0,0),
    heights   = seq(1000),
    fitFun=leastSquare.lvmrq,
    ...
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
#   fwhmTrans  full width at half maximum of bore-sight transmitter beam(s), in degrees
#   fwhmRec    full width at half maximum of bore-sight receiver beam(s), in degrees
#   fwhmRange  full width at half maximum of the gaussian range ambiguity function, in km
#   resNS      optional, north-south resolution of spatial integration voxels, in km
#   resEW      optional, east-west resolution of spatial integration voxels, in km
#   resH       optinal, height resolution of spatial integration voxels, in km
#   fradar     transmitter frequencies [Hz]
#   phArrTrans TRUE if transmitter antenna(s) is(are) phased arrays, a vector with own entry for each antenna
#   phArrRec   see above, for receiver antennas this time
#   absCalib   if TRUE, all sites are assumed to be absolutely calibrated, if FALSE, the ACF scalings of all other sites but the first one are assumed to be unknown and are fitted
#   TiIsotropic If true, an isotropic ion temperature is assumed, otherwise difference in between parallel and perpendicular temperatures is allowed
#
#   TperpTpar  ratio of perpendicular and parallel ion temperatures (with respect to local magnetic field direction)
#   vion       ion velocity vector in cartesian coordinates (z-axis upwards)
#   vele       electron velocity
#
#   cSite      ACF scaling factors for each bistatic combination (scaling with distance and beam shapes are not included in this test)
#   nlags      number of time-lags, the lag resolution is matched with fwhmRange
#   integrationTime  integration time of the measurements, in seconds
#   dutyCycle  duty cycle of the transmitter(s)
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
#   heights    heights at which iri parameters are originally being calculated, most probably the user does not need to change this value
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
      latlonTrans <- ISgeometry:::planarToSpherical.geographic(x=locTrans[[k]][1],y=locTrans[[k]][2],refPoint=refPoint)
    }else{
      latlonTrans <- list(lat=locTrans[[k]][1],lon=locTrans[[k]][2])
    }
    xyzTrans[[k]] <- ISgeometry:::sphericalToCartesian(c(latlonTrans$lat,latlonTrans$lon))
  }
  # receivers
  xyzRec   <- vector(mode='list',length=nRec)
  for(k in seq(nRec)){
    if(locxy){
      latlonRec  <- ISgeometry:::planarToSpherical.geographic(x=locRec[[k]][1],y=locRec[[k]][2],refPoint=refPoint)
    }else{
      latlonRec  <- list(lat=locRec[[k]][1],lon=locRec[[k]][2])
    }
    xyzRec[[k]]  <- ISgeometry:::sphericalToCartesian(c(latlonRec$lat,latlonRec$lon))
  }
  # target (here we need also latitude and longitude, because they are used by the IRI-model)
  if(locxy){
    latlonTarg   <- ISgeometry:::planarToSpherical.geographic(x=locTarg[1],y=locTarg[2],refPoint=refPoint)
  }else{
    latlonTarg   <- list(lat=locTarg[1],lon=locTarg[2])
  }
  xyzTarg        <- ISgeometry:::sphericalToCartesian(c(latlonTarg$lat,latlonTarg$lon,ISgeometry:::EarthRadius()+locTarg[3]))

  # lists for site parameters
  kSite <- vector(mode='list',length=nComb)
  aSite <- vector(mode='numeric',length=nComb)
  fSite <- rep(fradar,length.out=nComb)

  for(t in seq(nTrans)){
    for(r in seq(nRec)){
      # scattering wave vector, in a coordinate system with horizontal x- and y-axes
      kSite[[(t-1)*nRec+r]] <- ISgeometry:::rotateHorizontal.vector.cartesian(
                                 ISgeometry:::scatterPlaneNormal.cartesian(xyzTrans[[t]],xyzRec[[r]],xyzTarg),
                                 xyzTarg
                               )
      # scattering angle (angle between incident and scattered waves)
      aSite[(t-1)*nRec+r]   <- ISgeometry:::vectorAngle.cartesian((xyzTarg-xyzTrans[[t]]),(xyzRec[[r]]-xyzTarg),degrees=T)
    }
  }

  # parameters from iri model
  ptmp           <- iriParams(time=time,latitude=latlonTarg[['lat']],longitude=latlonTarg[['lon']],heights=heights)

  # height point closest to the user input value
  h              <- which(abs(heights-locTarg[3])==min(abs(heights-locTarg[3])))[1]

  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
  # the densities in outfmsis are in cm^-3
  ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,h])['NO+',] )

  # an approximation for electron-neutral collision frequency
  elecoll        <- ioncoll*.35714

  # asymmentric ion temperatures
  tion           <- ptmp['Ti',h]*c(1,TperpTpar)

  # electron parameters
  ele            <- c(ptmp[c('e-','Te','Te'),h],elecoll,vele)

  # ion parameters
  ion             <- list(
                         c(30.5,(ptmp['O2+',h]+ptmp['NO+',h])/ele[1],tion,ioncoll,vion),
                         c(16.0,ptmp['O+',h]/ele[1],tion,ioncoll,vion),
                         c(1.0,ptmp['H+',h]/ele[1],tion,ioncoll,vion)
                         )


  nIon           <- length(ion)
  mIon <- sapply(ion,FUN=function(x){x[1]})


  # parameters in the form used by ISparamfit
  if(general){
      par            <- ISparamList2Vec.general(ele,ion,rep(1,nComb))
      eleinit <- ele
      ioninit <- ion
      initParam <- ISparamList2Vec.general(eleinit,ioninit,rep(1,nComb))
      initParam[seq(5,(8*nIon+7),by=8)] <- 0
      initParam[seq(6,(8*nIon+7),by=8)] <- 0
      initParam[seq(7,(8*nIon+7),by=8)] <- 0
      # initially isotropic temperatures and Te=Ti
      initParam[c(10,11)]                     <- initParam[10]
      for(k in seq(0,nIon)) initParam[c(10,11)+((k-1)*8)] <- initParam[c(10,11)]

      # parameter scaling factors
      parScales      <- ISparamScales.general(initParam,nIon)

      # scale the initial parameter values
      initParam      <- scaleParams( initParam , parScales , inverse=F)

      # parameter value limits
      parLimits      <- ISparamLimits.general(nIon,nComb)

      # scale the parameter limits
      limitParam     <- parLimits
      limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
      limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)

      # apriori information
      apriori <- ISapriori.general( initParam , nIon , absCalib , TiIsotropic )
      directTheory <- ISdirectTheory.general

  }else{

      par            <- ISparamList2Vec(ele,ion,rep(1,nComb))
      eleinit <- ele
      ioninit <- ion
      initParam <- ISparamList2Vec(eleinit,ioninit,rep(1,nComb))
      initParam[7:9] <- 0

      # initially isotropic temperatures and Te=Ti
      initParam[3:5] <- initParam[2]

      # parameter scaling factors
      parScales      <- ISparamScales(initParam,nIon)

      # scale the initial parameter values
      initParam      <- scaleParams( initParam , parScales , inverse=F)

      # parameter value limits
      parLimits      <- ISparamLimits(nIon,nComb)

      # scale the parameter limits
      limitParam     <- parLimits
      limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
      limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)

      # apriori information
      apriori <- ISapriori( aprioriParam=initParam , nIon=nIon , absCalib=absCalib , TiIsotropic=TiIsotropic , TeIsotropic=TeIsotropic , refSite=refSite )

      directTheory <- ISdirectTheory
  }

  # magnetic field direction
  Btmp <- igrf(date=time[1:3],lat=latlonTarg[['lat']],lon=latlonTarg[['lon']],height=locTarg[3],isv=0,itype=1)
  B <- c(Btmp$x,-Btmp$y,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards north,
                                 # y towards west and z upwards


  # time-lags
  lagres <- fwhmRange*2/299792.458
  lags   <- (seq(nlags)-.5)*lagres


  # simulated ACF data
  simudata <- ISmeas.simu( refPoint=refPoint , locTrans=locTrans , locRec=locRec , locTarg=locTarg , locxy=locxy , fwhmTrans=fwhmTrans , fwhmRec=fwhmRec , fwhmRange=fwhmRange , fwhmIonSlab=fwhmIonSlab , resNS=resNS , resEW=resEW , resH=resH , Pt=Pt , Tnoise=Tnoise , fradar=fradar , phArrTrans=phArrTrans , phArrRec=phArrRec , ele=ele , ion=ion , freq=seq(-100000,100000,by=1000)*fradar/1e9 , lags=lags , integrationTime=integrationTime , dutyCycle=dutyCycle , time=time)

  # wave vectors, scattering angles, and site indices
  nData <- nlags * nComb
  isite <- rep(seq(nComb),each=nlags)
  acf   <- unlist(simudata$ACF)
  var   <- unlist(simudata$var)
  lags  <- rep(lags,nComb)

  print(
      system.time(
          fitpar   <- ISparamfit(
              acf             = acf,
              var             = var,
              lags            = lags,
              iSite           = isite,
              fSite           = fSite,
              aSite           = aSite,
              kSite           = kSite,
              B               = B,
              initParam       = initParam,
              aprioriTheory   = apriori$aprioriTheory,
              aprioriMeas     = apriori$aprioriMeas,
              invAprioriCovar = apriori$invAprioriCovar,
              nIon            = 3,
              paramLimits     = limitParam,
              directTheory    = directTheory,
              nData           = nData,
              absLimit        = absLimit,
              diffLimit       = diffLimit,
              scaleFun        = scaleParams,
              scale           = parScales,
              plotTest        = plotTest,
              plotFit        = plotFit,
              maxLambda       = maxLambda,
              maxIter         = maxIter,
              mIon=mIon,
              fitFun=fitFun,
              ...
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
}
