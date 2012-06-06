testfit.3D.beata <- function(
                             dataDirs  = list( TRO='/Volumes/R1V2/EISCAT_3D_tests/3D_paramfit/beata_cp1_2.0u_FI@uhf/20101207_00/',
                                               KIR='/Volumes/R1V2/EISCAT_3D_tests/3D_paramfit/beata_cp1_1.0r_FI@kir/20101207_00/',
                                               SOD='/Volumes/R1V2/EISCAT_3D_tests/3D_paramfit/beata_cp1_1.0r_FI@sod/20101207_00/'),
                             absCalib  = FALSE,
                             TiIsotropic=FALSE,
                             integrationTime = 60,
                             plotTest  = F,
                             absLimit  = 5,
                             diffLimit = 1e-2,
                             maxLambda = 1e30,
                             maxIter   = 100
                       ){
#
# Test the 3D plasma parameter fit with tristatic EISCAT beata data
#
# INPUT:
#   absCalib   if TRUE, all sites are assumed to be absolutely calibrated, if FALSE, the ACF scalings of all other sites but the first one are assumed to be unknown and are fitted
#   TiIsotropic If true, an isotropic ion temperature is assumed, otherwise difference in between parallel and perpendicular temperatures is allowed
#   integrationTime  integration time of the measurements, in seconds
#   plotTest   logical, if T, the measurements and direct theory values are plotted at each iteration step
#   absLimit   upper limit for fit residual, when the residual is reduced below absLimit, the solver starts to look for end condition
#              based on diffLimit
#   diffLimit  upper limit for the fractional change of residual in one iteration step. The iteration is stopped when residual
#              is below absLimit and the fractional change of residual is below diffLimit
#   maxLambda  Maximum value for the lambda parameter in levenberg-marquardt iteration
#   maxIter    Maximum number of iterations (Hess matrix calculations) in the levenberg-marquardt iteration
#
#
#
# I. Virtanen 2012  
#

  # load the GUISDAPIO package, not included in DESCRIPTION file on purpose
  require(GUISDAPIO)


  # read GUISDAP initialisations
  init.T <- readGUISDAPinit(expname='beata',site='T')
  init.R <- readGUISDAPinit(expname='beata',site='R')


  # read the data from files
  data.T <- readLagProfiles.GUISDAP( dfiles=dataDirs$T , init.T ) 
  data.K <- readLagProfiles.GUISDAP( dfiles=dataDirs$K , init.R ) 
  data.S <- readLagProfiles.GUISDAP( dfiles=dataDirs$S , init.R ) 


  # assume that all directories begin from the same integration period...
  # THIS WILL NEED TO BE REPLACED WITH SOMETHING BETTER!!!
  nint.T <- round(integrationTime/data.T$parbl[7,1])
  nint.K <- round(integrationTime/data.K$parbl[7,1])
  nint.S <- round(integrationTime/data.S$parbl[7,1])

  # go through all integration periods
  iper <- 1
  repeat{
    # columns of the data matrix to include in this integration period
    columns.T <- seq( (nint.T*(iper-1)+1) , (nint.T*iper) )
    columns.K <- seq( (nint.K*(iper-1)+1) , (nint.K*iper) )
    columns.S <- seq( (nint.S*(iper-1)+1) , (nint.S*iper) )

    if(max(columns.T)>dim(data.T$lp)[2]) break
    if(max(columns.K)>dim(data.K$lp)[2]) break
    if(max(columns.S)>dim(data.S$lp)[2]) break

    # beam ponting directions at the end of the integration period
    azel.T   <- data.T$parbl[c(10,9),max(columns.T)]
    azel.K   <- data.K$parbl[c(10,9),max(columns.K)]
    azel.S   <- data.S$parbl[c(10,9),max(columns.S)]

    # calculate the beam intersection location and distances from each site to the intersection
    intersect.TK <- beamIntersect.location( TRO , KIR , azel.T , azel.K )
    intersect.TS <- beamIntersect.location( TRO , SOD , azel.T , azel.S )

    # average of the two intersections calculated above
    intersect.xyz <- ( intersect.TK$intersect + intersect.TS$intersect ) / 2

    # latitude, longitude, and height of the intersection points
    intersect.latlon <- cartesianToSpherical( intersect.xyz , degrees=TRUE , r0=ISgeometry:::EarthRadius() )

    # scattering wave vectors and scattering angles
    kTT <- ISgeometry:::rotateHorizontal.vector.cartesian( ISgeometry:::scatterPlaneNormal.cartesian(intersect.TK$pdir1$site,intersect.TK$pdir1$site,intersect.xyz),
                                                           intersect.xyz
                                                          )
    kTK <- ISgeometry:::rotateHorizontal.vector.cartesian( ISgeometry:::scatterPlaneNormal.cartesian(intersect.TK$pdir1$site,intersect.TK$pdir2$site,intersect.xyz),
                                                           intersect.xyz
                                                          )
    kTS <- ISgeometry:::rotateHorizontal.vector.cartesian( ISgeometry:::scatterPlaneNormal.cartesian(intersect.TS$pdir1$site,intersect.TS$pdir2$site,intersect.xyz),
                                                           intersect.xyz
                                                          )
    aTT <- 180
    aTK <- ISgeometry:::vectorAngle.cartesian( (intersect.xyz-intersect.TK$pdir1$site) , (intersect.TK$pdir2$site-intersect.xyz),degrees=TRUE)
    aTS <- ISgeometry:::vectorAngle.cartesian( (intersect.xyz-intersect.TS$pdir1$site) , (intersect.TS$pdir2$site-intersect.xyz),degrees=TRUE)
    
    # magnetic field direction
    Btmp           <- igrf(date=data.T$parbl[1:3,1],lat=intersect.latlon[1],lon=intersect.latlon[2],height=intersect.latlon[3],isv=0,itype=1)
    B              <- c(Btmp$x,-Btmp$y,-Btmp$z) # the model has y-axis to east and z-axis downwards

    # find ranges to the beam intersection at Tromso,
    # SEEMS TO BE ALWAYS 1000 AT OTHER SITES, IS THIS TRUE!!!???!??!?
    rTT <- intersect.TK$R[1]*2/.299792458
    rTK <- 1000 #sum(intersect.TK$R)/.299792458
    rTS <- 1000 ##sum(intersect.TS$R)/.299792458

    # range-gate limits for the analysis
    rlimsTT <- round( rTT + c(-25,25) )
    rlimsTK <- round( rTK + c(-25,25) )
    rlimsTS <- round( rTS + c(-25,25) )

    # select data rows that have lag value larger than zero and range is within the above limits
    indTT <- which( ( data.T$lag > 10 ) & (data.T$lag%%10 == 0) & (data.T$range >= rlimsTT[1]) & (data.T$range <= rlimsTT[2]))
    indTK <- which( ( data.K$lag > 10 ) & (data.K$lag%%10 == 0) & (data.K$range >= rlimsTK[1]) & (data.K$range <= rlimsTK[2]))
    indTS <- which( ( data.S$lag > 10 ) & (data.S$lag%%10 == 0) & (data.S$range >= rlimsTS[1]) & (data.S$range <= rlimsTS[2]))

    dTT <- data.T$lp[indTT,columns.T]
    dTK <- data.K$lp[indTK,columns.K]
    dTS <- data.S$lp[indTS,columns.S]


    # create acf, variance, and lag vectors for each site
    lagTT <- data.T$lag[indTT]
    lagTK <- data.K$lag[indTK]
    lagTS <- data.S$lag[indTS]

    lagT <- unique(lagTT)
    lagK <- unique(lagTK)
    lagS <- unique(lagTS)

    acfT <- rep(0+0i,length(lagT))
    varT <- rep(0,length(lagT))
    for(k in seq(length(lagT))){
      indT <- which(lagTT==lagT[k])
      ndT  <- length(c(dTT[indT,]))
      acfT[k] <- mean(dTT[indT,]) * (2*intersect.TK$R[1]**2)
      varT[k] <- ( var(c(Re(dTT[indT,]))) + var(c(Im(dTT[indT,])) ) )  / ndT * (2*intersect.TK$R[1]**2)**2
    }
    acfK <- rep(0+0i,length(lagK))
    varK <- rep(0,length(lagK))
    for(k in seq(length(lagK))){
      indK <- which(lagTK==lagK[k])
      ndK  <- length(c(dTK[indK,]))
      acfK[k]   <- mean(dTK[indK,]) * sum(intersect.TK$R**2)
      varK[k] <- ( var(c(Re(dTK[indK,]))) + var(c(Im(dTK[indK,])))) / ndK * sum(intersect.TK$R**2)**2
    }
    acfS <- rep(0+0i,length(lagS))
    varS <- rep(0,length(lagS))
    for(k in seq(length(lagS))){
      indS <- which(lagTS==lagS[k])
      ndS  <- length(c(dTS[indS,]))
      acfS[k]   <- mean(dTS[indS,]) * sum(intersect.TS$R**2)
      varS[k] <- ( var(c(Re(dTS[indS,]))) + var(c(Im(dTS[indS,]))) ) / ndS * sum(intersect.TS$R**2)**2
    }

    
#    dev.new()
#    layout(matrix(seq(3),ncol=1))
#    plot(lagT,Re(acfT))
#    points(lagT,Im(acfT),col='red')
#    lines(lagT,sqrt(varT),col='blue')
#    plot(lagK,Re(acfK))
#    points(lagK,Im(acfK),col='red')
#    lines(lagK,sqrt(varK),col='blue')
#    plot(lagS,Re(acfS))
#    points(lagS,Im(acfS),col='red')
#    lines(lagS,sqrt(varS),col='blue')
#    dev.new()


    # initial parameters from IRI model
    #  # parameters from iri model
    ptmp           <- iriParams(time=data.T$parbl[1:6,max(columns.T)],latitude=intersect.latlon[1],longitude=intersect.latlon[2],heights=intersect.latlon[3])#
    # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
    # the densities in outfmsis are in cm^-3
    ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,1])['NO+',] )
    # an approximation for electron-neutral collision frequency
    elecoll        <- ioncoll*.35714
    # electron parameters
    ele            <- c(ptmp[c('e-','Te','Te'),1],elecoll,c(0,0,0))
    # ion parameters
    ion             <- list(
                            c(30.5,(ptmp['O2+',1]+ptmp['NO+',1])/ele[1],ptmp[c('Ti','Ti'),1],ioncoll,c(0,0,0)),
                            c(16.0,ptmp['O+',1]/ele[1],ptmp[c('Ti','Ti'),1],ioncoll,c(0,0,0)),
                            c(1.0,ptmp['H+',1]/ele[1],ptmp[c('Ti','Ti'),1],ioncoll,c(0,0,0))
                            )
    # number of ions
    nIon           <- length(ion)
    # parameters in the form used by ISparamfit
    initParam <- ISparamList2Vec(ele,ion,rep(1,3))
    # parameter scaling factors
    parScales      <- ISparamScales.default(initParam,nIon)
    # scale the initial parameter values
    initParam      <- scaleParams( initParam , parScales , inverse=F)
    # parameter value limits
    parLimits      <- ISparamLimits.default(nIon,3)
    # scale the parameter limits
    limitParam     <- parLimits
    limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
    limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
    # apriori information
    apriori        <- ISapriori.default.3D( initParam , nIon , absCalib , TiIsotropic )

    # wave vectors, scattering angles, and site indices
    nData <- length(acfT) + length(acfK) + length(acfS)
    isite <- c(rep(1,length(acfT)),rep(2,length(acfK)),rep(3,length(acfS)))
    acf   <- c(acfT/1e25,acfK/1e26,acfS/1e26)
    var   <- c(varT/1e50,varK/1e52,varS/1e52)
    lags  <- c(lagT,lagK,lagS)*1e-6
    fSite <- rep(933e6,3)
    aSite <- c(aTT,aTK,aTS)
    kSite <- list(kTT,kTK,kTS)

  
  print(
        system.time(
                    fitpar   <- ISparamfit(
                                           acf             = acf,
                                           var             = var,
                                           nData           = nData,
                                           iSite           = isite,
                                           fSite           = fSite,
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
                                           lags            = lags,
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
  initParam <- scaleParams(initParam,scale=parScales,inverse=TRUE)

  for (k in seq(length(initParam))){
    cat(sprintf("%25.10f %25.10f %25.10f´\n",initParam[k],parfit[k],fitstd[k]))
  }


    save(fitpar,parfit,fitstd,file=paste('iper',iper,'.Rdata',sep=''))
    
    iper <- iper + 1
  }












  

#  # if only one location is given as a vector, convert it into a list
#  if(!is.list(locTrans)) locTrans <- list(locTrans)
#  if(!is.list(locRec))   locRec   <- list(locRec)
#
#  # list lengths
#  nTrans <- length(locTrans)
#  nRec   <- length(locRec)
#  nComb  <- nTrans*nRec
#
#  # convert all coordinates to a cartesian system with origin at the centre of Earth, x-axis pointing to zero-meridian
#  # at equator, y-axis to 90-degrees east at equator, and z-axis towards north pole
#
#  # transmitters
#  xyzTrans <- vector(mode='list',length=nTrans)
#  for(k in seq(nTrans)){
#    if(locxy){
#      latlonTrans <- ISgeometry:::planarToSpherical.geographic(x=locTrans[[k]][1],y=locTrans[[k]][2],refPoint=refPoint)
#    }else{
#      latlonTrans <- list(lat=locTrans[[k]][1],lon=locTrans[[k]][2])
#    }
#    xyzTrans[[k]] <- ISgeometry:::sphericalToCartesian(c(latlonTrans$lat,latlonTrans$lon))
#  }
#  # receivers
#  xyzRec   <- vector(mode='list',length=nRec)
#  for(k in seq(nRec)){
#    if(locxy){
#      latlonRec  <- ISgeometry:::planarToSpherical.geographic(x=locRec[[k]][1],y=locRec[[k]][2],refPoint=refPoint)
#    }else{
#      latlonRec  <- list(lat=locRec[[k]][1],lon=locRec[[k]][2])
#    }
#    xyzRec[[k]]  <- ISgeometry:::sphericalToCartesian(c(latlonRec$lat,latlonRec$lon))
#  }
#  # target (here we need also latitude and longitude, because they are used by the IRI-model)
#  if(locxy){
#    latlonTarg   <- ISgeometry:::planarToSpherical.geographic(x=locTarg[1],y=locTarg[2],refPoint=refPoint)
#  }else{
#    latlonTarg   <- list(lat=locTarg[1],lon=locTarg[2])
#  }
#  xyzTarg        <- ISgeometry:::sphericalToCartesian(c(latlonTarg$lat,latlonTarg$lon,ISgeometry:::EarthRadius()+locTarg[3]))
#  
#  # lists for site parameters
#  kSite <- vector(mode='list',length=nComb)
#  aSite <- vector(mode='numeric',length=nComb)
#  fSite <- rep(fradar,length.out=nComb)#
#
#  for(t in seq(nTrans)){
#    for(r in seq(nRec)){
#      # scattering wave vector, in a coordinate system with horizontal x- and y-axes
#      kSite[[(t-1)*nRec+r]] <- ISgeometry:::rotateHorizontal.vector.cartesian(
#                                 ISgeometry:::scatterPlaneNormal.cartesian(xyzTrans[[t]],xyzRec[[r]],xyzTarg),
#                                 xyzTarg
#                               )                                                                 
#      # scattering angle (angle between incident and scattered waves)
#      aSite[(t-1)*nRec+r]   <- ISgeometry:::vectorAngle.cartesian((xyzTarg-xyzTrans[[t]]),(xyzRec[[r]]-xyzTarg),degrees=T)
#    }
#  }##
#
#  # parameters from iri model
#  ptmp           <- iriParams(time=time,latitude=latlonTarg[['lat']],longitude=latlonTarg[['lon']],heights=heights)#
#
#  # height point closest to the user input value
#  h              <- which(abs(heights-locTarg[3])==min(abs(heights-locTarg[3])))[1]
#
#  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
#  # the densities in outfmsis are in cm^-3
#  ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,h])['NO+',] )
#  
#  # an approximation for electron-neutral collision frequency
#  elecoll        <- ioncoll*.35714
#
#  # asymmentric ion temperatures
#  tion           <- ptmp['Ti',h]*c(1,TperpTpar)
#
#  # electron parameters
#  ele            <- c(ptmp[c('e-','Te','Te'),h],elecoll,vele)
#
#  # ion parameters
#  ion             <- list(
#                         c(30.5,(ptmp['O2+',h]+ptmp['NO+',h])/ele[1],tion,ioncoll,vion),
#                         c(16.0,ptmp['O+',h]/ele[1],tion,ioncoll,vion),
#                         c(1.0,ptmp['H+',h]/ele[1],tion,ioncoll,vion)
#                         )
#
#  
#  nIon           <- length(ion)
#
#
#  # parameters in the form used by ISparamfit
#  par            <- ISparamList2Vec(ele,ion,rep(1,nComb))
#
#  # initial (and apriori) parameter values
#  # set initial velocity to zero
##  initParam                         <- par
##  eleinit <- ele*(1+.1*rnorm(7))
##  ioninit <- list( ion[[1]]*(1+c(0,.1*rnorm(7))) , ion[[2]]*(1+c(0,.1*rnorm(7))) , ion[[3]]*(1+c(0,.1*rnorm(7))))
#  eleinit <- ele
#  ioninit <- ion
#  initParam <- ISparamList2Vec(eleinit,ioninit,rep(1,nComb))
#  initParam[seq(5,(8*nIon+7),by=8)] <- 0
#  initParam[seq(6,(8*nIon+7),by=8)] <- 0
#  initParam[seq(7,(8*nIon+7),by=8)] <- 0
#  # initially isotropic temperatures and Te=Ti
#  initParam[c(10,11)]                     <- initParam[10]
#  for(k in seq(0,nIon)) initParam[c(10,11)+((k-1)*8)] <- initParam[c(10,11)]
#  
#
#  # parameter scaling factors
#  parScales      <- ISparamScales.default(initParam,nIon)
#
#  # scale the initial parameter values
#  initParam      <- scaleParams( initParam , parScales , inverse=F)
#
#  # parameter value limits
#  parLimits      <- ISparamLimits.default(nIon,nComb)
#
#  # scale the parameter limits
#  limitParam     <- parLimits
#  limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
#  limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
#  
#  # apriori information
#  apriori        <- ISapriori.default.3D( initParam , nIon , absCalib , TiIsotropic )
#
#
#  # magnetic field direction
#  Btmp           <- igrf(date=time[1:3],lat=latlonTarg[['lat']],lon=latlonTarg[['lon']],height=locTarg[3],isv=0,itype=1)
#  B              <- c(Btmp$x,-Btmp$y,-Btmp$z) # the model has y-axis to east and z-axis downwards
#
#
#  # time-lags
#  lagres <- fwhmRange*2/299792.458
#  lags   <- (seq(nlags)-.5)*lagres
#
#
#  # simulated ACF data
#  simudata <- ISmeas.simu( refPoint=refPoint , locTrans=locTrans , locRec=locRec , locTarg=locTarg , locxy=locxy , fwhmTrans=fwhmTrans ,
#                          fwhmRec=fwhmRec , fwhmRange=fwhmRange , resNS=resNS , resEW=resEW , resH=resH , Pt=Pt , Tnoise=Tnoise , fradar=fradar ,
#                          phArrTrans=phArrTrans , phArrRec=phArrRec , ele=ele , ion=ion , freq=seq(-100000,100000,by=1000)*fradar/1e9 ,
#                          lags=lags , integrationTime=integrationTime , dutyCycle=dutyCycle , time=time)
#
#
#  # wave vectors, scattering angles, and site indices
#  nData <- nlags * nComb
#  isite <- rep(seq(nComb),each=nlags)
#  acf   <- unlist(simudata$ACF)
#  var   <- unlist(simudata$var)
#  lags  <- rep(lags,nComb)
#
#  
#  print(
#        system.time(
#                    fitpar   <- ISparamfit(
#                                           acf             = acf,
#                                           var             = var,
#                                           nData           = nData,
#                                           iSite           = isite,
#                                           fSite           = fSite,
#                                           aSite           = aSite,
#                                           initParam       = initParam,
#                                           invAprioriCovar = apriori$invAprioriCovar,
#                                           aprioriTheory   = apriori$aprioriTheory,
#                                           aprioriMeas     = apriori$aprioriMeas,
#                                           nIon            = 3,
#                                           paramLimits     = limitParam,
#                                           directTheory    = ISdirectTheory,
#                                           absLimit        = absLimit,
#                                           diffLimit       = diffLimit,
#                                           scaleFun        = scaleParams,
#                                           scale           = parScales,
#                                           lags            = lags,
#                                           plotTest        = plotTest,
#                                           B               = B,
#                                           kSite           = kSite,
#                                           maxLambda       = maxLambda,
#                                           maxIter         = maxIter
#                                           )
#                    )
#        )
#
#
#
#
#  parfit         <- scaleParams(fitpar$param,scale=parScales,inverse=T)
#  fitstd         <- scaleParams(sqrt(diag(fitpar$covar)),scale=parScales,inverse=T)
#
#  for (k in seq(length(par))){
#    cat(sprintf("%25.10f %25.10f %25.10f %25.10f \n",par[k],parfit[k],fitstd[k],parfit[k]/par[k]))
#  }
#  fitpar$parScales <- parScales
#  invisible(fitpar)
} # testfit.3D.beata



pointingDirection.cartesian <- function( locAnt=TRO , azel=c(185.79,77.39) ){

  # conversion from azimuth-elevation beam pointing and lattitude-longitude site location to cartesian coordinates

  # cartesian  coordinates of the site location
  xyzS <- ISgeometry:::sphericalToCartesian( locAnt , degrees=TRUE )

  # pointing direction in a local cartesian system, where x-axis points towards north, y-axis towards west, and z-axis to local zenith
  # azimuth multiplied by -1, because the function originally written for latitute-longitude conversions uses eastern longitudes, but the
  # EISCAT azimuths are positive clockwise
  xyzBloc <- ISgeometry:::sphericalToCartesian( c(azel[2],-azel[1],1) , degrees=TRUE)

  # find the axes of the local system
  #  y-axis is perpendicular to both the site position vector and c(0,0,1)
  ydir <- ISgeometry:::normalUnitVector.cartesian(xyzS,c(0,0,1))
  # x-axis is perpendicular to the x-axis and the position vector
  xdir <- ISgeometry:::normalUnitVector.cartesian(ydir,xyzS)
  # z-axis is parallel with the position vector
  zdir <- xyzS / sqrt(sum(xyzS**2))

  # conversion of the beam pointing to the cartesian system with origin at the centre of Earth, x-axis pointing towards equator at zero meridian, etc.
  xyzB <- xyzBloc[1] * xdir + xyzBloc[2] * ydir + xyzBloc[3] * zdir

  # return the site location and the beam pointing in cartesian coordinates
  return(list(site=xyzS,beam=xyzB))

} # pointingDirection.cartesian

beamIntersect.location <- function( site1 , site2 , azel1 , azel2 )
  {
    # estimate of the beam intersection point
    
    # pointing directions in cartesian coordinates
    pdir1 <- pointingDirection.cartesian( site1 , azel1 )
    pdir2 <- pointingDirection.cartesian( site2 , azel2 )

    # form a linear inverse problem, whose solution gives distances from the sites
    m <- pdir1$site - pdir2$site
    A <- matrix( c(-pdir1$beam,pdir2$beam),ncol=2,byrow=FALSE)

    R <- qr.solve(a=A,b=m)

    intersect <- pdir1$site + R[1]*pdir1$beam

    return(list(intersect=intersect,pdir1=pdir1,pdir2=pdir2,R=R))

  }

cartesianToSpherical <- function( xyz , degrees=TRUE , r0=0){

  r <- sqrt(sum(xyz**2)) - r0

  az <- atan2(x=xyz[1],y=xyz[2])
  el <- atan2(x=sqrt(sum(xyz[1:2]**2)),y=xyz[3])

  res <- c(el,az,r)
  
  if(degrees) res[1:2] <- res[1:2]*180/pi


  return(res)
  
}
