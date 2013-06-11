testfit.3D.beata <- function(dataFile = NULL,
                             dataDirs  = list( TRO='/media/raid/EISCATdata/2010/beata_cp1_2.0u_FI@uhf/20101207_00/',
                                               KIR='/media/raid/EISCATdata/2010/beata_cp1_1.0r_FI@kir/20101207_00/',
                                               SOD='/media/raid/EISCATdata/2010/beata_cp1_1.0r_FI@sod/20101207_00/'),
                             absCalib  = FALSE,
                             TiIsotropic=FALSE,
                             integrationTime = 60,
                             plotTest  = F,
                             plotFit   = T,
                             absLimit  = 5,
                             diffLimit = 1e-2,
                             maxLambda = 1e30,
                             maxIter   = 100,
                             aprioriFunction=ISapriori.3D,
                             sites=c('TRO','KIR','SOD'),
                             fitFun=leastSquare.lvmrq,
                             ...
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
  # Convert sites to lower-case one-character form
  sites <- substr(tolower(sites),1,1)
  tro <- any(sites=='t')
  kir <- any(sites=='k')
  sod <- any(sites=='s')


  # load the GUISDAPIO package, not included in DESCRIPTION file on purpose
  if(is.null(dataFile)){
      require(GUISDAPIO)


      # read GUISDAP initialisations
      init.T <- readGUISDAPinit(expname='beata',site='T')
      init.R <- readGUISDAPinit(expname='beata',site='R')


      # read the data from files
      data.T <- readLagProfiles.GUISDAP( dfiles=dataDirs$T , init.T ) 
      data.K <- readLagProfiles.GUISDAP( dfiles=dataDirs$K , init.R ) 
      data.S <- readLagProfiles.GUISDAP( dfiles=dataDirs$S , init.R ) 

  }else{
      load(dataFile)
  }

  # assume that data is continuous in all directories
  # THIS WILL NEED TO BE REPLACED WITH SOMETHING BETTER!!!
  nint.T <- round(integrationTime/data.T$parbl[7,1])
  nint.K <- round(integrationTime/data.K$parbl[7,1])
  nint.S <- round(integrationTime/data.S$parbl[7,1])

  btime.T <- data.T$parbl[11,1]
  btime.K <- data.K$parbl[11,1]
  btime.S <- data.S$parbl[11,1]
  btime <- min(c(btime.T,btime.K,btime.S))

  first.T <- min(which(data.T$parbl[11,]>=btime)) - 1
  first.K <- min(which(data.K$parbl[11,]>=btime)) - 1
  first.S <- min(which(data.S$parbl[11,]>=btime)) - 1

  # go through all integration periods
  iper <- 1
  repeat{
    # columns of the data matrix to include in this integration period
    columns.T <- seq( (nint.T*(iper-1)+1) , (nint.T*iper) ) + first.T
    columns.K <- seq( (nint.K*(iper-1)+1) , (nint.K*iper) ) + first.K
    columns.S <- seq( (nint.S*(iper-1)+1) , (nint.S*iper) ) + first.S

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
    if(azel.K[2]<0){
        intersect.xyz <- intersect.TS$intersect
    }else if(azel.S[2]<0){
        intersect.xyz <- intersect.TK$intersect
    }else{
        intersect.xyz <- ( intersect.TK$intersect + intersect.TS$intersect ) / 2
    }
    
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
    B              <- c(Btmp$x,-Btmp$y,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards north,
                                                # y towards west and z upwards


    # the "gain integral" assuming Gaussian beam-shapes + the factor from wave length + the factor from polarisation
    gainIntT <- ISgeometry:::bistaticResolutions.planar(refPoint=TRO,locTrans=TRO,locRec=TRO,locxy=FALSE,fwhmTrans=.6,fwhmRec=.6,fwhmRange=.001,x=0,y=0,height=intersect.latlon[3],phArrTrans=FALSE,phArrRec=FALSE)[["gainInt"]][1,1,1] * .32**2/4/pi * .5 * (1+cos(2*aTT*pi/180)**2)
    gainIntK <- ISgeometry:::bistaticResolutions.planar(refPoint=TRO,locTrans=TRO,locRec=KIR,locxy=FALSE,fwhmTrans=.6,fwhmRec=.6,fwhmRange=.001,x=0,y=0,height=intersect.latlon[3],phArrTrans=FALSE,phArrRec=FALSE)[["gainInt"]][1,1,1] * .32**2/4/pi * .5 * (1+cos(2*aTK*pi/180)**2)
    gainIntS <- ISgeometry:::bistaticResolutions.planar(refPoint=TRO,locTrans=TRO,locRec=SOD,locxy=FALSE,fwhmTrans=.6,fwhmRec=.6,fwhmRange=.001,x=0,y=0,height=intersect.latlon[3],phArrTrans=FALSE,phArrRec=FALSE)[["gainInt"]][1,1,1] * .32**2/4/pi * .5 * (1+cos(2*aTS*pi/180)**2)

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

    # read transmitter power from Tromso data
    txpow <- mean(data.T$parbl[8,columns.T])

    
    acfT <- rep(0+0i,length(lagT))
    varT <- rep(0,length(lagT))
    for(k in seq(length(lagT))){
      indT <- which(lagTT==lagT[k])
      ndT  <- length(c(dTT[indT,]))
      acfT[k] <- mean(dTT[indT,],na.rm=TRUE) / gainIntT / txpow
      varT[k] <- ( var(c(Re(dTT[indT,])),na.rm=TRUE) + var(c(Im(dTT[indT,])),na.rm=TRUE ) )  / ndT / gainIntT**2 / txpow**2
    }
    acfK <- rep(0+0i,length(lagK))
    varK <- rep(0,length(lagK))
    for(k in seq(length(lagK))){
      indK <- which(lagTK==lagK[k])
      ndK  <- length(c(dTK[indK,]))
      acfK[k]   <- mean(dTK[indK,],na.rm=TRUE) / gainIntK / txpow
      varK[k] <- ( var(c(Re(dTK[indK,])),na.rm=TRUE) + var(c(Im(dTK[indK,])),na.rm=TRUE)) / ndK / gainIntK**2 / txpow**2
    }
    acfS <- rep(0+0i,length(lagS))
    varS <- rep(0,length(lagS))
    for(k in seq(length(lagS))){
      indS <- which(lagTS==lagS[k])
      ndS  <- length(c(dTS[indS,]))
      acfS[k]   <- mean(dTS[indS,],na.rm=TRUE) / gainIntS / txpow
      varS[k] <- ( var(c(Re(dTS[indS,])),na.rm=TRUE) + var(c(Im(dTS[indS,])),na.rm=TRUE) ) / ndS / gainIntS**2 / txpow**2
    }

    

    # initial parameters from IRI model
    #  # parameters from iri model
    ptmp           <- iriParams(time=data.T$parbl[1:6,max(columns.T)],latitude=intersect.latlon[1],longitude=intersect.latlon[2],heights=intersect.latlon[3])
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
    # ion masses
    mIon <- sapply(ion,FUN=function(x){x[1]})
    # parameters in the form used by ISparamfit
    initParam <- ISparamList2Vec(ele,ion,rep(1,3))
    # parameter scaling factors
    parScales      <- ISparamScales(initParam,nIon)
    # scale the initial parameter values
    initParam      <- scaleParams( initParam , parScales , inverse=F)
    # parameter value limits
    parLimits      <- ISparamLimits(nIon,3)
    # scale the parameter limits
    limitParam     <- parLimits
    limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
    limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
    # apriori information
    apriori        <- aprioriFunction( initParam , nIon , absCalib , TiIsotropic )

    # cut off sites if they should not be used, stupid but works...
    if(!tro) acfT<-acfT[]*NA
    if(!kir) acfK<-acfK[]*NA
    if(!sod) acfS<-acfS[]*NA
    
    # wave vectors, scattering angles, and site indices
    nData <- length(acfT) + length(acfK) + length(acfS)
    isite <- c(rep(1,length(acfT)),rep(2,length(acfK)),rep(3,length(acfS)))
    acf   <- c(acfT,acfK,acfS)
    var   <- c(varT,varK,varS)
    lags  <- c(lagT,lagK,lagS)*1e-6
    fSite <- rep(933e6,3)
    aSite <- c(aTT,aTK,aTS)
    kSite <- list(kTT,kTK,kTS)
    inds <- which(is.na(acf)|is.na(var))
    acf[inds] <- 0+0i
    var[inds] <- Inf
    if(txpow>0){  
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
                    plotFit         = plotFit,
                    B               = B,
                    kSite           = kSite,
                    maxLambda       = maxLambda,
                    maxIter         = maxIter,
                    mIon=mIon,
                    fitFun=fitFun,
                    ...
                    )
                )
            )
        
        
        
        
#        parfit         <- scaleParams(fitpar$param,scale=parScales,inverse=T)
#        fitstd         <- scaleParams(sqrt(diag(fitpar$covar)),scale=parScales,inverse=T)
#        initParam <- scaleParams(initParam,scale=parScales,inverse=TRUE)
        
#        for (k in seq(length(initParam))){
#            cat(sprintf("%25.10f %25.10f %25.10f \n",initParam[k],parfit[k],fitstd[k]))
#        }
        
        
#        save(fitpar,parfit,fitstd,kTT,kTK,kTS,B,intersect.latlon,intersect.xyz,acf,var,file=paste('beatatest',as.integer(data.T$parbl[11,max(columns.T)]),'.Rdata',sep=''))
        print(fname <- paste('beatatest',as.integer(data.T$parbl[11,max(columns.T)]),'.Rdata',sep=''))
        save(fitpar,initParam,apriori,limitParam,parScales,kTT,kTK,kTS,B,intersect.latlon,intersect.xyz,acf,var,file=fname)
    }
    iper <- iper + 1
  }
  
  
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

readTestRes <- function(ddir='.'){
    f <- dir(ddir,pattern='beata',full.names=T)
    nf <- length(f)
    dmat <- emat <- matrix(ncol=15,nrow=nf)
    for(k in seq(nf)){
        load(f[k])
        dmat[k,] <- parfit
        emat[k,] <- fitstd
    }
    return(list(dmat=dmat,emat=emat))
}
