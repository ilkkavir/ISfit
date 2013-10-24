ISmeas.simu <- function(refPoint=KIR,locTrans=TRO,locRec=list(TRO,KIR,KIL),locTarg=c(TRO,200),locxy=FALSE,fwhmTrans=c(1),fwhmRec=c(1,2,4),fwhmRange=c(1),resNS,resEW,resH,Pt=c(1e6),Tnoise=c(200),fradar=233e6,phArrTrans=c(FALSE),phArrRec=c(FALSE,FALSE,TRUE),ele=c(1e11,300,300,0,0,0,0),ion=list(c(30.5,.7e11,300,300,0,0,0,0),c(16,.3e11,300,300,0,0,0,0),c(1,0,300,300,0,0,0,0)),freq=seq(-1000,1000)*100,lags=seq(50)*50e-6,integrationTime=10,dutyCycle=.12,time=c(2012,6,1))
  {
    # simulated incoherent scatter ACFs for a given measurement geometry and plasma parameters
    #
    #
    #
    #
    # I. Virtanen 2012

  # if only one location is given as a vector, convert it into a list
  if(!is.list(locTrans)) locTrans <- list(locTrans)
  if(!is.list(locRec))   locRec   <- list(locRec)
  

  nTrans <- length(locTrans)
  nRec   <- length(locRec)
  nComb <- nRec * nTrans
  
  # make sure that all locations are available as both (lat,lon,height) and (x,y,z)
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
  # target, we will need also the (x,y,h) coordinates
  if(locxy){
    latlonTarg   <- ISgeometry:::planarToSpherical.geographic(x=locTarg[1],y=locTarg[2],refPoint=refPoint)
  }else{
    latlonTarg   <- list(lat=locTarg[1],lon=locTarg[2])
  }
  xyhTarg        <- c(unlist(ISgeometry:::sphericalToPlanar.geographic(lat=latlonTarg$lat,lon=latlonTarg$lon,zeroLatitude=refPoint[1],zeroMeridian=refPoint[2]))[c('x','y')],locTarg[3])
  xyzTarg        <- ISgeometry:::sphericalToCartesian(c(latlonTarg$lat,latlonTarg$lon,ISgeometry:::EarthRadius()+locTarg[3]))

  # estimates for signal and noise levels, scattering wave vectors, etc. from ISgeometry:::multistaticNoiseLevels
  # if only one receiver, add a dummy copy of it to fulfil the requirements of multistaticNoiseLevel
  noiseLevel <- ISgeometry:::multistaticNoiseLevels( refPoint=refPoint , locTrans=locTrans , locRec=rep(locRec,length.out=max(2,length(locRec))) , locxy=locxy ,
                                                     fwhmTrans=fwhmTrans , fwhmRec=fwhmRec , fwhmRange=fwhmRange , x=xyhTarg[1] ,
                                                     y=xyhTarg[2] , heights=xyhTarg[3] , infinity=1e3 , resNS=resNS ,
                                                     resEW=resEW , resH=resH , Pt=Pt , Ne=ele[1] , Tnoise=Tnoise , fradar=fradar ,
                                                     tau0=100 , phArrTrans=phArrTrans , phArrRec=phArrRec )


  # baud length from fwhmRange
  blen <- fwhmRange*2/299792.458
  # number of samples integrated in one integration period
  nint <- integrationTime/blen*dutyCycle

  # scattering angles
  scattAngle <- vector(mode='numeric',length=nComb)

  for(t in seq(nTrans)){
    for(r in seq(nRec)){
      # scattering angle (angle between incident and scattered waves)
      scattAngle[(t-1)*nTrans+r]   <- ISgeometry:::vectorAngle.cartesian((xyzTarg-xyzTrans[[t]]),(xyzRec[[r]]-xyzTarg),degrees=T)
    }
  }

  
  # signal-to-noise ratios and scattering wave vectors for each transmitter - receiver combination
  snr <- rep(0,nComb)
  kscat <- rep(list(c(0,0,0),nComb))
  for(k in seq(nComb)){
    snr[k] <- noiseLevel$noiseLevel.site[[k]]$noiseLevel$Pr / ( noiseLevel$noiseLevel.site[[k]]$noiseLevel$Pn + noiseLevel$noiseLevel.site[[k]]$noiseLevel$Prn )
    kscat[[k]] <- noiseLevel$noiseLevel.site[[k]]$kscat[1,1,,]
  }

  # ratio of lag profile standard deviation and zero-lag power after time integration
  sdacf <- 1/(2*snr*sqrt(nint))

  # magnetic field direction
  Btmp           <- igrf(date=time[1:3],lat=latlonTarg[['lat']],lon=latlonTarg[['lon']],height=locTarg[3],isv=0,itype=1)
  B              <- c(Btmp$x,-Btmp$y,-Btmp$z) # the model has y-axis to east and z-axis downwards

  # a list of the simulated ACFs
  acflist <- vector(mode='list',length=nComb)
  for(k in seq(nComb)){
    acflist[[k]] <- simuACF( ele=ele , ion=ion , kdir=kscat[[k]] , fradar=fradar , scattAngle=scattAngle[k] , freq=freq, lags=c(0,lags) , Bdir=B )
  }

  # add noise and strip off the additional zero-lag
  nlags <- length(lags)
  varlist <- vector(mode='list',length=nComb)
  for( k in seq(nComb)){
    varlist[[k]] <- rep(Re(acflist[[k]][1])**2*sdacf[k]**2/2,nlags)
    acflist[[k]] <- acflist[[k]][2:(nlags+1)] + ( rnorm(nlags) + 1i*rnorm(nlags) ) * sqrt(varlist[[k]][1]/2)
  }


  return(list(ACF=acflist,var=varlist,sdacf=sdacf))
  
  }
