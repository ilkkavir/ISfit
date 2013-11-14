ISfit.3D <- function( ddirs='.' , odir='.' ,  heightLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=10 , plotTest=FALSE , plotFit=FALSE , caltable='caltable.Rdata' , absCalib=FALSE , TiIsotropic)
  { 
      
      # 3D incoherent scatter plasma parameter fit using LPI output files in ddirs
      # 
      # This is a GUISDAP-style fit, which uses a single ion temperature and collision frequency. Currently 3 ions ( 30.5 , 16.0 , 1.0 )
      # 
      #
      #
      # INPUT:
      #   ddirs           a vector of data directory paths
      #   odir            Output directory
      #   heightLimits.km analysis height-gate limits. If NA, range gates of the reference site are used
      #   timeRes.s       time resolution (integration time)
      #   beginTime       c(year,month,day,hour,minute,seconds) analysis start time
      #   endTime         c(year,month,day,hour,minute,seconds) analysis end time
      #   absLimit        limit for absolute value of the residual.
      #                   The iteration will not be stopped (unles maxIter is reached) before the residual is below absLimit
      #   diffLimit       Upper limit for fractional change in residual in an iteration step.
      #   maxLambda       maximum Lambda value in Levenberg-Marquardt iteration
      #   maxIter         maximum number of iterations
      #   plotTest        if TRUE, the direct theory and measurements are plotted at every single iteration step
      #   caltable        full path to a calibration table file
      #   absCalib        TRUE if the remotes are absolutely calibrated, FALSE to allow for scaling of their calibration coefficients
      #   TiIsotropic     TRUE if ion thermal velocity distribution is modeled as isotropic, FALSE if bi-maxwellian
      #
      #   a calibration table file contains a matrix named 'caltable' with following columns:
      #      TXlatitude  transmitter latitude in degrees
      #      TXlongitude transmitter longitude in degrees
      #      TXfwhm      full width at half maximum of the transmitter beam when pointed to zenith
      #      TXphArr     TRUE if the transmitter is a phased array, FALSE for a dish
      #      TXazimuth   transmitter beam azimuth in degrees [0,360], NA if same calibration for all azimuths
      #      TXelevation transmitter beam elevation in degrees [0,90], NA if same calibration for all elevations
      #      TXfreq      transmitter frequency in Hz
      #      RXlatitude  receiver latitude in degrees
      #      RXlongitude receiver longitude in degrees
      #      RXfwhm      full width at half maximum of the receiver beam *when pointed to zenith*
      #      RXphArr     TRUE if the receiver is a phased array, FALSE for a dish
      #      RXazimuth   receiver beam azimuth in degrees [0,360], NA if same calibration for all azimuths
      #      RXelevation receiver beam elevation in degrees [0,90], NA if same calibration for all elevations
      #      Conjugate   if TRUE, the ACF samples will be complex conjugated to compensate for spectrum mirroring in certain cases
      #      N           number of range points at which the calibration coefficients are given
      #      ranges      ranges of the point in metres, range is one half of signal path length from
      #                  the transmitter, via the target, to the receiver
      #      calconst    Numerical values of calibration constants
      #
      #   Linear interpolation will be used in between the range points, the same constant is used at all ranges if N==1
      #
      #   Radar site identification uses 0.1 degree accuracy in latitude, longitude, azimuth, and elevation,
      #   and 10 MHz accuracy in frequency.
      #
      #
      #   The calibration table is copied to odir 
      #
      # OUTPUT:
      #   None, the results are written to files in odir.
      #

      # load the calibration table
      load(caltable)


      # create the output directory
      dir.create( odir , recursive=TRUE , showWarnings=FALSE)
      
      # save a copy of the calibration table
      fname <- paste('caltable_',format(Sys.time(),"%Y-%m-%d_%H:%M"),'.Rdata',sep='')
      save( caltable , file=file.path(odir,fname))
    
      # list all data files
      dfiles <- c()
      for( dn in ddirs ) dfiles <- c( dfiles , dir( dn , full.names=T) )
      dfiles <- unique( dfiles[ grep( 'LP.Rdata' , dfiles ) ] )

      # stop if the data file vector is empty
      if ( length( dfiles ) == 0 ) stop( "Could not find any data files." )

      # read timestamps from file names. This is the unix time at end of integration period
      tstamps <- as.numeric(substr(dfiles,nchar(dfiles)-20,nchar(dfiles)-8)) / 1000

      # order the data files according to their timestamps
      torder  <- order(tstamps)
      dfiles  <- dfiles[ torder ]
      tstamps <- tstamps[ torder]

      # convert beginTime and endTime into unix time format
      bTime <- as.double(ISOdate(beginTime[1],beginTime[2], beginTime[3],beginTime[4], beginTime[5],beginTime[6])) + beginTime[6]%%1
      eTime <- as.double(ISOdate(endTime[1],endTime[2], endTime[3],endTime[4], endTime[5],endTime[6])) + endTime[6]%%1

      # first integration period from which we actually have data
      iperFirst <- max( 1, floor( ( tstamps[1] - bTime ) / timeRes.s ) )
    
      # last integration period
      iperLast <- ceiling( ( min( tstamps[length(tstamps)] , eTime ) - bTime ) / timeRes.s ) 

      # stop if analysis end time is before its start time
      if( iperLast < iperFirst ) stop()
    
      # integration period limits
      iperLimits <- seq( iperFirst , iperLast ) * timeRes.s + bTime

      # number of integration periods
      nIper <- length(iperLimits) - 1

      # walk through all integration periods
      for( k in seq( nIper ) ){

          fitPar <- list()
      
          # look for data files from the current integration period
          iperFiles <- dfiles[ which( ( tstamps > iperLimits[k] ) & ( tstamps <= iperLimits[ k + 1 ] ) ) ]

          # load all data and collect it in a list
          nFiles <- length(iperFiles)
          if( nFiles > 0 ){
              dlist <- vector( mode='list', length=nFiles)
              for( n in seq( nFiles ) ){
                  load( iperFiles[n] )
                  dlist[[n]] <- ACF
                  if(length(ACF[["nGates"]])!=length(ACF$lag.us)) dlist[[n]][["nGates"]] <- rep(length(ACF$range.km),length(ACF$lag.us))
              }

              # read acf, variance, lag, range, pointing directions, and TX / RX location of each data point
              acf    <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["ACF"]] , n=x[["nGates"]] ) ) ) } ) )
              var    <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["var"]] , n=x[["nGates"]] ) ) ) } ) )
              lag    <- 1e-6*unlist( lapply( dlist , function(x){ return( rep( x[["lag.us"]] , times=x[["nGates"]] ) ) } ) )
              ran    <- 1000*unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] ] ) } , x=x[["range.km"]] , n=x[["nGates"]] ) ) ) } ) )
              azTf   <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,1] , length.out=n[i] ) ) } , x=x[["azelT"]] , n=x[["nGates"]] ) ) ) } ) )
              elTf   <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,2] , length.out=n[i] ) ) } , x=x[["azelT"]] , n=x[["nGates"]] ) ) ) } ) )
              azRf   <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,1] , length.out=n[i] ) ) } , x=x[["azelR"]] , n=x[["nGates"]] ) ) ) } ) )
              elRf   <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,2] , length.out=n[i] ) ) } , x=x[["azelR"]] , n=x[["nGates"]] ) ) ) } ) )
              latTf  <- unlist( lapply( dlist , function(x){ return( rep( x[["llhT"]][1] , sum(x[["nGates"]] ) ) ) } ) )
              lonTf  <- unlist( lapply( dlist , function(x){ return( rep( x[["llhT"]][2] , sum(x[["nGates"]] ) ) ) } ) )
              hTf    <- unlist( lapply( dlist , function(x){ return( rep( x[["llhT"]][3] , sum(x[["nGates"]] ) ) ) } ) )
              latRf  <- unlist( lapply( dlist , function(x){ return( rep( x[["llhR"]][1] , sum(x[["nGates"]] ) ) ) } ) )
              lonRf  <- unlist( lapply( dlist , function(x){ return( rep( x[["llhR"]][2] , sum(x[["nGates"]] ) ) ) } ) )
              hRf    <- unlist( lapply( dlist , function(x){ return( rep( x[["llhR"]][3] , sum(x[["nGates"]] ) ) ) } ) )
              freqTf <- unlist( lapply( dlist , function(x){ return( rep( x[["radarFreq"]] , sum(x[["nGates"]]) ) ) } ) )
              
              # all necessary information about each data point
              # carrier frequency, TX latitude, TX longitude, TX elllipsoid height, TX azimuth, TX elevation, and same for RX
              # the first column will be site index
              datasites <- matrix( c( rep(0,length(freqTf)) , freqTf , latTf , lonTf , hTf , azTf , elTf , latRf , lonRf , hRf , azRf , elRf ) , ncol=12)

              # find unique combinations of llhT, llhR, azelT, azelR, and freqT. Each unique combination will be
              # be considered as a "site" in the analysis
              sites <- unique(datasites)
        
              # number of sites
              nsites <- dim(sites)[1]

              # generate the site indices
              for( s in seq(nsites)){
                  datasites[ which( apply( datasites , FUN=function(x,y){all(x==y)} , MARGIN=1 , y=sites[s,] ) ) , 1 ] <- s
              }
              sites[,1] <- seq(nsites)

              # select the reference site
              # if there is only one site this is trivial
              # if absCalib==TRUE we can just pick any site
              if( nsites==1 | absCalib ){ 
                  refsite <- 1
              }else{ # otherwise look for a monostatic site
                  refsite <- which( apply( sites[,3:7] == sites[,8:12] , FUN=all , MARGIN=1 ) )
                  # if there are none or several candidates we will need to ask the user
                  # ask only at the first period, the same choice will be used in all remaining periods!
                  if(k==1){
                      if( length(refsite) != 1){
                          if( length(refsite)==0){
                              cat('Could not find a monostatic site,\n')
                              cat('available sites are:\n')
                              cat('Site    TXfreq     TXlat     TXlon  TXheight      TXaz     TXele     RXlat     RXlon  RXheight      RXaz     RXele\n')
                              for(s in seq(nsites)){
                                  cat(sprintf("[%2i]",s),sprintf(" %10.2f",sites[s,2:12]),sprintf("\n"))
                              }
                              refsite <- as.numeric(readline("Please select number of the reference site: "))
                          }else{
                              cat('Found several monostatic sites,\n')
                              cat('available sites are:\n')
                              cat('Site    TXfreq     TXlat     TXlon  TXheight      TXaz     TXele     RXlat     RXlon  RXheight      RXaz     RXele\n')
                              for(s in refsite){
                                  cat(sprintf("[%2i]",s),sprintf(" %10.2f",sites[s,2:12]),sprintf("\n"))
                              }
                              refsite <- as.numeric(readline("Please select number of the reference site: "))
                          }
                      }
                  }
              }

              # scale all data with values tabulated in caltable
              dscales <- acfscales( datasites , ran , caltable )
              acf <- acf * dscales[,6]
              acf[dscales[,5]>0] <- Conj(acf[dscales[,5]>0])
              var <- var * dscales[,6]**2

              # a time vector converted from iperLimits
              t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='ut')
              date <- c(t$year+1900,t$mon,t$mday,t$hour,t$min,t$sec)
              
              # height gate limits
              if( all( is.na( heightLimits.km ) ) ){
                  rlims <- sort( unique( ran[which(datasites[,1]==refsite)] ) )
                  rlims <- c( rlims , max(rlims) + 1)
                  hlims <- rlims
                  for( h in seq(length(hlims))) hlims[h] <- range2llh(r=rlims[h],llhT=sites[refsite,3:5],llhR=sites[refsite,8:10],azelT=sites[refsite,6:7])['h']
              }else{
                  hlims <- unique( heightLimits.km )*1000
              }
              nh <- length( hlims ) - 1
              
              covar <- intersect <- list()
              for(h in seq(nh)) covar[[h]] <- matrix(ncol=12+nsites,nrow=12+nsites)
              
              model <- std <- param <- matrix(ncol=12+nsites,nrow=nh)
              latitude <- longitude <- height <- status <- chisqr <- rep(-1,nh)
              B <- matrix(ncol=3,nrow=nh)

              # convert all ranges to latitude, longitude, height
              llh <- matrix(nrow=length(ran),ncol=3)
              for( dind in seq(length(ran))){
                  llh[dind,] <- range2llh( r=ran[dind] , llhT=datasites[dind,3:5] , llhR=datasites[dind,8:10] , azelT=datasites[dind,6:7])
              }
              
              for( h in seq( nh ) ){
                  
                  # Initial values, these will be immediately updated if any data is found
                  height[h]      <- sum(hlims[h:(h+1)])/2000
                  latitude[h]    <- sites[refsite,3]
                  longitude[h]   <- sites[refsite,4]
                  Btmp           <- igrf(date=date[1:3],lat=latitude[h],lon=longitude[h],height=height[h],isv=0,itype=1)
                  B[h,]          <- c(Btmp$y,Btmp$x,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards east,
                  dimnames(covar[[h]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')))
                  
                  # if there are actual data, all initializations will be updated
                  
                  gateinds       <- which( ( llh[,3] >= hlims[h] ) & (llh[,3] < hlims[h+1]) )
                  
                  acf.gate       <- acf[ gateinds ]
                  var.gate       <- var[ gateinds ]
                  lag.gate       <- lag[ gateinds ]
                  ran.gate       <- ran[ gateinds ]
                  llh.gate       <- llh[ gateinds , ]
                  sites.gate     <- datasites[ gateinds , ]
            
                  # remove NA values
                  nainds         <- is.na(acf.gate) | is.na(var.gate)
          
                  if(any(!nainds)){
                                            
                      acf.gate       <- acf.gate[ !nainds ]
                      var.gate       <- var.gate[ !nainds ]
                      lag.gate       <- lag.gate[ !nainds ]
                      llh.gate       <- llh.gate[ !nainds , ]
                      sites.gate     <- sites.gate[ !nainds , ]
                      
                     # coordinates of the measurement volume
                      llhTarget      <- colMeans( llh.gate )
                      height[h]      <- llhTarget[3] / 1000
                      latitude[h]    <- llhTarget[1]
                      longitude[h]   <- llhTarget[2]

                      
                     # beam intersections
                      intersect[[h]] <- list()
                      gainR <- aSite <- nlags.site <-  rep(NA,nsites)
                      kSite <- list()
                      lag.site <- list()
                      acf.site <- list()
                      var.site <- list()
                      ind.site <- list()
                      for( s in seq(nsites) ){
                          # the beam widths and antenna types are read from caltable and returned to this function in dscales
                          ss <- which( datasites[,1] == s)[1]
                          intersect[[h]][[s]] <- beamIntersection( llhT=sites[s,3:5] , llhR=sites[s,8:10] , azelT=sites[s,6:7] , azelR=sites[s,11:12] , fwhmT=dscales[ss,1] , fwhmR=dscales[ss,3] , phArrT=dscales[ss,2]>0 , phArrR=dscales[ss,4]>0 , freq.Hz=sites[s,2] )
                          
                          # conversion from lat, lon, height to range in this gate
#                          rs1 <- (llhTarget2azelrBeam(llhSite=sites[s,3:5],llhTarget=c(llhTarget[1:2],hlims[h]))['r'] +
#                                  llhTarget2azelrBeam(llhSite=sites[s,8:10],llhTarget=c(llhTarget[1:2],hlims[h]))['r'] ) / 2 
#                          rs2 <- (llhTarget2azelrBeam(llhSite=sites[s,3:5],llhTarget=c(llhTarget[1:2],hlims[h+1]))['r'] +
#                                  llhTarget2azelrBeam(llhSite=sites[s,8:10],llhTarget=c(llhTarget[1:2],hlims[h+1]))['r'] ) / 2
                          rs1 <- height2range( llhT=sites[s,3:5] , azelT=sites[s,6:7] , llhR=sites[s,8:10] , h=hlims[h] )
                          rs2 <- height2range( llhT=sites[s,3:5] , azelT=sites[s,6:7] , llhR=sites[s,8:10] , h=hlims[h+1] )

                          # gain integral
                          gainR[s] <- gategain( intersect[[h]][[s]] , c(rs1,rs2) , maxdev=3)

                          # scattering angle
                          aSite[s] <- intersect[[h]][[s]][["phi"]]

                          # scattering wave vector
                          kSite[[s]] <- intersect[[h]][[s]][["k.ENU"]]
                          
                          if(is.na(gainR[s])){
                              lag.site[[s]] <- c()
                              nlags.site[[s]] <- 0
                              acf.site[[s]] <- c()
                              var.site[[s]] <- c()
                              ind.site[[s]] <- c()
                          }else{

                              # data points from this site
                              data.site <- which(sites.gate[,1]==s)

                              # lags measured at this site
                              lag.site[[s]] <- unique(lag.gate[data.site])

                              nlags.site[s] <- length(lag.site[[s]])

                              # average data points from each lag value
                              acf.site[[s]] <- rep(0+0i,nlags.site[s])
                              var.site[[s]] <- rep(0,nlags.site[s])
                              for( l in seq(nlags.site[s])){
                                  lagind <- which(lag.site[[s]][l]==lag.gate[data.site])
                                  for( lind in lagind ){
                                      acf.site[[s]][l] <- acf.site[[s]][l] + acf.gate[data.site][lind]/var.gate[data.site][lind]
                                      var.site[[s]][l] <- var.site[[s]][l] + 1/var.gate[data.site][lind]
                                  }
                              }
                              var.site[[s]] <- 1/var.site[[s]]
                              acf.site[[s]] <- acf.site[[s]] * var.site[[s]]
                              
                              acf.site[[s]] <- acf.site[[s]] / gainR[s]
                              var.site[[s]] <- var.site[[s]] / gainR[s]**2
                              
                              ind.site[[s]] <- rep(s,nlags.site[s])
                          }

                      }

                      # magnetic field direction
                      Btmp           <- igrf(date=date[1:3],lat=latitude[h],lon=longitude[h],height=height[h],isv=0,itype=1)
                      # the model has y-axis to east and z-axis downwards, we have x towards east and z upwards
                      B[h,]          <- c(Btmp$y,Btmp$x,-Btmp$z) 
                      
                      # parameters from iri model
                      ptmp           <- iriParams( time=date ,latitude=latitude[h],longitude=longitude[h],heights=height[h])

                      # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
                      # the densities in outfmsis are in cm^-3
                      ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,1])['NO+',] )

          
                      # initial plasma parameter values, there should not be negative ones, as ion velocity is set to zero
#                      parInit        <- pmax( c( ptmp['e-',1] , ptmp['Ti',1] , ptmp['Ti',1], ptmp['Te',1] , ptmp['Te',1] , ioncoll , 0 , 0 , 0 , ifelse( (sum(ptmp[c('O2+','NO+'),1])<0) , 0 , sum(ptmp[c('O2+','NO+'),1])/ptmp['e-',1]) , ifelse( (ptmp['O+',1]<0) , 0 , ptmp['O+',1]/ptmp['e-',1]) , ifelse( (ptmp['H+',1]<0) , 0 , ptmp['H+',1]/ptmp['e-',1] ) , rep(1,nsites) ) , 0 )
                      # switched to estimating the molecular ion abundance as 1 - O+ - H+. The above estimation has larger error at low altitudes where heavier ions actually exists
                      parInit        <- pmax( c( ptmp['e-',1] , ptmp['Ti',1] , ptmp['Ti',1], ptmp['Te',1] , ptmp['Te',1] , ioncoll , 0 , 0 , 0 , ifelse( (sum(ptmp[c('H+','O+'),1])/ptmp['e-',1]>=1) , 0 , 1-sum(ptmp[c('H+','O+'),1])/ptmp['e-',1]) , ifelse( (ptmp['O+',1]<0) , 0 , ptmp['O+',1]/ptmp['e-',1]) , ifelse( (ptmp['H+',1]<0) , 0 , ptmp['H+',1]/ptmp['e-',1] ) , rep(1,nsites) ) , 0 )
                      parInit[1]     <- max(parInit[1],1e6)
  
                      # parameter scaling factors
                      parScales      <- ISparamScales(parInit,3)
                          
                      # scale the initial parameter values
                      initParam      <- scaleParams( parInit , parScales , inverse=F)
                          
                      # parameter value limits
                      parLimits      <- ISparamLimits(3,nsites)
                          
                      # scale the parameter limits
                      limitParam     <- parLimits
                      limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
                      limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
                      
                      # apriori information
                      apriori        <- ISapriori( initParam , nIon=3 , absCalib=absCalib , TiIsotropic=TiIsotropic , refSite=refsite )

                      model[h,]      <- parInit

                      fitpar   <- ISparamfit(
                          acf             = unlist(acf.site),
                          var             = unlist(var.site),
                          lags            = unlist(lag.site),
                          nData           = sum(nlags.site),
                          fSite           = sites[,2],
                          aSite           = aSite,
                          kSite           = kSite,
                          iSite           = unlist(ind.site),
                          B               = B[h,],
                          initParam       = initParam,
                          invAprioriCovar = apriori$invAprioriCovar,
                          aprioriTheory   = apriori$aprioriTheory,
                          aprioriMeas     = apriori$aprioriMeas,
                          mIon            = c(30.5,16,1),
                          paramLimits     = limitParam,
                          directTheory    = ISdirectTheory,
                          absLimit        = absLimit,
                          diffLimit       = diffLimit,
                          scaleFun        = scaleParams,
                          scale           = parScales,
                          plotTest        = plotTest,
                          plotFit         = plotFit,
                          maxLambda       = maxLambda,
                          maxIter         = maxIter
                          )
                      
                     # scale back to physical units
                      param[h,] <- scaleParams( fitpar$param , scale=parScales , inverse=TRUE )
                      covar[[h]] <- scaleCovar( fitpar$covar , scale=parScales , inverse=TRUE)
                      std[h,]   <- sqrt(diag(covar[[h]]))
                      chisqr[h] <- fitpar[["chisqr"]]
                      status[h] <- fitpar[["fitStatus"]]
                      dimnames(covar[[h]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')))
                  }
              }
              
              time_sec <- iperLimits[k+1]
              POSIXtime <- as.POSIXlt(time_sec,origin='1970-01-01',tz='ut')
              PPI_param <- list(mi = c( 30.5 , 16.0 , 1.0 ) )
              std[is.na(std)] <- Inf
              

              dimnames(param) <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')))
              dimnames(std)   <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')))
              dimnames(model)   <- list(paste('gate',seq(nh),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),paste('Site',seq(nsites),sep='')))
              names(height) <- paste('gate',seq(nh),sep='')
              names(latitude) <- paste('gate',seq(nh),sep='')
              names(longitude) <- paste('gate',seq(nh),sep='')
              dimnames(B) <- list(paste('gate',seq(nh),sep=''),c('x','y','z'))
              
              
              # save the results to file
              PP <- list(param=param,std=std,model=model,chisqr=chisqr,status=status,time_sec=time_sec,date=date,POSIXtime=POSIXtime,height=height,latitude=latitude,longitude=longitude,sites=sites,intersect=intersect,covar=covar,B=B,heightLimits.km=hlims/1000)
              resFile <- file.path( odir , paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "PP.Rdata" , sep=''))
              save( PP , PPI_param, file=resFile )
              
              cat(iperLimits[k+1],'\n')       
              
              
          }
          
      }   
      
  }
