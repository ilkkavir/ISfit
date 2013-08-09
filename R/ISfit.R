ISfit <- function( ddirs='.' , odir='.' , llhT=NA , azelT=NA , llhR=NA , azelR=NA , freq.Hz=NA , rangeLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=100 , plotTest=FALSE , plotFit=FALSE , acfScale=1e-13)
  {
      #
      #
      # This function will be converted into ISfit.LPI that can handle all kinds
      # of fits with LPI input files.
      #
      #
      #
      #
      #
    # 1D incoherent scatter plasma parameter fit using LPI output files in ddirs
    # 
    # Only the line-of-sight velocity will be estimated.
    # This is a GUISDAP-style fit, which uses a single ion temperature and collision frequency. Currently 3 ions ( 30.5 , 16.0 , 1.0 )
    # 
    #
    #
    # INPUT:
    #   ddirs          a vector of data directory paths
    #   odir           Output directory
    #   llhT           latitude, longitude, and height of the transmitter site  ( if NA, values will be read from data files)
    #   azelT          azimuth and elevation of the transmitter beam            ( if NA, values will be read from data files)
    #   llhR           latitude, longitude, and height of the receiver site     ( if NA, values will be read from data files)
    #   freq.Hz        transmitter carrier frequency
    #   rangeLimits.km analysis range-gate limits. If NA, the range-resolution of LPI will be used
    #   timeRes.s      time resolution (integration time)
    #   beginTime      c(year,month,day,hour,minute,seconds) analysis start time
    #   endTime        c(year,month,day,hour,minute,seconds) analysis end time
    #   absLimit       limit for absolute value of the residual. The iteration will not be stopped (unles maxIter is reached) before the residual is below absLimit
    #   diffLimit      Upper limit for fractional change in residual in an iteration step.
    #   maxLambda      maximum Lambda value in Levenberg-Marquardt iteration
    #   maxIter        maximum number of iterations
    #   plotTest       if TRUE, the direct theory and measurements are plotted at every single iteration step
    #   acfScale       ACF scaling factor
    #
    # OUTPUT:
    #   None, the results are written to files in odir.
    #

    # create the output directory
    dir.create( odir , recursive=TRUE , showWarnings=FALSE) 
    
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
          dlist[[n]][["nGates"]] <- rep(length(ACF$range.km),length(ACF$lag.us))
         }
        
        # read acf, variance, lag, and range of each data point in a vector
        acf <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["ACF"]] , n=x[["nGates"]] ) ) ) } ) )
        var <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["var"]] , n=x[["nGates"]] ) ) ) } ) )
        lag <- unlist( lapply( dlist , function(x){ return( rep( x[["lag.us"]] , times=x[["nGates"]] ) ) } ) )
        ran <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] ] ) } , x=x[["range.km"]] , n=x[["nGates"]] ) ) ) } ) )

        # average positions and pointing directions during the integration period (well.. hope the antennas won't move.. but let's average anyway..)
        llhTf  <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["llhT"]] ) } ) )  , ncol=3 , byrow=TRUE ) )
        llhRf  <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["llhR"]] ) } ) )  , ncol=3 , byrow=TRUE ) )
        azelTf <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["azelT"]] ) } ) ) , ncol=2 , byrow=TRUE ) )
        azelRf <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["azelR"]] ) } ) ) , ncol=2 , byrow=TRUE ) )
        freqTf <- mean(unlist( lapply( dlist , function(x){ return( x[["radarFreq"]] ) } ) ) )


#        # positions and pointing directions 
#        llhTf  <- matrix( unlist( lapply( dlist , function(x){ return( x[["llhT"]] ) } ) )  , ncol=3 , byrow=TRUE ) 
#        llhRf  <- matrix( unlist( lapply( dlist , function(x){ return( x[["llhR"]] ) } ) )  , ncol=3 , byrow=TRUE ) 
#        azelTf <- matrix( unlist( lapply( dlist , function(x){ return( x[["azelT"]] ) } ) ) , ncol=2 , byrow=TRUE ) 
#
#        # replace with user input if that is given
#        if()
#
#        # unique values
#        llhTsite <- unique(llhTf)
#        llhRsite <- unique(llhRf)
#        azelTsite <- unique(azelTf)
#       
#
#
#
#       
#
#
#
        

        # a time vector converted from iperLimits
        t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='ut')
        date <- c(t$year+1900,t$mon,t$mday,t$hour,t$min,t$sec)

        # range-gate limits
        if( all( is.na( rangeLimits.km ) ) ){
          rlims <- c( sort( unique( ran ) ) , Inf )
        }else{
          rlims <- unique( rangeLimits.km )
        }

        nr <- length( rlims ) - 1


        # allow different ACF scale in each range gate
        acfScale <- rep(acfScale,length.out=nr)
        
        # select the values that will be actually used
        llhTu <- llhT
        llhRu <- llhR
        azelTu <- azelT
        azelRu <- azelR
        freqTu <- freq.Hz
        if( any( is.na( llhT ) ) ) llhTu <- llhTf
        if( any( is.na( llhR ) ) ) llhRu <- llhRf
        if( any( is.na( azelT ) ) ) azelTu <- azelTf
        if( any( is.na( azelR ) ) ) azelRu <- azelRf
        if( any( is.na( freq.Hz ) ) ) freqTu <- freqTf

        covar <- intersect <- list()
        for(r in seq(nr)) covar[[r]] <- matrix(ncol=13,nrow=13)

        model <- std <- param <- matrix(ncol=13,nrow=nr)
        latitude <- longitude <- range <- height <- status <- chisqr <- rep(-1,nr)
        B <- matrix(ncol=3,nrow=nr)
          
        for( r in seq( nr ) ){

          gateinds       <- which( ( ran >= rlims[r] ) & (ran < rlims[r+1]) )

          acf.gate       <- acf[ gateinds ]
          var.gate       <- var[ gateinds ]
          lag.gate       <- lag[ gateinds ]
          ran.gate       <- ran[ gateinds ]

          # remove NA values
          nainds         <- is.na(acf.gate) | is.na(var.gate)
          if(any(!nainds)){
            acf.gate       <- acf.gate[ !nainds ] * acfScale[r]
            var.gate       <- var.gate[ !nainds ] * acfScale[r]**2
            lag.gate       <- lag.gate[ !nainds ] * 1e-6
            ran.gate       <- ran.gate[ !nainds ]
            r.gate         <- mean( range( ran.gate ) )

            # exact coordinates of the measurement volume
            llhTarget      <- range2llh( azelT=azelTu , llhR=llhRu ,  llhT=llhTu ,  r=r.gate * 1000 )
            height[r]      <- llhTarget['h'] / 1000
            latitude[r]    <- llhTarget['lat']
            longitude[r]   <- llhTarget['lon']
            range[r]       <- r.gate

            # azelR is an emergency fix!!!...
            intersect[[r]] <- beamIntersection( llhT=llhTu , llhR=llhRu , azelT=azelTu , azelR=c(313.95,52+r*2) , fwhmT=1 , fwhmR=2 , phArrT=FALSE , phArrR=TRUE , freq.Hz=freqTu  )

            
print(            (rlims[r:(r+1)]*1000 - intersect[[r]][["range"]]["R"]) * sin( intersect[[r]][["phi"]]*pi/360)/1000)
print(            gainR <- gategain( intersect[[r]] , rlims[r:(r+1)]*1000))
            acf.gate <- acf.gate / gainR
            var.gate <- var.gate / gainR**2

            # magnetic field direction
            Btmp           <- igrf(date=date[1:3],lat=latitude[r],lon=longitude[r],height=height[r],isv=0,itype=1)
            B[r,]          <- c(Btmp$y,Btmp$x,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards east,

            # parameters from iri model
            ptmp           <- iriParams( time=date ,latitude=llhTarget[1],longitude=llhTarget[2],heights=llhTarget[3]/1000 )

            # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
            # the densities in outfmsis are in cm^-3
            ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,1])['NO+',] )

          
            # initial plasma parameter values, there should not be negative ones, as ion velocity is set to zero
            parInit        <- pmax( c( ptmp['e-',1] , ptmp['Ti',1] , ptmp['Ti',1], ptmp['Te',1] , ptmp['Te',1] , ioncoll , 0 , 0 , 0 , ifelse( (sum(ptmp[c('O2+','NO+'),1])<0) , 0 , sum(ptmp[c('O2+','NO+'),1])/ptmp['e-',1]) , ifelse( (ptmp['O+',1]<0) , 0 , ptmp['O+',1]/ptmp['e-',1]) , ifelse( (ptmp['H+',1]<0) , 0 , ptmp['H+',1]/ptmp['e-',1] ) , 1 ) , 0 )
            parInit[1]     <- max(parInit[1],1e6)

            # parameter scaling factors
            parScales      <- ISparamScales(parInit,3)
            # scale the initial parameter values
            initParam      <- scaleParams( parInit , parScales , inverse=F)
            # parameter value limits
            parLimits      <- ISparamLimits(3,1)
            # scale the parameter limits
            limitParam     <- parLimits
            limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
            limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
            # apriori information
            apriori        <- ISapriori( initParam , 3 , TRUE , TRUE )

            model[r,]      <- parInit

            fitpar   <- ISparamfit(
              acf             = acf.gate,
              var             = var.gate,
              nData           = length(acf.gate),
              fSite           = freqTu,
              aSite           = scattVector.llhazelr( llhTu , azelTu , llhRu , r.gate , freqTu )[["phi"]],
              kSite           = scattVector.llhazelr( llhTu , azelTu , llhRu , r.gate , freqTu )[["k.ENU"]],
              iSite           = 1,
              B               = B[r,],
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
              lags            = lag.gate,
              plotTest        = plotTest,
              plotFit         = plotFit,
              maxLambda       = maxLambda,
              maxIter         = maxIter
              )

            # scale back to physical units
            param[r,] <- scaleParams( fitpar$param , scale=parScales , inverse=TRUE )
            covar[[r]] <- scaleCovar( fitpar$covar , scale=parScales , inverse=TRUE)
            std[r,]   <- sqrt(diag(covar[[r]]))#scaleParams( sqrt(diag(fitpar$covar)) , scale=parScales , inverse=TRUE )
            chisqr[r] <- fitpar[["chisqr"]]
            status[r] <- fitpar[["fitStatus"]]
            dimnames(covar[[r]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'))
          }
        }

        time_sec <- iperLimits[k+1]
        POSIXtime <- as.POSIXlt(time_sec,origin='1970-01-01',tz='ut')
        PPI_param <- list(mi = c( 30.5 , 16.0 , 1.0 ) )
        std[is.na(std)] <- Inf

        # another bubble gum fix, will be replaced in final version...
        dimnames(param) <- list(paste('gate',seq(nr),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'))
        dimnames(std)   <- list(paste('gate',seq(nr),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'))
        dimnames(model)   <- list(paste('gate',seq(nr),sep=''),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'))
        names(range) <- paste('gate',seq(nr),sep='')
        names(height) <- paste('gate',seq(nr),sep='')
        names(latitude) <- paste('gate',seq(nr),sep='')
        names(longitude) <- paste('gate',seq(nr),sep='')
        dimnames(B) <- list(paste('gate',seq(nr),sep=''),c('x','y','z'))
        
        
        # save the results to file
        PP <- list(param=param,std=std,model=model,chisqr=chisqr,status=status,time_sec=time_sec,date=date,POSIXtime=POSIXtime,range=range,height=height,latitude=latitude,longitude=longitude,azelT=azelTu,llhT=llhTu,llhR=llhRu,acfScale=acfScale,intersect=intersect,covar=covar,B=B)
        resFile <- file.path( odir , paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "PP.Rdata" , sep=''))
        save( PP , PPI_param, file=resFile )

        cat(iperLimits[k+1],'\n')       

        
      }

    }


  }
