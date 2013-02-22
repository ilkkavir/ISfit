ISfit.1D.guisdap <- function( ddirs='.' , odir='.' , llhT=c(69.58,19.23,86.00) , azelT=c(0,90) , llhR=c(69.58,19.23,86.00) , freq.Hz=224e6 , rangeLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=100 , plotTest=FALSE , acfScale=1e-13)
  {
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
        }

        # read acf, variance, lag, and range of each data point into a vector
        acf <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["ACF"]] , n=x[["nGates"]] ) ) ) } ) )
        var <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] , i ] ) } , x=x[["var"]] , n=x[["nGates"]] ) ) ) } ) )
        lag <- unlist( lapply( dlist , function(x){ return( rep( x[["lag"]] , times=x[["nGates"]] ) ) } ) )
        ran <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( x[ 1 : n[i] ] ) } , x=x[["range"]] , n=x[["nGates"]] ) ) ) } ) )

        # average positions and pointing directions during the integration period (well.. hope the antennas won't move.. but let's average anyway..)
        llhTf  <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["llhT"]] ) } ) )  , ncol=3 , byrow=TRUE ) )
        llhRf  <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["llhR"]] ) } ) )  , ncol=3 , byrow=TRUE ) )
        azelTf <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["azelT"]] ) } ) ) , ncol=2 , byrow=TRUE ) )

        # select the values that will be actually used
        llhTu <- llhT
        llhRu <- llhR
        azelTu <- azelT
        if( any( is.na( llhT ) ) ) llhTu <- llhTf
        if( any( is.na( llhR ) ) ) llhRu <- llhRf
        if( any( is.na( azelT ) ) ) azelTu <- azelTf

        # a time vector converted from iperLimits
        t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='ut')
        date <- c(t$year+1900,t$mon,t$mday,t$hour,t$min,t$sec)
#        date <- as.numeric(c(substr(t,1,4),substr(t,6,7),substr(t,9,10),substr(t,12,13),substr(t,15,16),substr(t,18,19)))

        # range-gate limits
        if( all( is.na( rangeLimits.km ) ) ){
          rlims <- c( sort( unique( ran ) ) , Inf )
        }else{
          rlims <- unique( rangeLimits.km )
        }

        nr <- length( rlims ) - 1

        model <- std <- param <- matrix(ncol=7,nrow=nr)
        latitude <- longitude <- range <- height <- status <- chisqr <- rep(-1,nr)
          
        for( r in seq( nr ) ){

          gateinds       <- which( ( ran >= rlims[r] ) & (ran < rlims[r+1]) )

          acf.gate       <- acf[ gateinds ]
          var.gate       <- var[ gateinds ]
          lag.gate       <- lag[ gateinds ]
          ran.gate       <- ran[ gateinds ]

          # remove NA values
          nainds         <- is.na(acf.gate) | is.na(var.gate)
          acf.gate       <- acf.gate[ !nainds ] * acfScale
          var.gate       <- var.gate[ !nainds ] * acfScale**2
          lag.gate       <- lag.gate[ !nainds ] * 1e-6
          ran.gate       <- ran.gate[ !nainds ]
          r.gate         <- mean( range( ran.gate ) )

          # exact coordinates of the measurement volume
          llhTarget      <- range2llh( azelT=azelTu , llhR=llhRu ,  llhT=llhTu ,  r=r.gate * 1000 )
          height[r]      <- llhTarget['h'] / 1000
          latitude[r]    <- llhTarget['lat']
          longitude[r]   <- llhTarget['lon']
          range[r]       <- r.gate
          

          # parameters from iri model
          ptmp           <- iriParams( time=date ,latitude=llhTarget[1],longitude=llhTarget[2],heights=llhTarget[3]/1000 )

          # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
          # the densities in outfmsis are in cm^-3
          ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,1])['NO+',] )

          # initial plasma parameter values, there should not be negative ones, as ion velocity is set to zero
          parInit        <- pmax( c( ptmp['e-',1] , ptmp['Ti',1] , ptmp['Te',1]/ptmp['Ti',1] , ioncoll , 0 , ifelse( (ptmp['O+',1]<0) , 0 , ptmp['O+',1]/ptmp['e-',1]) , ifelse( (ptmp['H+',1]<0) , 0 , ptmp['H+',1]/ptmp['e-',1] )) , 0 )
          parInit[1]     <- max(parInit[1],1e6)

          model[r,]      <- parInit

          # parameter scaling factors
          parScales      <- parInit
          parScales[3]   <- 1
          parScales[5:7] <- 1

          # scale the initial parameter values
          parInit        <- scaleParams( parInit , parScales , inverse=F)

          # parameter value limits
          parLimits      <- matrix( c( 1e4 , 10 , 1e-2 , 0 , -1e4 , 0 , 0    ,     1e13 , 1e4 , 100 , 1e20 , 1e4 , 1 , 1) , nrow=2 , byrow=T )

          # scale the parameter limits
          limitParam     <- parLimits
          limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
          limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
  
          # apriori information
          apriori        <- ISapriori.default.guisdap( parInit , 3 )

          fitpar   <- ISparamfit.guisdap(
                        acf             = acf.gate,
                        var             = var.gate,
                        nData           = length(acf.gate),
                        fSite           = freq.Hz,
                        aSite           = scattVector.llhazelr( llhTu , azelTu , llhRu , r.gate , freq.Hz )[["phi"]],
                        initParam       = parInit,
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
                        lags            = lag.gate,
                        plotTest        = plotTest,
                        maxLambda       = maxLambda,
                        maxIter         = maxIter
                        )

          # scale back to physical units
          param[r,] <- scaleParams( fitpar$param , scale=parScales , inverse=TRUE )
          std[r,]   <- scaleParams( sqrt(diag(fitpar$covar)) , scale=parScales , inverse=TRUE )
          chisqr[r] <- fitpar[["chisqr"]]
          status[r] <- fitpar[["fitStatus"]]

        }

        time_sec <- iperLimits[k+1]
        POSIXtime <- as.POSIXlt(time_sec,origin='1970-01-01',tz='ut')
        PPI_param <- list(mi = c( 30.5 , 16.0 , 1.0 ) )
        std[is.na(std)] <- Inf
        
        # save the results to file
        PP <- list(param=param,std=std,model=model,chisqr=chisqr,status=status,time_sec=time_sec,date=date,POSIXtime=POSIXtime,range=range,height=height,latitude=latitude,longitude=longitude,azelT=azelTu,llhT=llhTu,llhR=llhRu)
        resFile <- file.path( odir , paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "PP.Rdata" , sep=''))
        save( PP , PPI_param, file=resFile )

        cat(iperLimits[k+1],'\n')       

        
      }

    }


  } # ISfit.1D.guisdap
