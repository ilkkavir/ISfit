ISfit <- function( ddirs='.' , odir='.' ,  rangeLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=100 , plotTest=FALSE , plotFit=FALSE , acfScale=1e-13)
  {
    # 1D incoherent scatter plasma parameter fit using LPI output files in ddirs
    # 
    # This is a GUISDAP-style fit, which uses a single ion temperature and collision frequency. Currently 3 ions ( 30.5 , 16.0 , 1.0 )
    # 
    #
    #
    # INPUT:
    #   ddirs          a vector of data directory paths
    #   odir           Output directory
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
        azTf <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,1] , length.out=n[i] ) ) } , x=x[["azelT"]] , n=x[["nGates"]] ) ) ) } ) )
        elTf <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,2] , length.out=n[i] ) ) } , x=x[["azelT"]] , n=x[["nGates"]] ) ) ) } ) )
        azRf <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,1] , length.out=n[i] ) ) } , x=x[["azelR"]] , n=x[["nGates"]] ) ) ) } ) )
        elRf <- unlist( lapply( dlist , function(x){ return( unlist( lapply( seq( ncol( x[["ACF"]] ) ) , function( i , n , x ){ return( rep( matrix( x , ncol=2 )[,2] , length.out=n[i] ) ) } , x=x[["azelR"]] , n=x[["nGates"]] ) ) ) } ) )
        azelTf <- matrix(c(azTf,elTf),ncol=2)
        azelRf <- matrix(c(azRf,elRf),ncol=2)

        # average positions (ok, I don't think that the sites will move but anyway...) and frequencies
        llhT  <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["llhT"]] ) } ) )  , ncol=3 , byrow=TRUE ) )
        llhR  <- colMeans(matrix( unlist( lapply( dlist , function(x){ return( x[["llhR"]] ) } ) )  , ncol=3 , byrow=TRUE ) )
        freqT <- mean(unlist( lapply( dlist , function(x){ return( x[["radarFreq"]] ) } ) ) )


        # a time vector converted from iperLimits
        t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='ut')
        date <- c(t$year+1900,t$mon,t$mday,t$hour,t$min,t$sec)

        # range-gate limits
        if( all( is.na( rangeLimits.km ) ) ){
          rlims <- sort( unique( ran ) )
          rlims <- c( rlims , max(rlims) + 1)
        }else{
          rlims <- unique( rangeLimits.km )
        }
        nr <- length( rlims ) - 1
        
        covar <- intersect <- list()
        for(r in seq(nr)) covar[[r]] <- matrix(ncol=13,nrow=13)

        model <- std <- param <- matrix(ncol=13,nrow=nr)
        latitude <- longitude <- range <- height <- status <- chisqr <- rep(-1,nr)
        B <- matrix(ncol=3,nrow=nr)
        azelT <- azelR <- matrix(ncol=2,nrow=nr)

        
        # allow different ACF scale in each range gate
        acfScale <- rep(acfScale,length.out=nr)

        for( r in seq( nr ) ){

            # Initial values, these will be immediately updated if any data is found
            range[r]       <- sum(rlims[r:(r+1)])/2
            azelT[r,]      <- colMeans(azelTf)
            azelR[r,]      <- colMeans(azelRf)
            llhTarget      <- range2llh( azelT=azelT[r,] , llhR=llhR ,  llhT=llhT ,  r=range[r] * 1000 )
            height[r]      <- llhTarget['h'] / 1000
            latitude[r]    <- llhTarget['lat']
            longitude[r]   <- llhTarget['lon']
            intersect[[r]] <- beamIntersection( llhT=llhT , llhR=llhR , azelT=azelT[r,] , azelR=azelR[r,] , fwhmT=1 , fwhmR=1 , phArrT=FALSE , phArrR=FALSE , freq.Hz=freqT  )
            Btmp           <- igrf(date=date[1:3],lat=latitude[r],lon=longitude[r],height=height[r],isv=0,itype=1)
            B[r,]          <- c(Btmp$y,Btmp$x,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards east,
            dimnames(covar[[r]])   <- list(c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'),c('Ne','Tipar','Tiperp','Tepar','Teperp','Coll','Vix','Viy','Viz',paste('Ion',seq(3),sep=''),'Site1'))
            
            # if there are actual data, all initializations will be updated
            
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
                azelT.gate     <- azelTf[gateinds,]
                azelR.gate     <- azelRf[gateinds,]
                range[r]       <- mean( range( ran.gate ) )
                azelT[r,]      <- colMeans(azelT.gate[!nainds,])
                azelR[r,]      <- colMeans(azelR.gate[!nainds,])

                # average ACF with identical lag values
                lag.unique     <- unique(lag.gate)
                nlag           <- length(lag.unique)
                acf.unique     <- rep(0+0i,nlag)
                var.unique     <- rep(0,nlag)
                for(l in seq(nlag)){
                    lagind <- which(lag.gate==lag.unique[l])
                    for( lind in lagind ){
                        acf.unique[l] <- acf.unique[l] + acf.gate[lind]/var.gate[lind]
                        var.unique[l] <- var.unique[l] + 1/var.gate[lind]
                    }
                }
                var.unique <- 1/var.unique
                acf.unique <- acf.unique * var.unique

                # exact coordinates of the measurement volume
                llhTarget      <- range2llh( azelT=azelT[r,] , llhR=llhR ,  llhT=llhT ,  r=range[r] * 1000 )
                height[r]      <- llhTarget['h'] / 1000
                latitude[r]    <- llhTarget['lat']
                longitude[r]   <- llhTarget['lon']

                # beam intersection
                intersect[[r]] <- beamIntersection( llhT=llhT , llhR=llhR , azelT=azelT[r,] , azelR=azelR[r,] , fwhmT=1 , fwhmR=1 , phArrT=FALSE , phArrR=FALSE , freq.Hz=freqT  )

            
                gainR <- gategain( intersect[[r]] , rlims[r:(r+1)]*1000)
                acf.unique <- acf.unique / gainR
                var.unique <- var.unique / gainR**2

                # magnetic field direction
                Btmp           <- igrf(date=date[1:3],lat=latitude[r],lon=longitude[r],height=height[r],isv=0,itype=1)
                B[r,]          <- c(Btmp$y,Btmp$x,-Btmp$z) # the model has y-axis to east and z-axis downwards, we have x towards east,

                # parameters from iri model
                ptmp           <- iriParams( time=date ,latitude=latitude[r],longitude=longitude[r],heights=height[r])

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
                    acf             = acf.unique,
                    var             = var.unique,
                    nData           = nlag,
                    fSite           = freqT,
                    aSite           = scattVector.llhazelr( llhT , azelT[r,] , llhR , range[r]*1000 , freqT )[["phi"]],
                    kSite           = list(scattVector.llhazelr( llhT , azelT[r,] , llhR , range[r]*1000 , freqT )[["k.ENU"]]),
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
                    lags            = lag.unique,
                    plotTest        = plotTest,
                    plotFit         = plotFit,
                    maxLambda       = maxLambda,
                    maxIter         = maxIter
                    )
                
                # scale back to physical units
                param[r,] <- scaleParams( fitpar$param , scale=parScales , inverse=TRUE )
                covar[[r]] <- scaleCovar( fitpar$covar , scale=parScales , inverse=TRUE)
                std[r,]   <- sqrt(diag(covar[[r]]))
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
        dimnames(azelT) <- list(paste('gate',seq(nr),sep=''),c('az','el'))
        dimnames(azelR) <- list(paste('gate',seq(nr),sep=''),c('az','el'))
        
        
        # save the results to file
        PP <- list(param=param,std=std,model=model,chisqr=chisqr,status=status,time_sec=time_sec,date=date,POSIXtime=POSIXtime,range=range,height=height,latitude=latitude,longitude=longitude,azelT=azelT,llhT=llhT,llhR=llhR,acfScale=acfScale,intersect=intersect,covar=covar,B=B)
        resFile <- file.path( odir , paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "PP.Rdata" , sep=''))
        save( PP , PPI_param, file=resFile )
        
        cat(iperLimits[k+1],'\n')       

        
    }

  }


}
