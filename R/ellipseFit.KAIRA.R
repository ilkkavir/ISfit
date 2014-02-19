ellipseFit.KAIRA  <- function( xxpath , yypath , xypath , yxpath , opath, heightLimits.km , beginTime=c(2012,1,1,0,0,0),endTime=c(2100,1,1,0,0,0) , timeRes.s , maxdev=1, azXpol=318 , azelBore=c(0,90) , handedness=-1 , plotFit=FALSE , absLimit=2)
    {
        #
        # fit general polarization ellipses to KAIRA data without gain and phase corrections
        #
        #
        # xxpath           path to x-polarization ACF results
        # yypath           y-pol acfs
        # xypath           xy-ccfs
        # yxpath           yx-ccfs
        # opath            output path
        # heightLimits.km  height gate limits in km
        # beginTime        analysis start time c(year,month,day,hour,min,sec)
        # endTime          analysis end time
        # timeRes.s        analysis time resolution in seconds, or a vector of
        #                  integration period limits
        # azXpol           x-polarizatio dipole azimuth angle in degrees (318 for KAIRA)
        # azelBore         bore sight direction (c(0,90) for KAIRA)
        # h                handedness of the transmitted circular polarization (-1 for EISCAT)
        #
        # 

        # create the output directory
        dir.create( opath , recursive=TRUE , showWarnings=FALSE)
      
        # list all data files
        dfiles <- list()
        dfiles[['xx']] <- dir( xxpath , full.names=T , pattern="[[:digit:]]LP.Rdata")
        dfiles[['yy']] <- dir( yypath , full.names=T , pattern="[[:digit:]]LP.Rdata")
        dfiles[['xy']] <- dir( xypath , full.names=T , pattern="[[:digit:]]LP.Rdata")
        dfiles[['yx']] <- dir( yxpath , full.names=T , pattern="[[:digit:]]LP.Rdata")

        # stop if there is no data
        if ( sum( sapply( dfiles , length ) ) == 0 ) stop( "Could not find any data files." )

        # read timestamps from file names. This is the unix time at end of integration period
        tstamps <- lapply( dfiles , function(x){as.numeric(substr(x,nchar(x)-20,nchar(x)-8)) / 1000 } )

        # convert beginTime and endTime into unix time format
        bTime <- as.double(ISOdate(beginTime[1],beginTime[2], beginTime[3],beginTime[4], beginTime[5],beginTime[6])) + beginTime[6]%%1
        eTime <- as.double(ISOdate(endTime[1],endTime[2], endTime[3],endTime[4], endTime[5],endTime[6])) + endTime[6]%%1
        
        if(length(timeRes.s)==1){
            # first integration period from which we actually have data
            iperFirst <- max( 1, floor( ( max( sapply( tstamps , min ) ) - bTime ) / timeRes.s ) )
    
            # last integration period
            iperLast <- ceiling( ( min( min( sapply( tstamps , max ) ) , eTime ) - bTime ) / timeRes.s ) 

            # stop if analysis end time is before its start time
            if( iperLast < iperFirst ) stop()
    
            # integration period limits
            iperLimits <- seq( iperFirst , iperLast ) * timeRes.s + bTime

        }else{
            iperLimits <- timeRes.s
        }

        # number of integration periods
        nIper <- length(iperLimits) - 1
      
        # walk through all integration periods
        for( k in seq( nIper ) ){

            # look for data files from the current integration period
            iperFiles <- lapply( tstamps , function(x,l1,l2){ which( ( x > l1 ) & ( x <= l2 ) ) } , l1=iperLimits[k] , l2=iperLimits[k+1] )

            # load all data and collect it in a list
            nFiles <- sum( sapply( iperFiles , length ) )
            if( nFiles > 0 ){
                
                dlist <- list()

                for( n in c('xx','yy','xy','yx') ){
                    dlist[[n]] <- readACF( dfiles[[n]][ iperFiles[[n]] ] )
                    if(length(dlist[[n]][["nGates"]])!=length(dlist[[n]][["lag.us"]])){
                        dlist[[n]][["nGates"]] <- rep(length(dlist[[n]][["range.km"]]),length(dlist[[n]][["lag.us"]]))
                    }
                }
                
                # a time vector converted from iperLimits
                t <- as.POSIXlt( iperLimits[k+1] , origin='1970-01-01' , tz='ut')
                date <- c(t$year+1900,t$mon,t$mday,t$hour,t$min,t$sec)
              
                # height gate limits
                if( all( is.na( heightLimits.km ) ) ){
                    rlims <- sort( unique( dlist[["xx"]][["range.km"]] ) )*1000
                    rlims <- c( rlims , max(rlims) + 1)
                    hlims <- rlims
                    for( h in seq(length(hlims))) hlims[h] <- range2llh(r=rlims[h],llhT=dlist[["xx"]][["llhT"]],llhR=dlist[["xx"]][["llhR"]],azelT=dlist[["xx"]][["azelT"]])['h']
                }else{
                    hlims <- unique( heightLimits.km )*1000
                }
                nh <- length( hlims ) - 1
                
                # convert all ranges to latitude, longitude, height
                llh <- matrix(nrow=length(dlist[["xx"]][["range.km"]]),ncol=3)
                for( dind in seq(length(dlist[["xx"]][["range.km"]]))){
                    llh[dind,] <- range2llh( r=dlist[["xx"]][["range.km"]][dind]*1000 , llhT=dlist[["xx"]][["llhT"]] , llhR=dlist[["xx"]][["llhR"]] , azelT=dlist[["xx"]][["azelT"]])
                }

                ACFout <- matrix(NA*(0+0i),nrow=nh,ncol=length(dlist[["xx"]][["lag.us"]]))
                VARout <- matrix(nrow=nh,ncol=length(dlist[["xx"]][["lag.us"]]))
                cosxx <- cosxy <- cosyx <- cosyy <- betaGeom <- beta <- varbeta <- latitude <- longitude <- height <- ran.gate <- phi <- varphi <- chisqr <- rep(NA,nh)

                intersect <- list()
                
                for( h in seq( nh ) ){
                  
                    # Initial values, these will be immediately updated if any data is found
                    height[h]      <- sum(hlims[h:(h+1)])/2000
                  
                    # if there are actual data, all initializations will be updated
                  
                    gateinds       <- which( ( llh[,3] >= hlims[h] ) & (llh[,3] < hlims[h+1]) )

                    acf.gate       <- list()
                    var.gate       <- list()
                    lag.gate       <- dlist[["xx"]][["lag.us"]]*1e-6
                    ran.gate[h]    <- mean(dlist[["xx"]][["range.km"]][gateinds]*1000)

                    if(length(gateinds)==1){
                        llh.gate <- llh[gateinds,]
                    }else{
                        llh.gate <- colMeans(llh[gateinds,])
                    }
                    for(n in c("xx","yy","xy","yx")){
                        if(length(gateinds)==1){
                            acf.gate[[n]] <- dlist[[n]][["ACF"]][ gateinds , ]
                            var.gate[[n]] <- dlist[[n]][["var"]][ gateinds , ]
                        }else{
                            acf.gate[[n]] <- colSums( dlist[[n]][["ACF"]][gateinds,] /
                                                     dlist[[n]][["var"]][gateinds,])
                            var.gate[[n]] <- 1/colSums( 1/dlist[[n]][["var"]][gateinds,])
                            acf.gate[[n]] <- acf.gate[[n]] * var.gate[[n]]
                        }
                    }

                  
                    # coordinates of the measurement volume
                    llhTarget      <- llh.gate
                    height[h]      <- llhTarget[3] / 1000
                    latitude[h]    <- llhTarget[1]
                    longitude[h]   <- llhTarget[2]

                      
                    # beam intersections
                    # the beam widths and antenna types are stored in dscales
                    intersect[[h]] <- beamIntersection( llhT=dlist[["xx"]][["llhT"]] , llhR=dlist[["xx"]][["llhR"]] , azelT=dlist[["xx"]][["azelT"]] , azelR=dlist[["xx"]][["azelR"]] , fwhmT=1.3 , fwhmR=2 , phArrT=F , phArrR=T , freq.Hz=dlist[["xx"]][["radarFreq"]] )
                          
                    # conversion from lat, lon, height to range in this gate
                    rs1 <- height2range( llhT=dlist[["xx"]][["llhT"]] , azelT=dlist[["xx"]][["azelT"]] , llhR=dlist[["xx"]][["llhR"]] , h=hlims[h] )
                    rs2 <- height2range( llhT=dlist[["xx"]][["llhT"]] , azelT=dlist[["xx"]][["azelT"]] , llhR=dlist[["xx"]][["llhR"]] , h=hlims[h+1] )

                    # gain integral
                    gainR <- gategain( intersect[[h]] , c(rs1,rs2) , maxdev=maxdev)


                    if(!is.na(gainR) & !is.na(ran.gate[h])){

                        acfscale <- max(abs(acf.gate[["xx"]]),na.rm=TRUE)
                        frot <- polEllipse.KAIRA( ACFxx=acf.gate[["xx"]]/acfscale,
                                                 ACFyy=acf.gate[["yy"]]/acfscale,
                                                 ACFxy=acf.gate[["xy"]]/acfscale,
                                                 ACFyx=acf.gate[["yx"]]/acfscale,
                                                 VARxx=var.gate[["xx"]]/acfscale**2,
                                                 VARyy=var.gate[["yy"]]/acfscale**2,
                                                 VARxy=var.gate[["xy"]]/acfscale**2,
                                                 VARyx=var.gate[["yx"]]/acfscale**2,
                                                 r=ran.gate[h],
                                                 azelT=dlist[["xx"]][["azelT"]],
                                                 llhR=dlist[["xx"]][["llhR"]],
                                                 llhT=dlist[["xx"]][["llhT"]],
                                                 azXpol=azXpol,
                                                 azelBore=azelBore,
                                                 h=handedness,
                                                 plotFit=plotFit,
                                                 absLimit=absLimit
                                                 )
                        ACFout[h,] <- frot[["ACF"]]*acfscale
                        VARout[h,] <- frot[["var"]]*acfscale**2
                        phi[h] <- frot[["phi"]]
                        varphi[h] <- frot[["varphi"]]
                        chisqr[h] <- frot[["chisqr"]]
                        beta[h] <- frot[["beta"]]
                        varbeta[h] <- frot[["varbeta"]]
                        betaGeom[h] <- frot[["betaGeom"]]
                        cosxx[h] <- frot[["cosxx"]]
                        cosxy[h] <- frot[["cosxy"]]
                        cosyx[h] <- frot[["cosyx"]]
                        cosyy[h] <- frot[["cosyy"]]
                        Sys.sleep(1)
                    }

                }

              
                time_sec <- iperLimits[k+1]
                POSIXtime <- as.POSIXlt(time_sec,origin='1970-01-01',tz='ut')


                ACF <- list()
                ACF[["ACF"]] <- ACFout
                ACF[["var"]] <- VARout
                ACF[["phi"]] <- phi
                ACF[["varphi"]] <- varphi
                ACF[["chisqrFR"]] <- chisqr
                ACF[["functionCall"]] <- match.call(expand.dots=TRUE)
                ACF[["nGates"]] <- nh
                ACF[["backgroundACF"]] <- dlist[["xx"]][["backgroundACF"]]*NA
                ACF[["backgroundvar"]] <- dlist[["xx"]][["backgroundvar"]]*NA
                ACF[["timeString"]] <- dlist[["xx"]][["timeString"]]
                ACF[["time.s"]] <- dlist[["xx"]][["time.s"]]
                ACF[["lag.us"]] <- dlist[["xx"]][["lag.us"]]
                ACF[["range.km"]] <- ran.gate/1000
                ACF[["llhT"]] <- dlist[["xx"]][["llhT"]]
                ACF[["llhR"]] <- dlist[["xx"]][["llhR"]]
                ACF[["azelT"]] <- dlist[["xx"]][["azelT"]]
                ACF[["azelR"]] <- dlist[["xx"]][["azelR"]]
                ACF[["radarFreq"]] <- dlist[["xx"]][["radarFreq"]]
                ACF[["maxRanges"]] <- rep(Inf,length(ACF[["lag.us"]]))
                ACF[["beta"]] <- beta
                ACF[["varbeta"]] <- varbeta
                ACF[["betaGeom"]] <- betaGeom
                ACF[["cosxx"]] <- cosxx
                ACF[["cosxy"]] <- cosxy
                ACF[["cosyx"]] <- cosyx
                ACF[["cosyy"]] <- cosyy
                ACF[["azXpol"]] <- azXpol
                ACF[["azelBore"]] <- azelBore
                ACF[["handedness"]] <- handedness
                ACF[["note"]] <- "ACFs from general polarization ellipse fit, do  not use for plasma parameter estimation"

                resFile <- file.path( opath , paste( sprintf( '%13.0f' , trunc( iperLimits[k+1]  * 1000 ) ) , "ELP.Rdata" , sep=''))
                
                save( ACF , file=resFile )

                cat(iperLimits[k+1],'\n')       

            }

        }

    }
