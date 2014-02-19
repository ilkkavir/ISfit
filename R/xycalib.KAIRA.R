xycalib.KAIRA  <- function( dpath , beginTime=c(2012,1,1,0,0,0),endTime=c(2100,1,1,0,0,0) , plotFit=FALSE )
    {
        #
        # Calibration of KAIRA gains and phases of the two linear receiver polarizations of KAIRA
        # The input data are general polarization ellipses from ellipseFit.KAIRA, and the output
        # is estimates of phcorr, syy, and scross for faradayFit.KAIRA
        #
        #
        # dpath           yx-ccfs
        # beginTime        analysis start time c(year,month,day,hour,min,sec)
        # endTime          analysis end time
        # azXpol           x-polarizatio dipole azimuth angle in degrees (318 for KAIRA)
        # azelBore         bore sight direction (c(0,90) for KAIRA)
        # h                handedness of the transmitted circular polarization (-1 for EISCAT)
        #
        # 

        # list all data files
        dfiles <- dir( dpath , full.names=T , pattern="[[:digit:]]ELP.Rdata")

        # stop if there is no data
        if ( sum( sapply( dfiles , length ) ) == 0 ) stop( "Could not find any data files." )

        # read timestamps from file names. This is the unix time at end of integration period
        tstamps <- lapply( dfiles , function(x){as.numeric(substr(x,nchar(x)-21,nchar(x)-9)) / 1000 } )

        # convert beginTime and endTime into unix time format
        bTime <- as.double(ISOdate(beginTime[1],beginTime[2], beginTime[3],beginTime[4], beginTime[5],beginTime[6])) + beginTime[6]%%1
        eTime <- as.double(ISOdate(endTime[1],endTime[2], endTime[3],endTime[4], endTime[5],endTime[6])) + endTime[6]%%1


        # use only those data files that are between bTime and eTime
        dInds <- which( (tstamps >= bTime) & (tstamps <= eTime))

        dfiles <- dfiles[ dInds ]


        # read scattering angles and rotation angles from files
        nfiles <- length(dfiles)
        if( nfiles == 0 ) stop("No data")

        load( dfiles[1] )
        nr         <- length( ACF[["range.km"]] )
        phiest     <- matrix( nrow=nr , ncol=nfiles )
        phiestvar  <- matrix( nrow=nr , ncol=nfiles )
        betaest    <- matrix( nrow=nr , ncol=nfiles )
        betaestvar <- matrix( nrow=nr , ncol=nfiles )

        beta <- acos(sqrt(cos(ACF[["betaGeom"]])**2))
        cosxx <- ACF[["cosxx"]]
        cosxy <- ACF[["cosxy"]]
        cosyx <- ACF[["cosyx"]]
        cosyy <- ACF[["cosyy"]]
        handedness <- ACF[["handedness"]]
        for( k in seq( nfiles ) ){
            load( dfiles[k] )
            phiest[,k]     <- ACF[["phi"]]%%pi
            phiestvar[,k]  <- ACF[["varphi"]]
            betaest[,k]    <- acos(sqrt(cos(ACF[["beta"]])**2))
            betaestvar[,k] <- ACF[["varbeta"]]
        }

        phcorr <- syy <- scross <- rep( NA , nr )
        phcorrv <- syyv <- scrossv <- rep( NA , nr )

        for( r in seq( nr ) ){
            dinds <- which( !is.na(phiest[r,]) &!is.na(betaest[r,]) &!is.na(betaestvar[r,]) )
            if( length(dinds)>0 ){
                plot(betaest[r,])
                abline(h=beta[r],col='red')
                Sys.sleep(10)
                # try grid search first
                phc <- seq(-pi,pi,length.out=10)
                yyc <- seq(1/10,1/5,1/2,1,2,5,10)
                xyc <- seq(1/10,1/5,1/2,1,2,5,10)
                nph <- length(phc)
                nyy <- length(yyc)
                nxy <- length(xyc)
                ssm <- array(dim=c(nph,nyy,nxy))
                for(k in seq(nph)){
                    for(l in seq(nyy)){
                        for(m in seq(nxy)){
                            ssm[k,l,m] <- sum((xycalibDirectTheory(p=c(phc[k],yyc[l],xyc[m]),betaest[r,dinds],phiest[r,dinds],cosxx[r],cosxy[r],cosyx[r],cosyy[r],handedness)-beta[r])**2/betaestvar[r,dinds])
                        }
                    }
                }
                ssind <- which(ssm==min(ssm),arr.ind=TRUE)
print(                p <- c( phc[ssind[1]] , yyc[ssind[2]] , xyc[ssind[3]] ))
                p2 <- leastSquare.lvmrq( measData=rep( beta[r] , length(dinds) ) , measVar=betaestvar[r,dinds] , initParam=p , aprioriTheory=matrix(rep(0,3),nrow=1) , aprioriMeas=0 , invAprioriCovar=1, paramLimits=matrix(c( -pi , .1 , .1 , pi , 10 , 10 ),nrow=2,byrow=TRUE),directTheory=xycalibDirectTheory,beta=betaest[r,dinds],phi=phiest[r,dinds],cosxx=cosxx[r],cosxy=cosxy[r],cosyx=cosyx[r],cosyy=cosyy[r],h=handedness,plotFit=plotFit)
                phcorr[r]  <- p2[["param"]][1]
                syy[r]     <- p2[["param"]][2]
                scross[r]  <- p2[["param"]][3]
                phcorrv[r] <- p2[["covar"]][1,1]
                syyv[r]    <- p2[["covar"]][2,2]
                scrossv[r] <- p2[["covar"]][3,3]
            }
        }

        return( list( phcorr=phcorr , syy=syy , scross=scross , phcorrv=phcorrv , syyv=syyv , scrossv=scrossv ) )
        
    }
