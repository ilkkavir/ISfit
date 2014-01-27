faradayRotation.KAIRA <- function( ACFxx , ACFyy, ACFxy , ACFyx , VARxx , VARyy, VARxy , VARyx, r , azelT=c(0,90) , llhR=radarSites()[["KIL"]] , llhT=radarSites()[["TRO"]] , azXpol=318 , azelBore=c(0,90) , h=-1)
    {
        #
        # Faraday rotation and autocovariance function estimation from
        # KAIRA cross-covariance measurements
        #
        #
        #
        #
        # I. Virtanen 2013
        #

        # stop if the ACF vectors have different lengths
        if( length( unique( c( length(ACFxx) , length(ACFyy) , length(ACFxy) , length(ACFyx) ) ) ) != 1 ) stop("ACF vectors must be of equal lengths")


        # transmitter location in EFEC coordinates
        xyzT <- LatLonH2xyz( llhT )

        # receiver location in EFEC coordinates
        xyzR <- LatLonH2xyz( llhR )

        # scattering volume location in EFEC coordinates
        xyzS <- LatLonH2xyz( range2llh(  r=r , llhT=llhT , azelT=azelT , llhR=llhR  ) )

        # unit vector pointing from the scatterer to the receiver
        xyzSR <- xyzR - xyzS
        xyzSR <- xyzSR / sqrt(sum(xyzSR**2))

        # unit vector pointing from the scatterer to the transmitter
        xyzST <- xyzT - xyzS
        xyzST <- xyzST / sqrt(sum(xyzST**2))

        # a right-handed cartesian coordinate system with
        # z-axis in the propagation direction of the scattered wave,
        # x-axis perpendicular to both incident and scattered waves, and
        # y-axis in the plane of both incident and scattered waves
        z <- xyzSR
        x <- radarPointings:::normalUnitVector.cartesian( z , xyzST )
        y <- radarPointings:::normalUnitVector.cartesian( z , x )

        # bore-sight direction in EFEC coordinates (zr-axis)
        zr <- ENU2EFEC.xyz( enuVec=azel2ENU( azelBore ) , xyzPos=xyzR )

        # x-polarization direction
        xp <- ENU2EFEC.xyz( enuVec=azel2ENU( c( azXpol , 0 ) ) , xyzPos=xyzR )

        # y-polarization direction (yr-axis), assuming that polarizations are orthogonal
        # with each other and the bore-sight direction
        yr <- radarPointings:::normalUnitVector.cartesian( zr , xp )

        # xr-axis
        xr <- radarPointings:::normalUnitVector.cartesian( yr , zr )


        # angles between x, y, xr, and yr vectors
        phixx <- radarPointings:::vectorAngle.cartesian( x , xr , degrees=FALSE)
        phixy <- radarPointings:::vectorAngle.cartesian( x , yr , degrees=FALSE)
        phiyx <- radarPointings:::vectorAngle.cartesian( y , xr , degrees=FALSE)
        phiyy <- radarPointings:::vectorAngle.cartesian( y , yr , degrees=FALSE)


        # angle between incident and scattered wave vectors
        beta <- radarPointings:::vectorAngle.cartesian( -xyzST , xyzSR , degrees=FALSE )

        # tabulate cosines of phi's
        cosxx <- abs(cos(phixx))
        cosxy <- abs(cos(phixy))
        cosyx <- abs(cos(phiyx))
        cosyy <- abs(cos(phiyy))

        
        # use only those lags from which we have all cross-correlations
        okinds <- which( !( is.na(ACFxx) & is.na(ACFyy) & is.na(ACFxy) & is.na(ACFyx) & is.na(VARxx) & is.na(VARyy) & is.na(VARxy) & is.na(VARyx) ) )

        # return only NA values if useful data was not given
        if( length(okinds)==0 ){
            acf <- rep(NA,length(ACFxx))
            rotAngle <- NA
            return(list(acf=acf,rotAngle=NA))
        }

        # number of usable data points in each ACF vector
        nlags <- length(okinds)

        # combine the ACF and variance vectors
        ACF <- c( ACFxx[okinds] , ACFyy[okinds] , ACFxy[okinds] , ACFyx[okinds] )
        VAR <- c( VARxx[okinds] , VARyy[okinds] , VARxy[okinds] , VARyx[okinds] )


        # initial values for ACF
        acfinit <- ACF[1:nlags] + ACF[(nlags+1):(2*nlags)]

        # initialize faraday rotation to zero, ACF is divided into real and imaginary parts
        p <- c( Re(acfinit) , Im(acfinit) , 0 ,  0)
        
        # simple grid search for iteration start values (with fixed ACF)
        ss <- matrix(ncol=10,nrow=10)
        eccss <- seq(0,.99,length.out=10)
        phiss <- seq(0,pi,length.out=10)
        for(k in seq(10) ){
            for(l in seq(10) ){
                pss <- p
                pss[2*nlags+1] <- eccss[k]
                pss[2*nlags+2] <- phiss[l]
                ss[k,l] <- sum(faradayDirectTheory(pss,ACF,VAR,nlags,cosxx,cosxy,cosyx,cosyy,h=h))
            }
        }

        ssinds <- which(ss==min(ss),arr.ind=TRUE)
#        layout(matrix(c(1,2),ncol=2))
#        image(ss,main=c(eccss[ssinds[1]],phiss[ssinds[2]]),col=rev(gray(seq(1000)/1000)))
#        plot(Re(ACF),type='l')
#        lines(Im(ACF),col='red')
#        lines(sqrt(VAR),col='blue')

        # initialize the eccentricity and faraday rotation with the minimum ss point
        # from grid search
        p[2*nlags+1] <- eccss[ssinds[1]]
        p[2*nlags+2] <- phiss[ssinds[2]]

        p2 <-modFit( f=faradayDirectTheory , p=p , ACF=ACF , VAR=VAR , nlags=nlags ,
                    cosxx=cosxx,cosxy=cosxy,cosyx=cosyx,cosyy=cosyy,
                    h=h,
                    lower=c(rep(-Inf,nlags*2), 0 , 0 ),
                    upper=c(rep(Inf,nlags*2), 1 , pi),
                    control=list(maxiter=500),
                    plot=FALSE
                    )
#cat('\n')        
        covar <- tryCatch( solve(p2$hessian) , error=function(x){print('err');matrix(NA,nrow=(2*nlags+2),ncol=(2*nlags+2))} )

        reslist <- list()
        reslist$ACF <- p2$par[1:nlags]+1i*p2$par[(nlags+1):(2*nlags)]
        reslist$var <- rep(0,nlags)
        for(k in seq(nlags)) reslist$var[k] <- covar[k,k] + covar[(k+nlags),(k+nlags)] + 2*covar[nlags,(k+nlags)]
        reslist$phi <- p2$par[2*nlags+2]
        reslist$varphi <- covar[(2*nlags+2),(2*nlags+2)]
        reslist$beta <- beta
        reslist$e <- p2$par[2*nlags+1]
        reslist$evar <- covar[(2*nlags+1),(2*nlags+1)]
        reslist$covar <- covar

        return(reslist)
    }
