polEllipse.KAIRA <- function( ACFxx , ACFyy, ACFxy , ACFyx , VARxx , VARyy, VARxy , VARyx, r , azelT=c(0,90) , llhR=radarSites()[["KIL"]] , llhT=radarSites()[["TRO"]] , azXpol=318 , azelBore=c(0,90) , h=-1 , plotFit=FALSE , absLimit=2)
    {
        #
        # 
        # Fit general polarization ellipses to KAIRA data
        #
        #
        #
        #
        # I. Virtanen 2013, 2014
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

        # projections of xr and yr to the plane of wave electric field
        xrp <- xr - sum( z*xr/sqrt(sum(z**2)))*z
        yrp <- yr - sum( z*yr/sqrt(sum(z**2)))*z

        # angles between x, y, xrp, and yrp vectors
        phixx <- radarPointings:::vectorAngle.cartesian( x , xrp , degrees=FALSE)
        phixy <- radarPointings:::vectorAngle.cartesian( x , yrp , degrees=FALSE)
        phiyx <- radarPointings:::vectorAngle.cartesian( y , xrp , degrees=FALSE)
        phiyy <- radarPointings:::vectorAngle.cartesian( y , yrp , degrees=FALSE)

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

        # initialize faraday rotation and scattering angles to zero and beta,
        # ACF is divided into real and imaginary parts
        p <- c( Re(acfinit) , Im(acfinit) , 0 , beta )
        
        p2 <- leastSquare.lvmrq( measData=ACF , measVar=VAR , initParam=p , aprioriTheory=matrix(rep(0,2*nlags+2),nrow=1) , aprioriMeas=0 , invAprioriCovar=1, paramLimits=matrix(c(rep(-Inf,nlags*2), -2*pi , -2*pi , rep(Inf,nlags*2), 2*pi , 2*pi ),nrow=2,byrow=TRUE),directTheory=ellipseDirectTheory,nlags=nlags,cosxx=cosxx,cosxy=cosxy,cosyx=cosyx,cosyy=cosyy,h=h,plotFit=plotFit,absLimit=absLimit)

        reslist <- list()
        reslist$ACF <- p2$param[1:nlags]+1i*p2$param[(nlags+1):(2*nlags)]
        reslist$var <- rep(0,nlags)
        for(k in seq(nlags)) reslist$var[k] <- p2$covar[k,k] + p2$covar[(k+nlags),(k+nlags)]
        reslist$phi <- p2$param[2*nlags+1]
        reslist$varphi <- p2$covar[(2*nlags+1),(2*nlags+1)]
        reslist$beta <- p2$param[2*nlags+2]
        reslist$varbeta <- p2$covar[(2*nlags+2),(2*nlags+2)]
        reslist$betaGeom <- beta
        reslist$covar <- p2$covar
        reslist$chisqr <- p2$chisqr
        reslist$cosxx <- cosxx
        reslist$cosxy <- cosxy
        reslist$cosyx <- cosyx
        reslist$cosyy <- cosyy

        return(reslist)
    }
