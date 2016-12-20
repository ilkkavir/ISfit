Nwind <- function(E,Ecov,vi,vicov,h,B,time,lat,lon){
    # neutral wind estimates
    #
    # INPUT:
    #  E     Electric field [V/m] in cartesian coordinates with x-axis to magnetic east and
    #        y-axis to magnetic north, perpendicular to B.
    #  Ecov  Covariance matrix of the electric field.
    #  vi    a nh x 3 matrix of plasma velocity estimates [m/s]. In geographic ENU coordinates.
    #  vicov covariance matrices of the velocity estimates (an nh x 3 x 3 array)
    #  h     heights [km]
    #  B     magnetic fields in each height gate [nT] (nh x 3 array)
    #  time  time in POSIX format
    #  lat   latitudes of each height gate
    #  lon   longitudes of each height gate
    #
    # OUTPUT:
    #  a list with elements
    #    nWind neutral wind estimates in geodetic ENU coordinates
    #    cov   covariance matrices for each 3D estimate
    #    nWindB neutral wind in geomagnetic coordinates
    #    covB   covariance matrices for each 3D estimate in geomagnetic coordinates
    #
    # I. Virtanen 2016


    date <- c(time[["year"]]+1900,time[["mon"]]+1,time[["mday"]],time[["hour"]],time[["min"]],time[["sec"]])
    echarge <- 1.60217662e-19
    mion <- c( 'i1'=30.5 , 'i2'=16 )
    u <- 1.66053904e-27

    # pad E and Ecov with zeros
    E2 <- c(E,0)
    Ecov2 <- matrix(0,ncol=3,nrow=3)
    Ecov2[1:2,1:2] <- Ecov

    nh <- length(h)
    nWind <- matrix(nrow=nh,ncol=3)
    nWindCov <- array(dim=c(nh,3,3))
    nWindB <- nWind
    nWindCovB <- nWindCov

    for(hh in seq(nh)){
        # the coordinate system will be different in each gate

        # horizontal component of B
        Bhor <- c(B[hh,c(1,2)],0)
        names(Bhor) <- c('x','y','z')
        # geomagnetic east
        Bx <- radarPointings:::vectorProduct.cartesian(B[hh,],Bhor)
        Bx <- Bx/sqrt(sum(Bx**2))
        names(Bx) <- c('x','y','z')
        # geomagnetic north, perpendicular to B
        By <- radarPointings:::vectorProduct.cartesian(Bx,B[hh,])
        By <- By/sqrt(sum(By**2))
        names(By) <- c('x','y','z')
        # up along B
        Bz <- -B[hh,]/sqrt(sum(B[hh,]**2))
        names(Bz) <- c('x','y','z')
        # field intensity (the input  is in nT)
        Bstrength <- sqrt(sum(B[hh,]**2)) * 1e-9
        names(Bstrength) <- NULL
        # Ion velocity projection to Bx, By, Bz coordinates
        rotMat <- matrix(c(Bx,By,Bz),ncol=3,byrow=TRUE)
        ViBxyz <- rotMat%*%vi[hh,]
        names(ViBxyz) <- NULL
        # Covariance
        ViBxyzCov <- rotMat%*%vicov[[hh]]%*%t(rotMat)
        dimnames(ViBxyzCov) <- NULL

        # compositions etc.
        ptmp <- iriParams( time=date ,latitude=lat[hh],longitude=lon[hh],heights=h[hh])

        # better collision frequencies and two ions
        ioncoll <- ionNeutralCollisionFrequency2(ptmp[,1])

        # gyrofrequencies for both ion species
        iongyro <- echarge * Bstrength / ( mion * u )
        # pedersen mobilities, use the element names to make sure that everything is in correct order
        kp1 <- ( iongyro[c('i1','i2')] * ioncoll[c('i1','i2')] ) /
            ( echarge * Bstrength * ( iongyro[c('i2','i2')]**2 + ioncoll[c('i1','i2')]**2 ) )
        # hall mobilities
        kh1 <- ( iongyro[c('i1','i2')]**2 ) /
            ( echarge * Bstrength * ( iongyro[c('i1','i2')]**2 + ioncoll[c('i1','i2')]**2 ) )
        # parallel mobilities
        kb1 <- 1 / ( mion[c('i1','i2')] * u * ioncoll[c('i1','i2')] )

        # fractions of i1 and i2 from IRI model
        comp <- ptmp['O+',1] / sum(ptmp[c('NO+','O2+','O+'),1])
        comp <- c( (1-comp) , comp )
        names(comp) <- c('i1','i2')
        comp <- comp / sum(comp)

        # total conductivities
        kp <- sum( comp * kp1 )
        kh <- sum( comp * kh1 )
        kb <- sum( comp * kb1)

        # inverse of ion mobility tensor, this is differenct from some
        # textbooks since we have z axis upwards
        mobi <- matrix(0,ncol=3,nrow=3)
        mobi[1,1] <- kp/(kh**2+kp**2)
        mobi[2,1] <- -kh/(kh**2+kp**2)
        mobi[1,2] <- -mobi[2,1]
        mobi[2,2] <- mobi[1,1]
        mobi[3,3] <- 1/kb

        # ViBxyz = mob%*%(echarge*E+mion*u*ioncoll*nWind)
        # multiply with mobi from left
        Vimob <- mobi%*%ViBxyz
        VimobCov <- mobi%*%ViBxyzCov%*%t(mobi)

        # subtract the electric field contribution and divide by mion*u*ioncoll
        nWindB[hh,] <- ( Vimob - echarge*E2 ) / sum( comp * mion[c('i1','i2')] * u * ioncoll[c('i1','i2')])
        nWindCovB[hh,,] <- ( VimobCov + echarge**2*Ecov2 ) /
            sum( comp * mion[c('i1','i2')] * u * ioncoll[c('i1','i2')] )**2

        # rotate to ENU coordinates
        rotMati <- solve(rotMat)
        nWind[hh,] <- rotMati%*%nWindB[hh,]
        nWindCov[hh,,] <- rotMati%*%nWindCovB[hh,,]%*%t(rotMati)

    }


    return(list(nWind=nWind,cov=nWindCov,nWindB=nWindB,covB=nWindCovB))

}
