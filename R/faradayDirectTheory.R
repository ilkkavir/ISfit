faradayDirectTheory <- function( p , ACF , VAR , nlags , cosxx , cosxy , cosyx , cosyy , h=-1 , plot=F)
    {
        # direct theory function for faraday rotation estimation
        # We are apparently unable to fit the known ellipsoid because
        # the x and y polarization are not in phase and have different gains.
        # We will thus simply estimate the polarization ellipsoid based on
        # the measured cross-covariances, ALTHOUGH WE KNOW THAT THE ESTIMATED
        # ELLIPSOID IS NOT THE CORRECT ONE!
        #
        # I. Virtanen 2013
        #

        # ellipsoid semi-major axis for unit power
        A <- 1/sqrt(2-p[2*nlags+1])

        # semi-minor axis
        B <- sqrt(1-A**2)

        # cross-covariance powers
        Pxx <- .5*( A**2 * ( 1 + cos( 2*p[2*nlags+2] ) ) + B**2 * ( 1 - cos( 2*p[2*nlags+2] ) ) )
        Pyy <- .5*( A**2 * ( 1 - cos( 2*p[2*nlags+2] ) ) + B**2 * ( 1 + cos( 2*p[2*nlags+2] ) ) )
        Pxy <- .5*( A**2 - B**2 ) * sin( 2*p[2*nlags+2] ) - 1i * A * B * h
        Pyx <- .5*( A**2 - B**2 ) * sin( 2*p[2*nlags+2] ) + 1i * A * B * h

        Pxxr <- Pxx * cosxx**2 + Pyy * cosyx**2 + ( Pxy + Pyx ) * cosxx * cosyx
        Pyyr <- Pxx * cosxy**2 + Pyy * cosyy**2 + ( Pxy + Pyx ) * cosxy * cosyy
        Pxyr <- Pxx * cosxx * cosxy + Pyy * cosyx * cosyy + Pxy * cosxx * cosyy + Pyx * cosyx * cosxy
        Pyxr <- Pxx * cosxx * cosxy + Pyy * cosyx * cosyy + Pxy * cosxy * cosyx + Pyx * cosyy * cosxx

        acfout  <- c( (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pxxr ,
                      (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pyyr ,
                      (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pxyr ,
                      (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pyxr )
        
        res <- abs( ACF - acfout )**2  / VAR 
#cat('                                           \r',Pxxr,Pyyr,Pxyr,Pyxr,p[2*nlags+1],sum(res)/(2*nlags+2))
#cat(Pxxr,Pyyr,Pxyr,Pyxr,p[2*nlags+1],p[2*nlags+2],sum(res),'\n')
        if(plot){
            layout(matrix(c(1)))
            plot(Re(ACF),ylim=c(-3,3))
            points(Im(ACF),col='red')
            lines(Re(acfout))
            lines(Im(acfout),col='red')
        }

        return( res )

    }
