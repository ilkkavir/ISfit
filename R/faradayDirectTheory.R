faradayDirectTheory <- function( p , nlags , cosxx , cosxy , cosyx , cosyy , beta , h=-1 , phcorr=0 , syy=1 , scross=1, ... )
    {
        #
        # direct theory function for faraday rotation estimation
        #
        # This function is for fitting a known ellipsoid when the
        # phase correction phcorr, the y-polarization scale syy, and
        # the cross-polarization scale scross are known. 
        #
        # Input:
        #   p      parameters, p[1:nlags] is real part of ACF,
        #                  p[(nlags+1):(2*nlags)] is imaginary part of ACF,
        #                  p[2*nlags+1] is ellipsoid rotation angle
        #   nlags  number of lags in ACf
        #   cosxx, cosxy, cosyy, cosines of angles between coordinate axes of the
        #                 system in which the wave propagation direction is positive
        #                 z-axis, and projections of x- and y-dipole directions to this
        #                 xy-plane
        #   beta   angle between incident and scattered waves
        #   h      ellipse handedness
        #   phcorr correction to relative phase in between the two receiver polarizations
        #   syy    correction to y-polarization gain
        #   scross correction to cross-polarization gain
        #
        # Returns:
        #
        # acfout  a vector of ACFs and CCFS. acfout[1:nlags] is the x-polarization ACF
        #         acfout[(nlags+1):(2*nlags)] is the y-polarization ACF
        #         acfout[(2*nlags+1):(3*nlags)] is the xy CCF
        #         acfout[(3*nlags+1):(4*nlags)] is the yx CCF
        #
        # I. Virtanen 2013, 2014
        #

        cossqr <- cos(beta)**2
        
        # ellipsoid semi-major axis for unit power
        A <- 1/sqrt(1+cossqr)

        # semi-minor axis
        B <- sqrt(cossqr/(1+cossqr))

        # cross-covariance powers
        Pxx <- .5*( A**2 * ( 1 + cos( 2*p[2*nlags+1] ) ) + B**2 * ( 1 - cos( 2*p[2*nlags+1] ) ) )
        Pyy <- .5*( A**2 * ( 1 - cos( 2*p[2*nlags+1] ) ) + B**2 * ( 1 + cos( 2*p[2*nlags+1] ) ) )
        Pxy <- .5*( A**2 - B**2 ) * sin( 2*p[2*nlags+1] ) - 1i * A * B * h
        Pyx <- .5*( A**2 - B**2 ) * sin( 2*p[2*nlags+1] ) + 1i * A * B * h

        # cross-covariance powers in dipoles
        Pxxr <- Pxx * cosxx**2 + Pyy * cosyx**2 + ( Pxy + Pyx ) * cosxx * cosyx
        Pyyr <- Pxx * cosxy**2 + Pyy * cosyy**2 + ( Pxy + Pyx ) * cosxy * cosyy
        Pxyr <- Pxx * cosxx * cosxy + Pyy * cosyx * cosyy + Pxy * cosxx * cosyy + Pyx * cosyx * cosxy
        Pyxr <- Pxx * cosxx * cosxy + Pyy * cosyx * cosyy + Pxy * cosxy * cosyx + Pyx * cosyy * cosxx

        # phase and amplitude corrections
        Pyyr <- Pyyr * syy
        Pxyr <- Pxyr * scross * exp(1i*phcorr)
        Pyxr <- Pyxr * scross * exp(-1i*phcorr)

        # the real part of acf is in p[1:nlags] and the imaginary part in p[(nlags+1):(2*nlags)]
        # the output vector contains is in order ACFxx, ACFyy, CCFxy, CCFyx
        acfout  <- c( (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pxxr ,
                      (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pyyr ,
                      (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pxyr ,
                      (p[1:nlags] + 1i*p[(nlags+1):(2*nlags)]) * Pyxr )

        return( acfout )

    }
