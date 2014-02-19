xycalibDirectTheory <- function( p , beta , phi , cosxx , cosxy , cosyx , cosyy , h=-1 , ... )
    {
        #
        # direct theory function for estimating gain and phase corrections for KAIRA
        #
        #
        # Input:
        #   p      p[1] is phase correction between x and y dipoles
        #          p[2] is gain correction for y dipoles
        #          p[3] is gain correction for cross covariances
        #   beta   scattering angle (estimated earlier without phase and gain corrections)
        #   phi    ellipse rotation angle (estimated earlier without phase and gain corrections)
        #   cosxx, cosxy, cosyx, cosyy, cosines of angles between coordinate axes of the
        #          system in which the wave propagation direction is positive
        #          z-axis, and projections of x- and y-dipole directions to this
        #          xy-plane
        #   h      ellipse handedness
        #
        # Returns:
        #   beta   angle between incident and scattered waves
        #
        # I. Virtanen 2013, 2014
        #

        cossqr <- cos(beta)**2
        
        # ellipsoid semi-major axis for unit power
        A <- 1/sqrt(1+cossqr)

        # semi-minor axis
        B <- sqrt(cossqr/(1+cossqr))

        # cross-covariance powers
        Pxx <- .5*( A**2 * ( 1 + cos( 2*phi ) ) + B**2 * ( 1 - cos( 2*phi ) ) )
        Pyy <- .5*( A**2 * ( 1 - cos( 2*phi ) ) + B**2 * ( 1 + cos( 2*phi ) ) )
        Pxy <- .5*( A**2 - B**2 ) * sin( 2*phi ) - 1i * A * B * h
        Pyx <- .5*( A**2 - B**2 ) * sin( 2*phi ) + 1i * A * B * h

        # cross-covariance powers in dipoles
        Pxxr <- Pxx * cosxx**2 + Pyy * cosyx**2 + ( Pxy + Pyx ) * cosxx * cosyx
        Pyyr <- Pxx * cosxy**2 + Pyy * cosyy**2 + ( Pxy + Pyx ) * cosxy * cosyy
        Pxyr <- Pxx * cosxx * cosxy + Pyy * cosyx * cosyy + Pxy * cosxx * cosyy + Pyx * cosyx * cosxy
        Pyxr <- Pxx * cosxx * cosxy + Pyy * cosyx * cosyy + Pxy * cosxy * cosyx + Pyx * cosyy * cosxx

        # apply the corrections
        Pyyr <- Pyyr / p[2]
        Pxyr <- Pxyr / p[3] * exp(-1i*p[1])
        Pyxr <- Pyxr / p[3] * exp(1i*p[1])
        
        # back to the plane perpendicular to the wave direction
        # these were solved by Mathematica...
        Pxx <- (cosyy**2 * Pxxr - cosyx * cosyy * (Pxyr + Pyxr) + cosyx**2 * Pyyr ) / (cosxy * cosyx - cosxx * cosyy)**2
        Pyy <- (cosxy**2 * Pxxr - cosxx * cosxy * (Pxyr + Pyxr) + cosxx**2 * Pyyr ) / (cosxy * cosyx - cosxx * cosyy)**2
        Pxy <- (-(cosxy * cosyy * Pxxr) + cosxx * cosyy * Pxyr + cosxy * cosyx * Pyxr - cosxx * cosyx * Pyyr) / (cosxy * cosyx - cosxx * cosyy)**2
#        Pyx <- (-(cosxy * cosyy * Pxxr) + cosxy * cosyx * Pxyr + cosxx * cosyy * Pyxr - cosxx * cosyx * Pyyr) / (cosxy * cosyx - cosxx * cosyy)**2

        # then solve the beta angle, first solve the Stokes parameters
#        I <- abs( Pxx ) + abs( Pyy )
        Q <- abs( Pxx ) - abs( Pyy )
        U <- 2 * Re( Pxy )
        V <- -2 * Im( Pxy )

        Ip <- sqrt( Q**2 + U**2 + V**2 )
        L <- Q + 1i*U

        A2 <- sqrt( .5 * ( Ip + abs(L) ) )
        B2 <- sqrt( .5 * ( Ip - abs(L) ) )

        beta2 <- acos( B2 / A2 )

        return( beta2 )

    }
