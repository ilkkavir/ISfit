projectParams <- function( par , llhT , azelT , llhR  ){
    #
    # Project vector parameters (Ion velocity and temperatures) to a(n) (imaginary) site
    # located at llhR, which receives signal scattered from transmission from a site located
    # at llhT, pointed azelT. The site is added to the input list as an additional
    # site "Rproj", which can be subsequently handled in the same way as the other receivers.
    #
    # INPUT:
    #   par    a parameter list from readPP.3D
    #   llhT   transmitter latitude (deg), longitude (deg) and ellipsoid height (m),
    #          with respect to the WGS84 ellipsoid
    #   azelT  azimuth (deg east) and elevation (deg) of the transmitter beam
    #   llhR   latitude (deg), longitude (deg) and ellipsoid height (m), with respect to the WGS84 ellipsoid
    # OUTPUT:
    #   parout modified parameters list, with the site Rproj added
    #
    #
    # Details:
    #   The system calculates azimuth and elevation for TX and RX beams based on llhT, llhR,
    #   and scattering volume locations in the data list. The TR pair does not need to exist in
    #   reality.
    #
    # IV 2017
    #
    parout <- par

    # dimensions of the parater array etc.
    dimparam <- dim(par[["param"]])
    dimsites <- dim(par[["sites"]])
    dimcovar <- dim(par[["covar"]])

    # dimension names
    namesparam <- dimnames(par[["param"]])
    namescovar <- dimnames(par[["covar"]])

    # create larger arrays
    parout[["param"]] <- array(NA,dim=dimparam+c(0,4,0))
    parout[["std"]] <- array(NA,dim=dimparam+c(0,4,0))
    parout[["model"]] <- array(NA,dim=dimparam+c(0,4,0))
    parout[["sites"]] <- array(NA,dim=dimsites+c(1,0,0))
    parout[["covar"]] <- array(NA,dim=dimcovar+c(0,4,4,0))

    # copy the data
    parout[["param"]][,1:dimparam[2],] <- par[["param"]]
    parout[["std"]][,1:dimparam[2],] <- par[["std"]]
    parout[["model"]][,1:dimparam[2],] <- par[["model"]]
    parout[["sites"]][1:dimsites[1],,] <- par[["sites"]]
    parout[["covar"]][,1:dimcovar[2],1:dimcovar[3],] <- par[["covar"]]

    # udpdate dimnmames
    projnames <- c('ViRproj','ViRprojhor','TiRproj','TeRproj')
    dimnames(parout[["param"]]) <- list(namesparam[[1]],c(namesparam[[2]],projnames),namesparam[[3]])
    dimnames(parout[["std"]]) <- list(namesparam[[1]],c(namesparam[[2]],projnames),namesparam[[3]])
    dimnames(parout[["model"]]) <- list(namesparam[[1]],c(namesparam[[2]],projnames),namesparam[[3]])
    dimnames(parout[["covar"]]) <- list(namescovar[[1]],c(namescovar[[2]],projnames),c(namescovar[[3]],projnames),namescovar[[4]])


    # the actual projections. For each height-gate, we
    # 1. Calculate azelR based on llhT, azelT, and height from the parameter list
    # 2. Use llhT, azelT, llhR, and azelR to calculate the k-vector
    # 3. Project Ti, Te, Vi, and their covariances to k
    # These are done at every single height gate and integration
    for(dd in seq(dim(par[["param"]])[3])){

        # add the new site to parout$sites, we do not have frequency and azel
        parout[["sites"]][dimsites[1]+1,,dd] <- c(dimsites[1]+1,NA,llhT,azelT,llhR,NA,NA)

        for(hh in seq(dim(par[["height"]])[1])){

            # monostatic range from transmitter to target
            rangeT <- height2range( llhT=llhT , azelT=azelT , h=par[["height"]][hh,dd]*1000)

            # target coordinates
            llhS <- range2llh(r=rangeT, llhT=llhT, azelT=azelT)

            # azimuth, elevation of the receiver beam
            azelR <- llhTarget2azelrBeam( llhTarget=llhS , llhSite=llhR )

            #  beam intersection with guessed widths etc.
            isect <- beamIntersection(llhT=llhT, llhR=llhR, azelT=azelT, azelR=azelR,
                                      fwhmT=1, fwhmR=1, phArrT=TRUE, phArrR=TRUE,
                                      freq.Hz=224e6, stdXY = FALSE, rMonost = 1e+05)
            # pick the k-vector
            kvec <- isect[["k.ENU"]]

            # project velocity to the new site
            parout[["param"]][hh,'ViRproj',dd] <- par[["param"]][hh,7:9,dd]%*%kvec/sqrt(sum(kvec**2))
            parout[["std"]][hh,'ViRproj',dd] <- sqrt(kvec%*%par[["covar"]][hh,7:9,7:9,dd]%*%kvec/sum(kvec**2))
            # horizontal velocity component
            khor <- c(kvec[c(1,2)],0)
            khor <- khor / sqrt(sum(khor**2))
            parout[["param"]][hh,'ViRprojhor',dd] <- par[["param"]][hh,7:9,dd]%*%khor
            parout[["std"]][hh,'ViRprojhor',dd] <- sqrt(khor%*%par[["covar"]][hh,7:9,7:9,dd]%*%khor)

            # ion temperature
            phi <- radarPointings:::vectorAngle.cartesian(kvec,par[["B"]][hh,,dd],degrees=FALSE)
            prvec <- c( cos(phi)**2 , sin(phi)**2 )
            parout[["param"]][hh,'TiRproj',dd] <- par[["param"]][hh,2:3,dd]%*%prvec
            parout[["std"]][hh,'TiRproj',dd] <- sqrt(prvec%*%par[["covar"]][hh,2:3,2:3,dd]%*%prvec)

            # electron temperature
            parout[["param"]][hh,'TeRproj',dd] <- par[["param"]][hh,4:5,dd]%*%prvec
            parout[["std"]][hh,'TeRproj',dd] <- sqrt(prvec%*%par[["covar"]][hh,4:5,4:5,dd]%*%prvec)

        }
    }

    return(parout)
}
