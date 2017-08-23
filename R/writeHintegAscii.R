writeHintegAscii <- function( PPlist , par='ViRproj' , hmin=200 , hmax=400 , fname=paste(par,'_height-integrated_',hmin,'-',hmax,'km','.dat',sep='')){
    #
    # Write a time series of a height-integrated plasma parameter to file
    #
    # INPUT:
    #  PPlist  an oridinal list of plasma parameters from readPP.3D (or from projectParams)
    #  par     name of the parameter to integrate in height
    #  hmin    lower boundary for integration [km]
    #  hmax    upper boundary for integration [km]
    #  fname   output file name
    #
    #
    # IV 2017
    #

    # the post-integration
    cat('Integrating data in height...\n')

    # first remove any data points from below hmin and above hmax
    for(dd in seq(dim(PPlist[["height"]])[2])){
        PPlist[["param"]][PPlist[["height"]][,dd]<hmin,,dd] <- NA
        PPlist[["param"]][PPlist[["height"]][,dd]>hmax,,dd] <- NA
    }
    # then actually integrate
    ppinteg <- postInteg.height(PPlist,par)

    # then write the integrated data to file

    fcon <- file( fname , 'w' )

    dims <- dim(PPlist[["param"]])
    nt <- dims[3]

    cat(sprintf("%-1s%13s %10s %10s\n","#","unixtime",par,"std"),file=fcon,append=FALSE)

    for( t in seq(nt) ){
        cat(sprintf("%14.3f %10.0f %10.0f\n",as.numeric(PPlist[["POSIXtime"]][[t]]),ppinteg[["p"]][t],ppinteg[["std"]][t]),file=fcon,append=TRUE)
    }

    close(fcon)

    return(invisible(ppinteg))

}
