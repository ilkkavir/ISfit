plotEvectorsPolar <- function( Elist , tlims=NULL, maxE=100 , cex=1 , width=8.27, height=8.27,paper='special',bg='transparent',fg='black',title=NA,dirname=NULL,...){
    #
    # Plot the F-region electric fields as arrows
    #
    # I. Virtanen 2016
    #

    EmV <- Elist[["E"]]*1000


    # times
    time <- as.numeric(Elist[["time"]])
    if(is.null(tlims)){
        tLim <- range(time)
    }else{
        tLim <- (floor(as.numeric(Elist[["time"]][1])/3600/24)*24+tlims)*3600
    }

    tInds <- which( (time>=tLim[1]) & (time<=tLim[2]) )

    if( length(tInds)==0){
        warning("No data from the given time period")
        return(invisible(data2))
    }

    if (is.null(dirname)){
        outpath <- paste('EvectorsPolar_',format(Elist$POSIXtime[[tInds[1]]],format="%Y%m%dT%H%M%S"),sep='')
    }else{
        outpath <- dirname
    }

    dir.create(outpath,recursive=TRUE,showWarnings=FALSE)

    Eabs <- sqrt(EmV[,1]**2+EmV[,2]**2)
    Eangle <- atan2(EmV[,1],EmV[,2])*180/pi

    rlims <- seq(0,maxE,by=25)

    for (tt in seq(1,length(tInds))){
        # open the figure file
        fname <- file.path(outpath,paste( 'frame' , paste(rep('0',(6-ceiling(log10(tt+.1)))) , sep='',collapse='') , tt , '.png' , sep='',collapse=''))

        png(fname)

        par(mar=c(3,1,2,1)*cex,mgp=c(1.5,.5,0)*cex)

        polar.plot(Eabs[tt],Eangle[tt],radial.lim=rlims)


        if(is.na(title)){
            mtext( paste(substr( as.character( as.POSIXlt(Elist[["time"]][tInds[tt]],origin='1970-01-01',tz='utc') ) , 1 , 20 )) , side = 3, line = -1.5, outer = TRUE , cex = cex )
        }else{
            mtext( title , side = 3, line = -1.5, outer = TRUE , cex = cex )
        }

        dev.off()
        cat('\r',tt)
    }


} # plotEfieldF



