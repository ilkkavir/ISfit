plotEvectorsF <- function( Elist,xlim=NULL,tmV=1000, cex=1 , pdf=NULL , figNum=NULL , width=8.27, height=2.9225,paper='special',tickRes=NULL,bg='white',fg='black',title=NA,opendev=TRUE,closedev=TRUE,...){
    #
    # Plot the F-region electric fields as arrows
    #
    # I. Virtanen 2016
    #

    EmV <- Elist[["E"]]*1000


    # times
    time <- as.numeric(Elist[["time"]])
    if(is.null(xlim)){
        tLim <- range(time)
    }else{
        tLim <- (floor(as.numeric(Elist[["time"]][1])/3600/24)*24+xlim)*3600
    }

    tInds <- which( (time>=tLim[1]) & (time<=tLim[2]) )

    if( length(tInds)==0){
        warning("No data from the given time period")
        return(invisible(data2))
    }

    if(opendev){
        # open the  proper device
        figList <- c(is.null(figNum),is.null(pdf))
        if(sum(figList) < 1 ) stop('Only one output device can be selected at a time')
        # a new x11 by defaul
        if(sum(figList) == 2) x11(width=width,height=height)
        # new plot to an existing x11 window
        if(!is.null(figNum)) {dev.set(figNum);plot.new()}
        # a new pdf file
        if(!is.null(pdf)) pdf(file=paste(pdf,'.pdf',sep=''),paper=paper,width=width,height=height)
    }

    par(mar=c(3,1,2,1)*cex,mgp=c(1.5,.5,0)*cex)

    ticks <- timeTicks(tLim,tickRes)



    # The full x-axis scale is 1 V/m, scale the times accordingly
    tscale <- tmV/diff(tLim)

    # plot zeros, these will be the vector start points
    plot(time[tInds]*tscale,EmV[tInds,1]*0,xlim=(tLim*tscale+c(-.1,.25)*tmV),xaxt='n',xlab='UTC',ylab='',yaxt='n',asp=1,pch=20,cex=.5,bty='n')
    lines(tLim*tscale,c(0,0),lwd=1)
    arrows(time[tInds]*tscale,EmV[tInds,1]*0,time[tInds]*tscale+EmV[tInds,1],EmV[tInds,2],length=.05)
    
    arrows((tLim[2]*tscale+.18*tmV),.05*tmV,(tLim[2]*tscale+.23*tmV),.05*tmV,length=.05)
    arrows((tLim[2]*tscale+.18*tmV),.05*tmV,(tLim[2]*tscale+.13*tmV),.05*tmV,length=.05)
    arrows((tLim[2]*tscale+.18*tmV),.05*tmV,(tLim[2]*tscale+.18*tmV),.1*tmV,length=.05)
    arrows((tLim[2]*tscale+.18*tmV),.05*tmV,(tLim[2]*tscale+.18*tmV),0,length=.05)
    text((tLim[2]*tscale+.25*tmV),.05*tmV,'E')
    text((tLim[2]*tscale+.11*tmV),.05*tmV,'W')
    text((tLim[2]*tscale+.18*tmV),.12*tmV,'N')
    text((tLim[2]*tscale+.18*tmV),-.02*tmV,'S')

    arrows((tLim[2]*tscale+.13*tmV),-.08*tmV,(tLim[2]*tscale+.23*tmV),-.08*tmV,length=.05)
    text((tLim[2]*tscale+.18*tmV),-.06*tmV,paste(round(tmV/10),' mV/m'),cex=1)
    abline(v=ticks$tick*tscale,lty=3,col="#707070")

    axis(1,at=ticks$tick*tscale,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)
    

    if(is.na(title)){
        mtext( paste("F region electric field ",substr( as.character( as.POSIXlt(Elist[["time"]][tInds[1]],origin='1970-01-01',tz='utc') ) , 1 , 10 )) , side = 3, line = -1.5, outer = TRUE , cex = cex )
    }else{
        mtext( title , side = 3, line = -1.5, outer = TRUE , cex = cex )
    }

    if(closedev){
        if((!is.null(pdf))&is.null(figNum)) dev.off()
    }

} # plotEfieldF



