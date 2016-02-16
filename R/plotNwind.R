plotNwind <- function(nWind,timeRes=NULL,xlim=NULL,ylim=NULL,zlimE=c(-1,1)*500,zlimN=c(-1,1)*500,zlimU=c(-1,1)*100,plotEfield=TRUE,ylimEfield=NULL,pdf=NULL,figNum=NULL,width=8.27,height=2.9225,paper='special',tickRes=NULL,bg='white',fg='black',cex=1.0,col.regions=guisdap.colors,...){
    #
    # Plot the electric field and neutral wind returned by NeutralWinds
    # 
    # INPUT:
    #  nWind   output list from NeutralWinds
    #  timeRes time resolution to which we will integrate the neutral wind [s]
    #  xlim    x axis limits, UT hours counted from 00 UT on the day of the first data point
    #  ylim    y axis limits [UT]
    #  zlimE   z axis limits for the East velocity component
    #  zlimN   z axis limits for the North component
    #  zlimU   z axis limits for the upward velocity component
    #  plotEfield logical, should the electric field be plotted?
    #  ylimEfield y axis limits for the electric field plot [mV/m]
    #
    #  pdf     name of a pdf output file
    #  figNum  device number to use for plotting
    #  width   width of the figure panels
    #  height  height of a single figure panel
    #  paper   paper setting for pdf output
    #  tickRes time tick resolution in minutes
    #  bg      figure background color
    #  fg      figure foreground color
    #  cex     text scaling factor
    #  col.regions col.regions to use in the color plots
    #  
    # I. Virtanen 2016
    #

    # optional integration in time
    if(!is.null(timeRes)){
        nWindF <- filterNwind( nWind , timeRes , ... )
    }else{
        nWindF <- nWind
    }


    # time limits
    if(is.null(xlim)){
        tLim <- range(nWindF[["time"]])
    }else{
        tLim <- (floor(nWindF[["time"]][1] / 3600 / 24) * 24 + xlim) * 3600
    }
    
    # data points within tLim
    tInds <- which(nWindF[["time"]]>=tLim[1] & nWindF[["time"]]<=tLim[2])
    
    if( length(tInds)==0){
        warning("No data from the given time period")
        return()
    }
    
    # height limits
    if(is.null(ylim)){
        hLim <- range( nWindF[["height"]] , na.rm=TRUE )
    }else{
        hLim <- ylim
    }
    

    # number of figure panels (3 velocity components, plus the optional Efield)
    nFig <- ifelse(plotEfield,5,3)

    # height of the final plot window
    wHeight <- nFig*height

    # open the  proper device
    figList <- c(is.null(figNum),is.null(pdf))
    if(sum(figList) < 1 ) stop('Only one output device can be selected at a time')
    # a new x11 by defaul
    if(sum(figList) == 2) x11(width=width,height=wHeight)
    # new plot to an existing x11 window
    if(!is.null(figNum)) {dev.set(figNum);plot.new()}
    # a new pdf file
    if(!is.null(pdf)) pdf(file=paste(pdf,'.pdf',sep=''),paper=paper,width=width,height=wHeight)

    # plot layout
    layout(matrix(seq(2*nFig),ncol=2,byrow=T),widths=rep(c(.9/sqrt(cex),.1*sqrt(cex)),2*nFig))
    par(mar=c(3,3,2,1)*cex,mgp=c(1.5,.5,0)*cex)

    # tick marks in the time axis
    ticks <- timeTicks(tLim,tickRes)

    # the optional E field plots
    if(plotEfield){
        EmV <- nWind[["E"]]*1000
        stdmV <- t(sqrt(apply(nWind[["Ecov"]],FUN=diag,MARGIN=1)))*1000
        errLims1 <- EmV-stdmV
        errLims2 <- EmV + stdmV

        # adjustment to time axis for plotting
        ttE <- nWind[["time"]] - median(diff(nWind[["time"]]))/2
        
        plot(ttE,EmV[,1],xlim=tLim,xaxt='n',xlab='',ylab=expression(paste("E"[E]^{}," [mVm"[]^{-1},"]")),ylim=ylimEfield,type='n',cex.axis=cex,cex.lab=cex)
        abline(h=0,lwd=2)
        arrows(ttE,errLims1[,1],ttE,errLims2[,1],code=3,length=0,col='red',lwd=2)
        lines(ttE,EmV[,1],lwd=2)
        axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)

        plot.new()
        
        plot(ttE,EmV[,2],xlim=tLim,xaxt='n',xlab='',ylab=expression(paste("E"[N]^{}," [mVm"[]^{-1},"]")),ylim=ylimEfield,type='n',cex.axis=cex,cex.lab=cex)
        abline(h=0,lwd=2)
        arrows(ttE,errLims1[,2],ttE,errLims2[,2],code=3,length=0,col='red',lwd=2)
        lines(ttE,EmV[,2],lwd=2)
        axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)

        plot.new()
    }


    tt <- nWindF[["time"]] - median(diff(nWindF[["time"]]))/2
    
    # East velocity component
    image(tt[tInds],nWindF[["height"]],nWindF[["nWind"]][tInds,,1],xlim=tLim,ylim=hLim,zlim=zlimE,col=col.regions(1000),xaxt='n',ylab='Height [km]',xlab='',cex=cex,cex.lab=cex,cex.axis=cex)
    axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)
    image(c(0,1),seq(zlimE[1],zlimE[2],length.out=1000),t(matrix(rep(seq(zlimE[1],zlimE[2],length.out=1000),2),ncol=2)),col=col.regions(1000),ylab=expression(paste("V"[n]," East [ms"[]^{-1},"]")),xaxt='n',xlab='',cex=cex,cex.lab=cex,cex.axis=cex)

    # North velocity component
    image(tt[tInds],nWindF[["height"]],nWindF[["nWind"]][tInds,,2],xlim=tLim,ylim=hLim,zlim=zlimN,col=col.regions(1000),xaxt='n',ylab='Height [km]',xlab='',cex=cex,cex.lab=cex,cex.axis=cex)
    axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)
    image(c(0,1),seq(zlimN[1],zlimN[2],length.out=1000),t(matrix(rep(seq(zlimE[1],zlimE[2],length.out=1000),2),ncol=2)),col=col.regions(1000),ylab=expression(paste("V"[n]," North [ms"[]^{-1},"]")),xaxt='n',xlab='',cex=cex,cex.lab=cex,cex.axis=cex)

    # Up velocity component
    image(tt[tInds],nWindF[["height"]],nWindF[["nWind"]][tInds,,3],xlim=tLim,ylim=hLim,zlim=zlimU,col=col.regions(1000),xaxt='n',ylab='Height [km]',xlab='UTC',cex=cex,cex.lab=cex,cex.axis=cex)
    axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)
    image(c(0,1),seq(zlimU[1],zlimU[2],length.out=1000),t(matrix(rep(seq(zlimE[1],zlimE[2],length.out=1000),2),ncol=2)),col=col.regions(1000),ylab=expression(paste("V"[n]," Up [ms"[]^{-1},"]")),xaxt='n',xlab='',cex=cex,cex.lab=cex,cex.axis=cex)

    if(plotEfield){
        mtext( paste("Electric field and neutral wind ",substr( as.character( as.POSIXlt(nWindF[["time"]][tInds[1]],origin='1970-01-01',tz='utc') ) , 1 , 10 )) , side = 3, line = -1.5, outer = TRUE , cex = cex )
    }else{
        mtext( paste("Neutral wind ",substr( as.character( as.POSIXlt(nWindF[["time"]][tInds[1]],origin='1970-01-01',tz='utc') ) , 1 , 10 )) , side = 3, line = -1.5, outer = TRUE , cex = cex )
    }


    # THE X AXES OF E FIELD AND NEUTRAL WIND PLOTS DO NOT MATCH, AND THE XLIM ARGUMENT OF IMAGE DOES NOT SEEM TO WORK AS EXPECTED...+


    
    
    # if we did not plot on an x11 device, we must close the device properly
    if((sum(figList)==1)&is.null(figNum)) dev.off()

    return(invisible(nWindF))
}
