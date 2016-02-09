plotEfieldF <- function( Elist,xlim=NULL,ylim=c(-1,1)*max(abs(Elist[["E"]]),na.rm=T)*1000, cex=1 , pdf=NULL , figNum=NULL , width=8.27, height=5.845,paper='special',tickRes=NULL,bg='white',fg='black',main='NA',...){
    #
    # Plot the F-region electric field components returned by ElectricFieldsF
    #
    # I. Virtanen 2016
    #

    EmV <- Elist[["E"]]*1000
    stdmV <- t(sqrt(apply(Elist[["Ecov"]],FUN=diag,MARGIN=1)))*1000

    # lower limits of the errorbars
    errLims1 <- EmV-stdmV
    # upper limits
    errLims2 <- EmV + stdmV

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

    # open the  proper device
    figList <- c(is.null(figNum),is.null(pdf))
    if(sum(figList) < 1 ) stop('Only one output device can be selected at a time')
    # a new x11 by defaul
    if(sum(figList) == 2) x11(width=width,height=height)
    # new plot to an existing x11 window
    if(!is.null(figNum)) {dev.set(figNum);plot.new()}
    # a new pdf file
    if(!is.null(pdf)) pdf(file=paste(pdf,'.pdf',sep=''),paper=paper,width=width,height=height)


    layout(matrix(seq(2),ncol=1,byrow=T))
    par(mar=c(3,3,2,1)*cex,mgp=c(1.5,.5,0)*cex)

    ticks <- timeTicks(tLim,tickRes)

    
    # plot the data as lines
    plot(time,EmV[,1],xlim=tLim,xaxt='n',xlab='',ylab=expression(paste("E"[E]^{}," [mVm]"[]^{-1})),ylim=ylim,type='n',cex.axis=cex,cex.lab=cex)
    abline(h=0,lwd=2)
    arrows(time,errLims1[,1],time,errLims2[,1],code=3,length=0,col='red',lwd=2)
    lines(time,EmV[,1],lwd=2)
    axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)
    
    plot(time,EmV[,2],xlim=tLim,xaxt='n',xlab='UTC',ylab=expression(paste("E"[N]^{}," [mVm]"[]^{-1})),ylim=ylim,type='n',cex.axis=cex,cex.lab=cex)
    abline(h=0,lwd=2)
    arrows(time,errLims1[,2],time,errLims2[,2],code=3,length=0,col='red',lwd=2)
    lines(time,EmV[,2],lwd=2)
    axis(1,at=ticks$tick,labels=ticks$string,cex=cex,cex.lab=cex,cex.axis=cex)
    
    mtext( paste("F region electric field ",substr( as.character( as.POSIXlt(Elist[["time"]][tInds[1]],origin='1970-01-01',tz='utc') ) , 1 , 10 )) , side = 3, line = -1.5, outer = TRUE , cex = cex )


    if((sum(figList)==1)&is.null(figNum)) dev.off()

} # plotEfieldF



