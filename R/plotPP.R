plotPP <- function(fpath,par=list(Ne=c(10,12),Ti=c(0,3000),Te=c(0,4000),ViR1=c(-200,200),azel=T),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,ytype='height'){
# 
# read plasma parameters from file(s) and plot them
# 
# I. Virtanen 2010, 2013
#

  warning("this function has not yet been fully modified to work with the 3D parameter fit results")
  
  # read all data
  datas <- readPP(fpath)

  # take a copy that will be returned
  datas2 <- datas

  # if no data, return NULL
  if(is.null(datas)) invisible(NULL)
  
  # cut off large variances (only Ne is studied)
  if(!is.null(stdThreshold)){
    thrInds <- (datas$std[,1,] / datas$param[,1,]) >  stdThreshold
    for(k in seq(dim(datas[["param"]])[2])) datas$param[,k,][thrInds] <- NA
  }

  # cut off all data where Ne is below NeMin
  if(!is.null(NeMin)){
    NeInds <- datas$param[,1,] < NeMin
    for(k in seq(dim(datas$param)[2])) datas$param[,k,][NeInds] <- NA
  }

  # cut off all failed iterations
  statInds <- datas$status != 0
  for(k in seq(dim(datas$param)[2])) datas$param[,k,][statInds] <- NA

  # cut off very large chi-squared values
  if(!is.null(chisqrLim)){
    chiInds <- datas$chisqr >  chisqrLim
    for(k in seq(dim(datas$param)[2])) datas$param[,k,][chiInds] <- NA
  }


  # plot the initial(model) parameters instead of the fit results
  if(model){
    datas$param   <- datas$model
    datas$std[,,] <- 0
  }

  # number of figures panels
  nFig <- length(par)

  # time limits
  if(is.null(xlim)){
    tLim <- range(datas[["time_sec"]])
  }else{
    tLim <- (floor(datas$time_sec[1] / 3600 / 24) * 24 + xlim) * 3600
  }

  # height limits
  if(is.null(ylim)){
    hLim <- range(datas[["height"]])
  }else{
    hLim <- ylim
  }
  
  # height of the full plot window
  wHeight <- min(nFig,4)*height

  # open the  proper device
  figList <- c(is.null(figNum),is.null(pdf),is.null(jpg))
  if(sum(figList) < 2 ) stop('Only one output device can be selected at a time')
  # a new x11 by defaul
  if(sum(figList) == 3) x11(width=width,height=wHeight)
  # new plot to an existing x11 window
  if(!is.null(figNum)) {dev.set(figNum);plot.new()}
  # a new pdf file
  if(!is.null(pdf)) pdf(file=paste(pdf,'.pdf',sep=''),paper=paper,width=width,height=wHeight)
  # a new jpeg file
  if(!is.null(jpg)){
    jpeg(file=paste(jpg,'.jpg',sep=''),width=width,height=wHeight,units='in',res=res)
    cex = cex*res/72
  }
  
  treold <- trellis.par.get()
  trenew <- treold
  trenew$bacground$col <- bg
  trenew$layout.heights$strip <- 0
  trenew$layout.heights$xlab <- 0
  trenew$layout.heights$sub <- 0
  trenew$layout.heights$between <- 0
  trenew$layout.heights$main <- .1
  trenew$layout.heights$top.padding <- 0
  trenew$layout.heights$main.key.padding <- 0
  trenew$layout.heights$key.top <- 1
  trenew$layout.heights$bottom.padding <- 0
  trenew$layout.heights$xlab.key.padding <- 0
  trenew$layout.widths$panel <- 5
  trenew$layout.widths$right.padding <- 1
  trenew$layout.heights$key.axis.padding<-1
  trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights,layout.widths=trenew$layout.widths))
  
  # tick marks in the time axis
  ticks <- timeTicks(tLim,tickRes)

  # actual plotting
  curFig <- 1

  for(p in seq(length(par))){

    # select whether x-axis is logarithmic or not
    xlog <- switch( names(par[p]) , Ne=TRUE , Coll=TRUE , FALSE )

    # plot titles, there are spare ones for future systems.. :)
    main <- switch( names(par[p]) ,
                   Ne=list(expression(paste("Log","(N"[e]^{},") [m"[]^{-3},"]")),col=fg,cex=cex,vjust=1.2),
                   Ti=list(expression(paste("T"[i]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                   Tipar=list(expression(paste("T"["i||"]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                   Tiperp=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                   Te=list(expression(paste("T"[e]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                   Tepar=list(expression(paste("T"["e||"]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                   Teperp=list(expression(paste("T"[paste("e",symbol("\136"))]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                   Coll=list(expression(paste("Log(",nu['in']^{},") [s"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   Vix=list(expression(paste("V"[iX]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   Viy=list(expression(paste("V"[iY]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   Viz=list(expression(paste("V"[iZ]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViB=list(expression(paste("V"["i||"]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR1=list(expression(paste("V"[iR1]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR2=list(expression(paste("V"[iR2]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR3=list(expression(paste("V"[iR3]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR4=list(expression(paste("V"[iR4]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR5=list(expression(paste("V"[iR5]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR6=list(expression(paste("V"[iR6]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR7=list(expression(paste("V"[iR7]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR8=list(expression(paste("V"[iR8]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR9=list(expression(paste("V"[iR9]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   ViR10=list(expression(paste("V"[iR10]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                   Ion1=list(sprintf("Fraction of ion mass %.1f u",datas$mi[1]),col=fg,cex=cex,vjust=1.2),
                   Ion2=list(sprintf("Fraction of ion mass %.1f u",datas$mi[2]),col=fg,cex=cex,vjust=1.2),
                   Ion3=list(sprintf("Fraction of ion mass %.1f u",datas$mi[3]),col=fg,cex=cex,vjust=1.2),
                   Ion4=list(sprintf("Fraction of ion mass %.1f u",datas$mi[4]),col=fg,cex=cex,vjust=1.2),
                   Ion5=list(sprintf("Fraction of ion mass %.1f u",datas$mi[5]),col=fg,cex=cex,vjust=1.2)
                   )
    # azimuth and elevation
    if(names(par[p])=="azel"){
      if(par[[p]]){
        trenew$layout.widths$right.padding <- 10.5
        trenew$layout.heights$key.axis.padding<--1
        trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights,layout.widths=trenew$layout.widths))
        curFig <- addAzElPlot(
          az     = datas$azelT[1,]%%360,
          el     = datas$azelT[2,],
          t      = datas$time_sec,
          xlim   = tLim,
          cex    = cex,
          ticks  = ticks,
          nFig   = nFig,
          curFig = curFig
          )
        trenew$layout.widths$right.padding <- 1
        trenew$layout.heights$key.axis.padding<-1
        trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights,layout.widths=trenew$layout.widths))
      }
    }else{
      d <- datas[["param"]][,names(par[p]),]
      if( xlog ) d <- log10( d )
      
      curFig <- addPPcolorPlot(
        d = d,
        h = datas[["height"]],
        t = datas[["time_sec"]],
        xlim = tLim,
        ylim = hLim,
        zlim = par[[p]],
        main = main,
        cex = cex,
        ticks = ticks,
        nFig = nFig,
        curFig = curFig,
        col.regions = col.regions
        )
    }
    
  }
  
  # if we did not plot on an x11 device, we must close the device properly
  if((sum(figList)==2)&is.null(figNum)) dev.off()

  invisible(datas2)

} # plotPP


