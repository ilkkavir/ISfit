plotPlasmaParams <- function(fpath,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',Ne=T,Ti=T,Te=T,Vi=T,Ve=F,coll=F,pIon=c(F,F,F),NeLim=c(9,12),TiLim=c(0,2000),TeLim=c(0,2000),ViLim=c(-200,200),VeLim=c(-200,200),collLim=c(0,5),pIonLim=c(0,1),rLim=NULL,hLim=NULL,tLim=NULL,stdThreshold=.5,tickRes=NULL,plotModel=F,padModel=F,chisqrLim=10,NeMin=9,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors){
# 
# read plasma parameters from file(s) and plot them
# 
# I. Virtanen 2010
# 
  # read all data
  datas <- readPPIdir(fpath)

  # if no data, return NULL
  if(is.null(datas)) invisible(NULL)
  

  # find the parameters the user wants to plot
  rInds <- rep(T,datas$nRange)
  pInds <- c(Ne,Ti,Te,coll,Vi,Ve,pIon)
  tInds <- rep(T,datas$n)

  hLim2 <- range(datas$height,na.rm=T)
  if(!is.null(hLim)){ rInds <- (datas$height >= hLim[1]) & (datas$height <= hLim[2]); hLim2 <- hLim}
  if(!is.null(rLim)){ rInds <- rInds & ((datas$range >= rLim[1]) & (datas$range <= rLim[2])); hLim2 <- range(datas$height[rInds])}

  tLim2 <- range(datas$time_sec,na.rm=T)
  if(!is.null(tLim)){
    tLim2 <- (floor(datas$time_sec[1] / 3600 / 24) * 24 + tLim) * 3600
    tInds <- (datas$time_sec >= tLim2[1]) & (datas$time_sec <= tLim2[2])
  }

  if(!any(rInds)) stop('No data from the given range interval')
  if(!any(pInds)) stop('No parameters to fit')
  if(!any(tInds)) stop('No data from the given time interval')

  # cut off large variances (only Ne is studied)
  if(!is.null(stdThreshold)){
    thrInds <- (datas$std[,1,] / datas$param[,1,]) >  stdThreshold
    for(k in seq(dim(datas[["param"]])[2])) datas$param[,k,][thrInds] <- NA
  }

  # cut off all data where Ne is below NeMin
  if(!is.null(NeMin)){
    NeInds <- log10(datas$param[,1,]) < NeMin
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
  if(plotModel){
    datas$param   <- datas$model
    datas$std[,,] <- 0
  }

  # if only the missing values are taken fromt the model
  if(padModel){
    datas$param[is.na(datas$std)] <- datas$model[is.na(datas$std)]
    datas$std[is.na(datas$std)]   <- 0
  }

  # how many frames do we actuallly have?
  nFig <- sum(pInds)


  # if there were several parameter files, make a color plot
  if(datas$n > 1){

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
    trenew$layout.heights$strip=0
    trenew$layout.heights$xlab=0
    trenew$layout.heights$sub = 0
    trenew$layout.heights$between=0
    trenew$layout.heights$main=.1
    trenew$layout.heights$top.padding=0
    trenew$layout.heights$main.key.padding=0
    trenew$layout.heights$bottom.padding=0
    trenew$layout.heights$xlab.key.padding=0
    trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights))

    # tick marks in the time axis
    ticks <- timeTicks(tLim2,tickRes)

    # actual plotting
    curFig <- 1

    # electron density
    if(Ne) curFig <- addPPcolorPlot(
                                      d      = log10(datas$param[rInds,1,tInds]),  # the data
                                      h      = datas$height[rInds],                # heights
                                      t      = datas$time_sec[tInds],              # times
                                      xlim   = tLim2,                              # x-axis (time) limits
                                      ylim   = hLim2,                              # y-axis (height) limits
                                      zlim   = NeLim,                              # z-axis (colorcoded Ne) limits
                                      # title of the plot
                                      main   = list(expression(paste("Log","(N"[e]^{},") [m"[]^{-3},"]")),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,                                # maginification of axis lines and axis labels
                                      ticks  = ticks,                              # tick mark locations and label texts
                                      nFig   = nFig,                               # total number of figure panels
                                      curFig = curFig,                             # number of the current panel
                                 col.regions = col.regions
                                    )

    # electron temperature (datas$param[,2,] is the ratio Te/Ti)
    if(Te) curFig <- addPPcolorPlot(
                                      d      = datas$param[rInds,2,tInds]*datas$param[rInds,3,tInds],
                                      h      = datas$height[rInds],
                                      xlim   = tLim2,
                                      ylim   = hLim2,
                                      t      = datas$time_sec[tInds],
                                      zlim   = TeLim,
                                      main   = list(expression(paste("T"[e]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,
                                      ticks  = ticks,
                                      nFig   = nFig,
                                      curFig = curFig,
                                 col.regions = col.regions
                                    )

    # ion temperature
    if(Ti) curFig <- addPPcolorPlot(
                                      d      = datas$param[rInds,2,tInds],
                                      h      = datas$height[rInds],
                                      t      = datas$time_sec[tInds],
                                      zlim   = TiLim,
                                      xlim   =  tLim2,
                                      ylim   = hLim2,
                                      main   = list(expression(paste("T"[i]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,
                                      ticks  = ticks,
                                      nFig   = nFig,
                                      curFig = curFig,
                                 col.regions = col.regions
                                    )

    # ion-neutral collision frequency
    if(coll) curFig <- addPPcolorPlot(
                                      d      = log10(datas$param[rInds,4,tInds]),
                                      h      = datas$height[rInds],
                                      xlim   = tLim2,
                                      ylim   = hLim2,
                                      t      = datas$time_sec[tInds],
                                      zlim   = collLim,
                                      main   = list(expression(paste("Log(",nu['in']^{},") [s"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,
                                      ticks  = ticks,
                                      nFig   = nFig,
                                      curFig = curFig,
                                 col.regions = col.regions
                                    )

    # line-of-sight ion velocity
    if(Vi) curFig <- addPPcolorPlot(
                                      d      = datas$param[rInds,5,tInds],
                                      h      = datas$height[rInds],
                                      t      = datas$time_sec[tInds],
                                      zlim   = ViLim,
                                      xlim   = tLim2,
                                      ylim   = hLim2,
                                      main   = list(expression(paste("V"[i]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,
                                      ticks  = ticks,
                                      nFig   = nFig,
                                      curFig = curFig,
                                 col.regions = col.regions
                                   )

    # line-of-sight electron velocity (datas$param[,6,] is the difference Ve-Vi)
    if(Ve) curFig <- addPPcolorPlot(
                                      d      = (datas$param[rInds,6,tInds]+datas$param[rInds,5,tInds]),
                                      h      = datas$height[rInds],
                                      t      = datas$time_sec[tInds],
                                      zlim   = VeLim,
                                      xlim   = tLim2,
                                      ylim   = hLim2,
                                      main   = list(expression(paste("V"[e]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,
                                      ticks  = ticks,
                                      nFig   = nFig,
                                      curFig = curFig,
                                 col.regions = col.regions
                                   )

    # ion abundances
    for(k in seq(length(pIon))){

      # pIon[1] is the abundance of ion mass 16u and pIon[2] that of ion mass 1u (H+)
      # the remaining fraction is given for the ion mass 30.5u (NO+ and O2+)
      if(pIon[k]){
        curFig <- addPPcolorPlot(
                                      d      = datas$param[rInds,(6+k),tInds],h=datas$height[rInds],
                                      t      = datas$time_sec[tInds],
                                      zlim   = pIonLim,
                                      xlim   = tLim2,
                                      ylim   = hLim2,
                                      main   = list(sprintf("Fraction of ion mass %.1f u",datas$mi[k+1]),col=fg,cex=cex,vjust=1.2),
                                      cex    = cex,
                                      ticks  = ticks,
                                      nFig   = nFig,
                                      curFig = curFig,
                                 col.regions = col.regions
                                )
      }
    }

  # if only one result file was found, then plot the individual points
  }else{

    
    # try to make equal number of rows and columns, more columns than rows if necessary
    nRows   <- min(2,ceiling(nFig/2))
    nCols   <- ceiling(nFig/nRows)


    # height of the full plot window
    wHeight <- nRows*height*2
    wWidth  <- nCols*width/2

    # open the  proper device
    figList <- c(is.null(figNum),is.null(pdf),is.null(jpg))
    if(sum(figList) < 2 ) stop('Only one output device can be selected at a time')

    # a new x11 by defaul
    if(sum(figList) == 3) x11(width=wWidth,height=wHeight)

    # new plot to an existing x11 window
    if(!is.null(figNum)) {dev.set(figNum);plot.new()}

    # a new pdf file
    if(!is.null(pdf)) pdf(file=paste(pdf,'.pdf',sep=''),paper=paper,width=wWidth,height=wHeight)

    # a new jpeg file
    if(!is.null(jpg)) jpeg(file=paste(jpg,'.jpg',sep=''),width=wWidth,height=wHeight,units='in',res=res)

    # set margins and maginification, plus divide the plot area
    par(cex.lab=cex,cex.axis=cex,bg=bg,mar=c(5,5,5,1),fg=fg,bg=bg,col.lab=fg,col.axis=fg)
    layout(matrix(seq(nRows*nCols),nrow=nRows))

    # actual plotting

    if(Ne) addPPlinePlot(
                          d    = datas$param[rInds,1,tInds],                                               # the data
                          err  = 3*datas$std[rInds,1,tInds],                                               # 3-sigma errorbars
                          log  = T,                                                                        # logarithmic scale in x-axis
                          h    = datas$height[rInds],                                                      # heigts
                          xlim = NeLim,                                                                    # limits of the x-axis
                          ylim = hLim2,                                                                    # limits of the y-axis
                          xlab = list(expression(paste("Log","(N"[e]^{},") [m"[]^{-3},"]")),col=fg,cex=cex)# x-axis label
                        )


    if(Te) addPPlinePlot(
                          d    = datas$param[rInds,2,tInds]*datas$param[rInds,3,tInds],
                          err  = 3*(datas$param[rInds,2,tInds]*datas$std[rInds,3,tInds] + 
                                    datas$param[rInds,3,tInds]*datas$std[rInds,2,tInds]),
                          log  = F,
                          h    = datas$height[rInds],
                          xlim = TeLim,
                          ylim = hLim2,
                          xlab = list(expression(paste("T"[e]^{}," [K]")),col=fg,cex=cex)
                        )

    if(Ti) addPPlinePlot(
                          d    = datas$param[rInds,2,tInds],
                          err  = 3*datas$std[rInds,2,tInds],
                          log  = F,
                          h    = datas$height[rInds],
                          xlim = TiLim,
                          ylim = hLim2,
                          xlab = list(expression(paste("T"[i]^{}," [K]")),col=fg,cex=cex)
                        )

    if(coll) addPPlinePlot(
                          d    = datas$param[rInds,4,tInds],
                          err  = 3*datas$std[rInds,4,tInds],
                          log  = T,
                          h    = datas$height[rInds],
                          xlim = collLim,
                          ylim = hLim2,
                          xlab = list(expression(paste("Log(",nu['in']^{},") [s"[]^{-1},"]")),col=fg,cex=cex)
                          )

    if(Vi) addPPlinePlot(
                          d    = datas$param[rInds,5,tInds],
                          err  = 3*datas$std[rInds,5,tInds],
                          log  = F,
                          h    = datas$height[rInds],
                          xlim = ViLim,
                          ylim = hLim2,
                          xlab = list(expression(paste("V"[i]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex)
                        )

    if(Ve) addPPlinePlot(
                          d    = (datas$param[rInds,6,tInds]+datas$param[rInds,5,tInds]),
                          err  = 3*(datas$std[rInds,6,tInds]+datas$std[rInds,5,tInds]),
                          log  = F,
                          h    = datas$height[rInds],
                          xlim = VeLim,
                          ylim = hLim2,
                          xlab = list(expression(paste("V"[e]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex)
                        )

    for(k in seq(length(pIon))){
      if(pIon[k]){
        addPPlinePlot(
                          d    = datas$param[rInds,(6+k),tInds],
                          err  = 3*datas$std[rInds,(5+k),tInds],
                          log  = F,
                          h    = datas$height[rInds],
                          xlim = pIonLim,
                          ylim = hLim2,
                          xlab =list(sprintf("Fractional abundance of ion mass %.1f u",datas$mi[k+1]),col=fg,cex=cex)
                        )
      }
    }

  }

  # if we did not plot on an x11 device, we must close the device properly
  if((sum(figList)==2)&is.null(figNum)) dev.off()

  invisible(datas)

} # plotPlasmaParams


