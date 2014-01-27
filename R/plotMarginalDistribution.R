plotMarginalDistribution <- function(fpath,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',Ne=T,Ti=T,Te=T,Vi=T,Ve=F,coll=F,pIon=c(F,F,F),NeLim=c(9,12),TiLim=c(0,2000),TeLim=c(0,2000),ViLim=c(-200,200),VeLim=c(-200,200),collLim=c(0,5),pIonLim=c(0,1),rLim=NULL,hLim=NULL,tLim=NULL,bg='white',fg='black',res=300,cex=1.0,allPoints=F,confLimits=c(.682,.954,.996),burnin=0.1,cutModel=T){
# 
# Plot the marginal posterior distributions from a MCMC analysis. 
# Either individual points (X11 only to avoid large files) or color-coded probability regions
# The possibility for large files with individual points plotted is easy to add, if someone really needs that
# 
# I. Virtanen 2010
# 

  # load the data
  load(fpath)


  # find the parameters the user wants to plot
  rInds <- rep(T,length(PP$range))
  pInds <- c(Ne,Ti,Te,coll,Vi,Ve,pIon)

  hLim2 <- range(PP$height,na.rm=T)
  if(!is.null(hLim)){ rInds <- (PP$height >= hLim[1]) & (PP$height <= hLim[2]); hLim2 <- hLim}
  if(!is.null(rLim)){ rInds <- rInds & ((PP$range >= rLim[1]) & (PP$range <= rLim[2])); hLim2 <- range(PP$height[rInds])}

  tLim2 <- range(PP$time_sec,na.rm=T)
  if(!is.null(tLim)){
    tLim2 <- (floor(PP$time_sec[1] / 3600 / 24) * 24 + tLim) * 3600
    tInds <- (PP$time_sec >= tLim2[1]) & (PP$time_sec <= tLim2[2])
  }

  if(!any(rInds)) stop('No data from the given range interval')
  if(!any(pInds)) stop('No parameters to fit')

  # cut off the model values
  if(cutModel){
    for(r in seq(length(PP$range))){
      for(p in seq(dim(PP$param)[2])){
        if(is.na(PP$param[r,p])) PP$mcmc[r,p,] <- NA
      }
    }
  }

  # how many frames do we actuallly have?
  nFig <- sum(pInds)


  # try to make equal number of rows and columns, more columns than rows if necessary
  nRows   <- min(2,ceiling(nFig/2))
  nCols   <- ceiling(nFig/nRows)

  # height of the full plot window
  wHeight <- nRows*height*2
  wWidth  <- nCols*width/2

  # open the  proper device
  figList <- c(is.null(figNum),is.null(pdf),is.null(jpg))
  # all points can be plotted only to x11 window
  if(allPoints & ((!is.null(pdf))|(!is.null(jpg)))) cat('The points can be plotted only to x11 device, set allPoints=F to plot to a file.')
  if(allPoints) figList <- rep(T,3)

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
  par(cex.lab=cex,cex.axis=cex,bg=bg,mar=c(5,5,5,1))
  layout(matrix(seq(nRows*nCols),nrow=nRows))

  # number of burnin points
  Nburnin <- burnin*PP$Nmcmc
  # actual plotting
  if(Ne) addDistPlot(
                        d    = PP$mcmc[rInds,1,(Nburnin+1):PP$Nmcmc],                                   # the data
                        log  = T,                                                                          # logarithmic scale in x-axis
                        h    = PP$height[rInds],                                                           # heigts
                        xlim = NeLim,                                                                      # limits of the x-axis
                        ylim = hLim2,                                                                      # limits of the y-axis
                        xlab = list(expression(paste("Log","(N"[e]^{},") [m"[]^{-3},"]")),col=fg,cex=cex), # x-axis label
                      points = allPoints,
                  confLimits = confLimits
                    )

  if(Ti) addDistPlot(
                        d    = PP$mcmc[rInds,2,(Nburnin+1):PP$Nmcmc],
                        log  = F,
                        h    = PP$height[rInds],
                        xlim = TiLim,
                        ylim = hLim2,
                        xlab = list(expression(paste("T"[i]^{}," [K]")),col=fg,cex=cex),
                      points = allPoints,
                  confLimits = confLimits      
                      )

  if(Te) addDistPlot(
                        d    = PP$mcmc[rInds,2,(Nburnin+1):PP$Nmcmc]*PP$mcmc[rInds,3,(Nburnin+1):PP$Nmcmc],
                        log  = F,
                        h    = PP$height[rInds],
                        xlim = TeLim,
                        ylim = hLim2,
                        xlab = list(expression(paste("T"[e]^{}," [K]")),col=fg,cex=cex),
                      points = allPoints,
                  confLimits = confLimits      
                      )

    if(coll) addDistPlot(
                        d    = PP$mcmc[rInds,4,(Nburnin+1):PP$Nmcmc],
                        log  = T,
                        h    = PP$height[rInds],
                        xlim = collLim,
                        ylim = hLim2,
                        xlab = list(expression(paste("Log(",nu['in']^{},") [s"[]^{-1},"]")),col=fg,cex=cex),
                      points = allPoints,
                  confLimits = confLimits      
                      )

    if(Vi) addDistPlot(
                        d    = PP$mcmc[rInds,5,(Nburnin+1):PP$Nmcmc],
                        log  = F,
                        h    = PP$height[rInds],
                        xlim = ViLim,
                        ylim = hLim2,
                        xlab = list(expression(paste("V"[i]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex),
                      points = allPoints,
                  confLimits = confLimits      
                      )

    if(Ve) addDistPlot(
                        d    = (PP$mcmc[rInds,6,(Nburnin+1):PP$Nmcmc]+PP$mcmc[rInds,5,(Nburnin+1):PP$Nmcmc]),
                        log  = F,
                        h    = PP$height[rInds],
                        xlim = VeLim,
                        ylim = hLim2,
                        xlab = list(expression(paste("V"[e]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex),
                      points = allPoints,
                  confLimits = confLimits      
                      )

  for(k in seq(length(pIon))){
    if(pIon[k]){
      addDistPlot(
                        d    = PP$mcmc[rInds,(6+k),(Nburnin+1):PP$Nmcmc],
                        log  = F,
                        h    = PP$height[rInds],
                        xlim = pIonLim,
                        ylim = hLim2,
                        xlab =list(sprintf("Fractional abundance of ion mass %.1f u",PP$mi[k+1]),col=fg,cex=cex),
                      points = allPoints,
                  confLimits = confLimits      
                      )
    }

  }

  # if we did not plot on an x11 device, we must close the device properly
  if(sum(figList)==2) dev.off()

  invisible(PP)

} # plotMarginalDistribution

