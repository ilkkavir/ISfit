plotMarginalDistribution <- function(fpath,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',par=list(Ne=c(10,12),Tepar=c(0,4000),Teperp=c(0,4000),Tipar=c(0,3000),Tiperp=c(0,3000),Vix=c(-1000,1000),Viy=c(-1000,1000),Viz=c(-400,400)),hLim=NULL,tLim=NULL,bg='white',fg='black',res=300,cex=1.0,allPoints=F,confLimits=c(.682,.954,.996)){
#
# Plot the marginal posterior distributions from an MCMC analysis.
# Either individual points (X11 only to avoid large files) or color-coded probability regions
# The possibility for large files with individual points plotted is easy to add, if someone really needs that
#
# I. Virtanen 2010
#

  # load the data
    if(!file.exists(fpath)) stop(paste("Path",fpath,"does not exist"))
    if(file.info(fpath)$isdir) stop("Only contents of individual files can be plotted")
    load(fpath)


  # find the parameters the user wants to plot
  rInds <- rep(T,length(PP$height))

  hLim2 <- range(PP$height,na.rm=T)
  hInds <- seq(length(PP$height))
  if(!is.null(hLim)){ hInds <- (PP$height >= hLim[1]) & (PP$height <= hLim[2]); hLim2 <- hLim}

  if(!any(hInds)) stop('No data from the given height interval')

  # how many frames do we actuallly have?
  nFig <- length(par)

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

  # actual plotting
        for(p in seq(length(par))){

            # select whether x-axis is logarithmic or not
            xlog <- switch( names(par[p]) , Ne=TRUE , Coll=TRUE , FALSE )

            # x-axis labels
            xlab <- switch( names(par[p]) ,
                           Ne=list(expression(paste("Log","(N"[e]^{},") [m"[]^{-3},"]")),col=fg,cex=cex,vjust=1.2),
                           Ti=list(expression(paste("T"[i]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           Tipar=list(expression(paste("T"["i||"]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           Tiperp=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                           Tiratio=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"/","T"["i||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                           Teratio=list(expression(paste("T"[paste("e",symbol("\136"))]^{},"/","T"["e||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                           Te=list(expression(paste("T"[e]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           Tepar=list(expression(paste("T"["e||"]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           Teperp=list(expression(paste("T"[paste("e",symbol("\136"))]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                           Coll=list(expression(paste("Log(",nu['in']^{},") [s"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           Vix=list(expression(paste("V"[i]^{}," East [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           Viy=list(expression(paste("V"[i]^{}," North [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           Viz=list(expression(paste("V"[i]^{}," Up [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViB=list(expression(paste("V"["i||"]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViBx=list(expression(paste("V"[paste("i",symbol("\136"))]^{}," East [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViBy=list(expression(paste("V"[paste("i",symbol("\136"))]^{}," North [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
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
                           ViR11=list(expression(paste("V"[iR11]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR12=list(expression(paste("V"[iR12]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR13=list(expression(paste("V"[iR13]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR14=list(expression(paste("V"[iR14]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR15=list(expression(paste("V"[iR15]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR16=list(expression(paste("V"[iR16]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR17=list(expression(paste("V"[iR17]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR18=list(expression(paste("V"[iR18]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR19=list(expression(paste("V"[iR19]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR20=list(expression(paste("V"[iR20]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR21=list(expression(paste("V"[iR21]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR22=list(expression(paste("V"[iR22]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR23=list(expression(paste("V"[iR23]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR24=list(expression(paste("V"[iR24]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR25=list(expression(paste("V"[iR25]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR26=list(expression(paste("V"[iR26]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR27=list(expression(paste("V"[iR27]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR28=list(expression(paste("V"[iR28]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR29=list(expression(paste("V"[iR29]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),

                           ViR1hor=list(expression(paste("V"[iR1hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR2hor=list(expression(paste("V"[iR2hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR3hor=list(expression(paste("V"[iR3hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR4hor=list(expression(paste("V"[iR4hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR5hor=list(expression(paste("V"[iR5hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR6hor=list(expression(paste("V"[iR6hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR7hor=list(expression(paste("V"[iR7hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR8hor=list(expression(paste("V"[iR8hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR9hor=list(expression(paste("V"[iR9hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR10hor=list(expression(paste("V"[iR10hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR11hor=list(expression(paste("V"[iR11hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR12hor=list(expression(paste("V"[iR12hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR13hor=list(expression(paste("V"[iR13hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR14hor=list(expression(paste("V"[iR14hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR15hor=list(expression(paste("V"[iR15hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR16hor=list(expression(paste("V"[iR16hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR17hor=list(expression(paste("V"[iR17hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR18hor=list(expression(paste("V"[iR18hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR19hor=list(expression(paste("V"[iR19hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR20hor=list(expression(paste("V"[iR20hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR21hor=list(expression(paste("V"[iR21hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR22hor=list(expression(paste("V"[iR22hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR23hor=list(expression(paste("V"[iR23hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR24hor=list(expression(paste("V"[iR24hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR25hor=list(expression(paste("V"[iR25hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR26hor=list(expression(paste("V"[iR26hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR27hor=list(expression(paste("V"[iR27hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR28hor=list(expression(paste("V"[iR28hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                           ViR29hor=list(expression(paste("V"[iR29hor]^{}," [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),

                           TiR1=list(expression(paste("T"[iR1]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR2=list(expression(paste("T"[iR2]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR3=list(expression(paste("T"[iR3]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR4=list(expression(paste("T"[iR4]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR5=list(expression(paste("T"[iR5]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR6=list(expression(paste("T"[iR6]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR7=list(expression(paste("T"[iR7]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR8=list(expression(paste("T"[iR8]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR9=list(expression(paste("T"[iR9]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR10=list(expression(paste("T"[iR10]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR11=list(expression(paste("T"[iR11]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR12=list(expression(paste("T"[iR12]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR13=list(expression(paste("T"[iR13]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR14=list(expression(paste("T"[iR14]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR15=list(expression(paste("T"[iR15]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR16=list(expression(paste("T"[iR16]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR17=list(expression(paste("T"[iR17]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR18=list(expression(paste("T"[iR18]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR19=list(expression(paste("T"[iR19]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR20=list(expression(paste("T"[iR20]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR21=list(expression(paste("T"[iR21]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR22=list(expression(paste("T"[iR22]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR23=list(expression(paste("T"[iR23]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR24=list(expression(paste("T"[iR24]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR25=list(expression(paste("T"[iR25]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR26=list(expression(paste("T"[iR26]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR27=list(expression(paste("T"[iR27]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR28=list(expression(paste("T"[iR28]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TiR29=list(expression(paste("T"[iR29]^{}," [K]")),col=fg,cex=cex,vjust=1.2),

                           TeR1=list(expression(paste("T"[eR1]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR2=list(expression(paste("T"[eR2]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR3=list(expression(paste("T"[eR3]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR4=list(expression(paste("T"[eR4]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR5=list(expression(paste("T"[eR5]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR6=list(expression(paste("T"[eR6]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR7=list(expression(paste("T"[eR7]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR8=list(expression(paste("T"[eR8]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR9=list(expression(paste("T"[eR9]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR10=list(expression(paste("T"[eR10]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR11=list(expression(paste("T"[eR11]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR12=list(expression(paste("T"[eR12]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR13=list(expression(paste("T"[eR13]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR14=list(expression(paste("T"[eR14]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR15=list(expression(paste("T"[eR15]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR16=list(expression(paste("T"[eR16]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR17=list(expression(paste("T"[eR17]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR18=list(expression(paste("T"[eR18]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR19=list(expression(paste("T"[eR19]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR20=list(expression(paste("T"[eR20]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR21=list(expression(paste("T"[eR21]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR22=list(expression(paste("T"[eR22]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR23=list(expression(paste("T"[eR23]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR24=list(expression(paste("T"[eR24]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR25=list(expression(paste("T"[eR25]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR26=list(expression(paste("T"[eR26]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR27=list(expression(paste("T"[eR27]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR28=list(expression(paste("T"[eR28]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                           TeR29=list(expression(paste("T"[eR29]^{}," [K]")),col=fg,cex=cex,vjust=1.2),

                           Ion1=list(sprintf("Fraction of ion mass %.1f u",data$mi[1]),col=fg,cex=cex,vjust=1.2),
                           Ion2=list(sprintf("Fraction of ion mass %.1f u",data$mi[2]),col=fg,cex=cex,vjust=1.2),
                           Ion3=list(sprintf("Fraction of ion mass %.1f u",data$mi[3]),col=fg,cex=cex,vjust=1.2),
                           Ion4=list(sprintf("Fraction of ion mass %.1f u",data$mi[4]),col=fg,cex=cex,vjust=1.2),
                           Ion5=list(sprintf("Fraction of ion mass %.1f u",data$mi[5]),col=fg,cex=cex,vjust=1.2)
                           )


            # need to somehow combine list elements from different heights.. the MCMC output is not a huge array anymore!
            nheights <- length(hInds)
            nmcmc <- PP$functionCall$MCMCsettings$niter - PP$functionCall$MCMCsettings$burninlength
            d <- matrix(NA,nrow=nheights,ncol=nmcmc)
            for(rii in seq(nheights)){
                if(!is.null(PP[["MCMC"]][[hInds[rii]]])){
                    d[ rii , ] <- PP[["MCMC"]][[hInds[rii]]][["pars"]][ , names(par)[p]]
                }
            }

            addDistPlot( d=d , log=xlog , h=PP$height[hInds] , xlim=par[[p]] , ylim=hLim2 , xlab=xlab[[1]] , points=allPoints , confLimits=confLimits )

        }


  # if we did not plot on an x11 device, we must close the device properly
  if(sum(figList)==2) dev.off()

  invisible(PP)

} # plotMarginalDistribution

