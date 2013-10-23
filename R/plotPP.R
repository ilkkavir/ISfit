## file:plotPP.R
## (c) 2010- University of Oulu, Finland
## Written by Ilkka Virtanen <ilkka.i.virtanen@oulu.fi>
## Licensed under FreeBSD license.
##

##
## Color-coded image of plasma parametr profiles as function of time and height
##
## Arguments:
##
##   data     either a character vector of data files and/or directories, or an output list
##            of either PPplot, readPP, or readPP.3D
##   par      a list of parameter names, corresponding z-axis limits, and limits for standard deviations
##   xlim     x-axis limits. In decimal UT hours with zero in the beginning of the day from which
##            first data point is found. Use values larger than 24 to cross the midnight
##   ylim     y-axis (height) limits in km
##   jpg      file name for optional jpg file output
##   figNum   number of a graphics device to plot on. Used for avoiding opening numerous devices in
##            repeated calls.
##   width    width of a single parameter panel in  the plot in inches
##   height   height of a single parameter panel
##   paper    paper type selection
##   tickRes  time tick mark resolution in minutes
##   model    logical, plot model values instead of the actual data
##   stdThreshold threshold for relative error in electron density.
##            All data points at which std(Ne)/Ne > stdThreshold are replaced  with NA
##   NeMin    lower limit for electron density, works in the same way as stdThreshold
##   chisqrLim Limit for normalized residual
##   bg       graphics device background color
##   fg       graphics device foreground color
##   res      resolution in jpg plot
##   cex      magnification in axis labels etc.
##   col.regions color palette selection
##   ytype    y-axis type, 'height' in 3D analysis, 'height' or 'range' in 1D
##   multistatic Logical, readPP.3D is used for reading the result files if
##            multistatic is TRUE, otherwise readPP is used
##
##   All arguments except 'data' are optional
##


plotPP <- function(data,par=list(Ne=c(10,12),Ti=c(0,3000),Te=c(0,4000),ViR1=c(-200,200)),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,ytype='height',multistatic=TRUE)
    {

        UseMethod("plotPP")

    }

plotPP.character <- function(data,par=list(Ne=c(10,12),Ti=c(0,3000),Te=c(0,4000),ViR1=c(-200,200),azel=T),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,ytype='height',multistatic=TRUE)
    {

        # read all data
        if(multistatic){
            data <- readPP.3D(data)
        }else{
            data <- readPP(data)
        }

        args <- formals()
        argnames <- names(args)
        args <- lapply( names( args ) , FUN=function(x){ eval( as.name( x ) ) } )
        names(args) <- argnames
    
        do.call( plotPP , args )

    }
    

plotPP.list <- function(data,par=list(Ne=c(10,12),Ti=c(0,3000),Te=c(0,4000),ViR1=c(-200,200),azel=T),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,ytype='height',multistatic=TRUE)
    {
        # take a copy that will be returned
        data2 <- data

        # if no data, return NULL
        if(is.null(data)) invisible(NULL)
  
        # cut off large variances (only Ne is studied)
        if(!is.null(stdThreshold)){
            thrInds <- (data$std[,1,] / data$param[,1,]) >  stdThreshold
            for(k in seq(dim(data[["param"]])[2])) data$param[,k,][thrInds] <- NA
        }

        # cut off all data where Ne is below NeMin
        if(!is.null(NeMin)){
            NeInds <- data$param[,1,] < NeMin
            for(k in seq(dim(data$param)[2])) data$param[,k,][NeInds] <- NA
        }

        # cut off all failed iterations
        statInds <- data$status != 0
        for(k in seq(dim(data$param)[2])) data$param[,k,][statInds] <- NA

        # cut off very large chi-squared values
        if(!is.null(chisqrLim)){
            chiInds <- data$chisqr >  chisqrLim
            for(k in seq(dim(data$param)[2])) data$param[,k,][chiInds] <- NA
        }

        # plot the initial(model) parameters instead of the fit results
        if(model){
            data$param   <- data$model
            data$std[,,] <- 0
        }

        # number of figures panels
        nFig <- length(par)

        # time limits
        if(is.null(xlim)){
            tLim <- range(data[["time_sec"]])
        }else{
            tLim <- (floor(data$time_sec[1] / 3600 / 24) * 24 + xlim) * 3600
        }

        # height limits
        if(is.null(ylim)){
            hLim <- range(data[["height"]])
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
                           Tiratio=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"/","T"["i||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
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
                           Ion1=list(sprintf("Fraction of ion mass %.1f u",data$mi[1]),col=fg,cex=cex,vjust=1.2),
                           Ion2=list(sprintf("Fraction of ion mass %.1f u",data$mi[2]),col=fg,cex=cex,vjust=1.2),
                           Ion3=list(sprintf("Fraction of ion mass %.1f u",data$mi[3]),col=fg,cex=cex,vjust=1.2),
                           Ion4=list(sprintf("Fraction of ion mass %.1f u",data$mi[4]),col=fg,cex=cex,vjust=1.2),
                           Ion5=list(sprintf("Fraction of ion mass %.1f u",data$mi[5]),col=fg,cex=cex,vjust=1.2)
                           )
            # azimuth and elevation in 1D fits
            if(names(par[p])=="azel"){
                if(par[[p]]){
                    trenew$layout.widths$right.padding <- 10.5
                    trenew$layout.heights$key.axis.padding<--1
                    trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights,layout.widths=trenew$layout.widths))
                    curFig <- addAzElPlot(
                        az     = data$azelT[,1,]%%360,
                        el     = data$azelT[,2,],
                        t      = data$time_sec,
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
                d <- data[["param"]][,names(par[p]),]
                if( xlog ){
                    d <- log10( d )
                    d[data[["std"]][,names(par[p]),]>10**(par[[p]][3])] <- NA
                }else{
                    d[data[["std"]][,names(par[p]),]>par[[p]][3]] <- NA
                }                
                curFig <- addPPcolorPlot(
                    d = d,
                    h = data[["height"]],
                    t = data[["time_sec"]],
                    xlim = tLim,
                    ylim = hLim,
                    zlim = par[[p]][1:2],
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

        invisible(data2)

    }


