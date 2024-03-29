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
##   multistatic Logical, readPP.3D is used for reading the result files if
##            multistatic is TRUE, otherwise readPP is used
##
##   All arguments except 'data' are optional
##


plotPP <- function(data,par=list(Ne=c(10,12),TeR1=c(0,4000),TiR1=c(0,3000),ViR1=c(-400,400)),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,multistatic=TRUE,trellis=FALSE,cutgaps=TRUE,main=NA,title=NA,measuredOnly=T,nSiteVi=3)
    {

        UseMethod("plotPP")

    }

plotPP.character <- function(data,par=list(Ne=c(10,12),TeR1=c(0,4000),TiR1=c(0,3000),ViR1=c(-400,400)),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,multistatic=TRUE,trellis=FALSE,cutgaps=TRUE,main=NA,title=NA,measuredOnly=T,nSiteVi=3)
    {

        # read all data
        if(multistatic){
            data <- readPP.3D(data,measuredOnly,nSiteVi)
        }else{
            data <- readPP(data)
        }

        args <- formals()
        argnames <- names(args)
        args <- lapply( names( args ) , FUN=function(x){ eval( as.name( x ) ) } )
        names(args) <- argnames

        do.call( plotPP , args )

    }


plotPP.list <- function(data,par=list(Ne=c(10,12),TeR1=c(0,4000),TiR1=c(0,3000),ViR1=c(-400,400)),xlim=NULL,ylim=NULL,pdf=NULL,jpg=NULL,figNum=NULL,width=8.27,height=2.9225,paper='a4',tickRes=NULL,model=F,stdThreshold=.5,NeMin=1e9,chisqrLim=10,bg='white',fg='black',res=300,cex=1.0,col.regions=guisdap.colors,multistatic=TRUE,trellis=FALSE,cutgaps=TRUE,main=NA,title=NA,...)
    {
        # take a copy that will be returned
        data2 <- data

        # if no data, return NULL
        if(is.null(data)){
            warning("No data")
            return(invisible(NULL))
        }

        # if the data matrices have only one row, we add NA rows above and below to allow plotting as usual
        if(dim(data[["param"]])[1] == 1 ){
            tmparr <- array(dim=dim(data[["param"]])+c(2,0,0))
            tmparr[2,,] <- data[["param"]]
            dimnames(tmparr) <- list(c('gate1','gate2','gate3'),dimnames(data[["param"]])[[2]],dimnames(data[["param"]])[[3]])
            data[["param"]] <- tmparr
            tmparr <- array(dim=dim(data[["std"]])+c(2,0,0))
            tmparr[2,,] <- data[["std"]]
            dimnames(tmparr) <- list(c('gate1','gate2','gate3'),dimnames(data[["std"]])[[2]],dimnames(data[["std"]])[[3]])
            data[["std"]] <- tmparr
            tmparr <- array(dim=dim(data[["model"]])+c(2,0,0))
            tmparr[2,,] <- data[["model"]]
            dimnames(tmparr) <- list(c('gate1','gate2','gate3'),dimnames(data[["model"]])[[2]],dimnames(data[["model"]])[[3]])
            data[["model"]] <- tmparr
            tmparr <- array(dim=dim(data[["height"]])+c(2,0))
            tmparr[2,] <- data[["height"]]
            tmparr[1,] <- data[["height"]] - 30
            tmparr[3,] <- data[["height"]] + 30
            data[["height"]] <- tmparr
        }

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

        # set  NA values to edges of data gaps
        if(cutgaps){
            td <- diff(data[["time_sec"]])
            tdmed <- median(td)
            tdlarge <- which(td>3*tdmed)
            data[["param"]][,,tdlarge] <- NA
            data[["param"]][,,(tdlarge+1)] <- NA
        }


        # number of figures panels
        nFig <- length(par)

        # time limits
        if(is.null(xlim)){
            tLim <- range(data[["time_sec"]])
        }else{
            tLim <- (floor(data$time_sec[1] / 3600 / 24) * 24 + xlim) * 3600
        }

        # data points within tLim
        tInds <- which(data[["time_sec"]]>=tLim[1] & data[["time_sec"]]<=tLim[2])

        if( length(tInds)==0){
            warning("No data from the given time period")
            return(invisible(data2))
        }

        # height limits
        if(is.null(ylim)){
            hLim <- range( data[["height"]] , na.rm=TRUE )
        }else{
            hLim <- ylim
        }


        # height of the full plot window
        wHeight <- min(nFig,4)*height

        if(length(tInds)==1){
            width <- ceiling(nFig/2)*width/2
        }
        # open the  proper device
        figList <- c(is.null(figNum),is.null(pdf),is.null(jpg))
        if(sum(figList) < 2 ) stop('Only one output device can be selected at a time')
        # a new x11 by defaul
        if(sum(figList) == 3) dev.new(width=width,height=wHeight)
        # new plot to an existing x11 window
        if(!is.null(figNum)) {dev.set(figNum);plot.new()}
        # a new pdf file
        if(!is.null(pdf)) pdf(file=paste(pdf,'.pdf',sep=''),paper=paper,width=width,height=wHeight,bg=bg)
        # a new jpeg file
        if(!is.null(jpg)){
            jpeg(file=paste(jpg,'.jpg',sep=''),width=width,height=wHeight,units='in',res=res)
            cex = cex*res/72
        }

        if(length(tInds)==1){
            layout( matrix(seq(ceiling(length(par)/2)*2) , nrow=2 ) )
            par(mar=c(3,3,3,1)*cex,mgp=c(1.5,.5,0)*cex)
            trellis <- FALSE
        }else{
            if(trellis){
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
            }else{
                layout(matrix(seq(2*length(par)+4),ncol=2,byrow=T),widths=c(.9/sqrt(cex),.1*sqrt(cex)),heights=c(.2,rep(1,length(par)),.4))
                par(mar=c(0,3,0,1)*cex,mgp=c(1.5,.5,0)*cex)
                plot.new()
                plot.new()
            }
            # tick marks in the time axis
            ticks <- timeTicks(tLim,tickRes)
            # tick marks for the height axis
            hticks <- heightTicks(hLim)

        }

        # actual plotting
        curFig <- 1

        maincpy <- main

        for(p in seq(length(par))){

            # select whether x-axis is logarithmic or not
            xlog <- switch( names(par[p]) , Ne=TRUE , Coll=TRUE , FALSE )

            # plot titles, there are spare ones for future systems.. :)
            if(is.na(maincpy[p])){
                main <- switch( names(par[p]) ,
                               Ne=list(expression(paste("Log","(N"[e]^{},") [m"[]^{-3},"]")),col=fg,cex=cex,vjust=1.2),
                               Ti=list(expression(paste("T"[i]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                               Tipar=list(expression(paste("T"["i||"]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                               Tiperp=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                               Tiratio=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"/","T"["i||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                               Teratio=list(expression(paste("T"[paste("e",symbol("\136"))]^{},"/","T"["e||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                               Tidiff=list(expression(paste("T"[paste("i",symbol("\136"))]^{},"-","T"["i||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                               Tediff=list(expression(paste("T"[paste("e",symbol("\136"))]^{},"-","T"["e||"]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                               Te=list(expression(paste("T"[e]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                               Tepar=list(expression(paste("T"["e||"]^{}," [K]")),col=fg,cex=cex,vjust=1.2),
                               Teperp=list(expression(paste("T"[paste("e",symbol("\136"))]^{},"[K]")),col=fg,cex=cex,vjust=1.2),
                               Coll=list(expression(paste("Log(",nu['in']^{},") [s"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               Vix=list(expression(paste("V"[i]^{}," East [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               Viy=list(expression(paste("V"[i]^{}," North [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               Viz=list(expression(paste("V"[i]^{}," Up [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViB=list(expression(paste("V"["i||"]^{}," Up [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViBx=list(expression(paste("V"[paste("i",symbol("\136"))]^{}," East [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViBy=list(expression(paste("V"[paste("i",symbol("\136"))]^{}," North [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR1=list(expression(paste("V"[iR1]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR2=list(expression(paste("V"[iR2]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR3=list(expression(paste("V"[iR3]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR4=list(expression(paste("V"[iR4]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR5=list(expression(paste("V"[iR5]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR6=list(expression(paste("V"[iR6]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR7=list(expression(paste("V"[iR7]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR8=list(expression(paste("V"[iR8]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR9=list(expression(paste("V"[iR9]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR10=list(expression(paste("V"[iR10]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR11=list(expression(paste("V"[iR11]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR12=list(expression(paste("V"[iR12]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR13=list(expression(paste("V"[iR13]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR14=list(expression(paste("V"[iR14]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR15=list(expression(paste("V"[iR15]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR16=list(expression(paste("V"[iR16]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR17=list(expression(paste("V"[iR17]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR18=list(expression(paste("V"[iR18]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR19=list(expression(paste("V"[iR19]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR20=list(expression(paste("V"[iR20]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR21=list(expression(paste("V"[iR21]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR22=list(expression(paste("V"[iR22]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR23=list(expression(paste("V"[iR23]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR24=list(expression(paste("V"[iR24]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR25=list(expression(paste("V"[iR25]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR26=list(expression(paste("V"[iR26]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR27=list(expression(paste("V"[iR27]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR28=list(expression(paste("V"[iR28]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR29=list(expression(paste("V"[iR29]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),

                               ViR1hor=list(expression(paste("V"[iR1hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR2hor=list(expression(paste("V"[iR2hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR3hor=list(expression(paste("V"[iR3hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR4hor=list(expression(paste("V"[iR4hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR5hor=list(expression(paste("V"[iR5hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR6hor=list(expression(paste("V"[iR6hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR7hor=list(expression(paste("V"[iR7hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR8hor=list(expression(paste("V"[iR8hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR9hor=list(expression(paste("V"[iR9hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR10hor=list(expression(paste("V"[iR10hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR11hor=list(expression(paste("V"[iR11hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR12hor=list(expression(paste("V"[iR12hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR13hor=list(expression(paste("V"[iR13hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR14hor=list(expression(paste("V"[iR14hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR15hor=list(expression(paste("V"[iR15hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR16hor=list(expression(paste("V"[iR16hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR17hor=list(expression(paste("V"[iR17hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR18hor=list(expression(paste("V"[iR18hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR19hor=list(expression(paste("V"[iR19hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR20hor=list(expression(paste("V"[iR20hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR21hor=list(expression(paste("V"[iR21hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR22hor=list(expression(paste("V"[iR22hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR23hor=list(expression(paste("V"[iR23hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR24hor=list(expression(paste("V"[iR24hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR25hor=list(expression(paste("V"[iR25hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR26hor=list(expression(paste("V"[iR26hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR27hor=list(expression(paste("V"[iR27hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR28hor=list(expression(paste("V"[iR28hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),
                               ViR29hor=list(expression(paste("V"[iR29hor]^{}," (away) [ms"[]^{-1},"]")),col=fg,cex=cex,vjust=1.2),

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
            }else{
                main <- maincpy[p]
            }

            # azel fit removed, this is not trivial with several transmitters...


#            # azimuth and elevation in 1D fits
            if(names(par[p])=="azel"){
                stop("Sorry, the azel fit option has been removed.")
#                if(par[[p]]){
#                    if(trellis){
#                        trenew$layout.widths$right.padding <- 10.5
#                        trenew$layout.heights$key.axis.padding<--1
#                        trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights,layout.widths=trenew$layout.widths))
#                        curFig <- addAzElPlot(
#                            az     = data$azelT[,1,tInds]%%360,
#                            el     = data$azelT[,2,tInds],
#                            t      = data$time_sec[tInds],
#                            xlim   = tLim,
#                            cex    = cex,
#                            ticks  = ticks,
#                            nFig   = nFig,
#                            curFig = curFig
#                            )
#                        trenew$layout.widths$right.padding <- 1
#                        trenew$layout.heights$key.axis.padding<-1
#                        trellis.par.set(list(background=trenew$background,layout.heights=trenew$layout.heights,layout.widths=trenew$layout.widths))
#                    }else{
#                        plot(data$time_sec[tInds],data$azelT[,1,tInds]%%360,ylim=c(0,360),xlim=tLim,xaxt='n',xlab='',ylab='Degrees')
#                        points(data$time_sec[tInds],data$azelT[,2,tInds],col='red')
#                        axis(1,at=ticks$tick,labels=ticks$string)
#                        plot.new()
#                    }
#                }
            }else{
                d <- data[["param"]][,names(par[p]),tInds]
                if(length(tInds)==1){
                    err <- data[["std"]][,names(par[p]),tInds]
                    if( xlog ){
                        d[data[["std"]][,names(par[p]),tInds]>10**(par[[p]][3])] <- NA
                        err[data[["std"]][,names(par[p]),tInds]>10**(par[[p]][3])] <- NA
                    }else{
                        d[data[["std"]][,names(par[p]),tInds]>par[[p]][3]] <- NA
                        err[data[["std"]][,names(par[p]),tInds]>par[[p]][3]] <- NA
                    }
                    addPPlinePlot( d=d , err=err , log=xlog , h=data[["height"]][,tInds] , xlim=par[[p]][1:2] , ylim=hLim , xlab=main[[1]] , cex=cex )
                }else{
                    if( xlog ){
                        d <- log10( d )
                        d[data[["std"]][,names(par[p]),tInds]>10**(par[[p]][3])] <- NA
                    }else{
                        d[data[["std"]][,names(par[p]),tInds]>par[[p]][3]] <- NA
                    }
                    d[d<par[[p]][1]] <- par[[p]][1]
                    d[d>par[[p]][2]] <- par[[p]][2]
                    if(trellis){
                        curFig <- addPPcolorPlot(
                            d = d,
                            h = data[["height"]][,tInds],
                            t = colMeans(data[["timeLimits"]][,tInds]),
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
                    }else{
                        if (p==length(par)){
                            tickstr <- ticks$string
                        }else{
                            tickstr <- rep('',length(ticks$tick))
                        }
                        image(colMeans(data$timeLimits[,tInds]),data$height[,tInds[1]],t(d),xlim=tLim,ylim=hLim,zlim=par[[p]][1:2],col=col.regions(1000),xlab='',xaxt='n',ylab='Height [km]',cex=cex,cex.lab=cex,cex.axis=cex,yaxt='n')
                        axis(1,at=ticks$tick,labels=tickstr,cex=cex,cex.lab=cex,cex.axis=cex)
                        axis(2,at=hticks$tick,labels=hticks$tick,cex=cex,cex.lab=cex,cex.axis=cex)
                        marold <- par()$mar
                        par(mar=c(1,3,1,1)*cex)
                        image(c(0,1),seq(par[[p]][1],par[[p]][2],length.out=1000),t(matrix(rep(seq(par[[p]][1],par[[p]][2],length.out=1000),2),ncol=2)),col=col.regions(1000),ylab=main[[1]],xaxt='n',xlab='',cex=cex,cex.lab=cex,cex.axis=cex)
                        par(mar=marold)
                    }
                }
            }
        }
        if(!trellis){
            if(length(tInds)>1){
                plot.new()
                title(xlab=list("UTC",cex=cex),line=-2,outer=FALSE)
                plot.new()
            }
            if(is.na(title)){
                if(length(tInds)==1){
                    mtext( paste( as.character( data[["POSIXtime"]][[tInds]] ) , "UTC" ), side = 3, line = -2, outer = TRUE , cex = cex )
                }else{
                    mtext( substr( as.character( data[["POSIXtime"]][[tInds[1]]] ) , 1 , 10 ) , side = 3, line = -1.5, outer = TRUE , cex = cex )
                }
            }else{
                mtext( title , side = 3, line = -1.5, outer = TRUE , cex = cex )
            }
        }
        # if we did not plot on an x11 device, we must close the device properly
        if((sum(figList)==2)&is.null(figNum)) dev.off()

        return(invisible(data2))

    }


