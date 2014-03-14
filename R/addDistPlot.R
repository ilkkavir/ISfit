addDistPlot <- function(d,log=F,h,xlim,ylim=range(h,na.rm=T),xlab,points,confLimits){
#
# plot of the posteriori distributions from an MCMC analysis
#
# I. Virtanen 2010
#

  if(log) d <- log10(d)

  n <- length(d[1,])

  if(points){
    # plot the first range and set the axis limits
    plot(d[1,],rep(h[1],n),pch='.',xlim=xlim,ylim=ylim,xlab=xlab,ylab='Height [km]')

    # plot all the remaining ranges
    for(k in seq(1,length(h))){
      points(d[k,],rep(h[k],n),pch='.')
    }

  }else{

    # polygons around the confidence areas
    confPolygons <- confidencePolygons(d,h,confLimits,xlim + c(-1,1)*10*max(abs(xlim)))

    # an idiotic but simple way to set the proper margins
    plot(c(0,0),col='transparent',xlim=xlim,ylim=ylim,xlab=xlab,ylab='Height [km]')

    # plot the polygons, the list contains both the limits and colors
    for(k in seq(length(confPolygons))){
      polygon(confPolygons[[k]]$x,confPolygons[[k]]$y,col=confPolygons[[k]]$col,border=NA)
    }

    # line of the median value
    lines(apply(d,FUN=median,MARGIN=c(1)),h,lwd=2,col=rgb(0,0,0))

  }

} # addDistPlot

