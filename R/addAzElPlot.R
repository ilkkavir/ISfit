addAzElPlot <- function(az,el,t,xlim,ticks,cex,nFig,curFig){
# 
# add a plot with antenna azimuth and elevation
# 
# 
# I. Virtanen 2013
# 

  # the x-axis label 'UT' is added only on the bottom panel
  xlab = ''
  if(nFig==curFig) xlab = 'UT'

  # a separate print command is needed with xy to actually print the plot the current device
  print(

    # form the grid graphics object
    xyplot(
      az+el       ~ t,
      data        = data.frame(az=az,el=el),
      scales      = list(x=list(at=ticks$tick,labels=ticks$string),cex=cex),
      xlab        = list(xlab,cex=cex),
      ylab        = list("Degrees",cex=cex),
      xlim        = xlim,
      ylim        = c(0,360),
      pch         = 20,
      col         = c('black','red'),
      key=list(space="top",text=list(c('TX azimuth','TX elevation')),points=list(pch=20,col=c('black','red')),columns=2)
    ),

    split=c(1,curFig,1,nFig),
    more=T

  )
  
  return((curFig+1))

} #addAzElPlot

