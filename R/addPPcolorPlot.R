addPPcolorPlot <- function(d,h,t,xlim=range(t),ylim=range(h),zlim,main,ticks,cex,nFig,curFig,col.regions=guisdap.colors){
# 
# add a pseudo-color plot (or what is it called?!?) of the given data
# 
# 
# I. Virtanen 2010
# 


  # make a proper grid
  grid <- expand.grid(x=t,y=h)
  # add the data to the grid (there must be another (faster) way to use grid graphics..)
  grid$z <- as.vector(t(d))

  # the x-axis label 'UT' is added only on the bottom panel
  xlab = ''
  if(nFig==curFig) xlab = 'UT'

  # a separate print command is needed with levelplot to actually print the plot the current device
  print(

    # form the grid graphics object
    levelplot(
      z           ~ x*y,
                    grid,
      col.regions = col.regions,
      at          = seq(zlim[1],zlim[2],length.out=100),
      scales      = list(x=list(at=ticks$tick,labels=ticks$string),cex=cex),
      xlab        = list(xlab,cex=cex),
      ylab        = list("Height [km]",cex=cex),
      xlim        = xlim,
      ylim        = ylim,
      main        = main,
      colorkey    = list(labels=list(cex=cex))
    ),

    split=c(1,curFig,1,nFig),
    more=T

  )
  
  return((curFig+1))

} #addPPcolorPlot

