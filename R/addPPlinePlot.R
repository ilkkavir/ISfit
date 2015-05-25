addPPlinePlot <- function(d,err,log=F,h,xlim,ylim=range(h,na.rm=T),xlab , cex=1 ){
# OK, actually this should be addPPpointplot.. :)
#
# plots the data points as function of height, plus the given errorbars as lines
#
# I. Virtanen 2010
#

  # lower limits of the errorbars
  errLims1 <- d-err
  # upper limits of the errorbars
  errLims2 <- d+err

  # if a logarithmic x-axis is used
  if(log){

    # because the x-axis is logarithmic, negative values are not allowed
    errLims1 <- pmax(errLims1,0)
    errLims2 <- pmax(errLims2,0)

    # take a 10-based logarithm of the errorbars
    errLims1 <- log10(errLims1)
    errLims2 <- log10(errLims2)

    # now the originally negative error limits have the value -Inf, replace them with a negative value smaller than the x-axis limit
    errLims1[is.infinite(errLims1)] <- -2*max(xlim)

    # a similar treatment for +Inf values, though they were not produced here
    errLims2[is.infinite(errLims2)] <- 2*max(xlim)

    # plot the logarithm of the actual data
    plot(log10(d),h,xlim=xlim,xlab=xlab,ylab='Height [km]',ylim=ylim,cex.axis=cex,cex.lab=cex)

    # plot the errorbars as lines
    arrows(errLims1,h,errLims2,h,code=3,length=0)

  # a linear scale
  }else{

    # plot the data as points
    plot(d,h,xlim=xlim,xlab=xlab,ylab='Height [km]',ylim=ylim,cex.axis=cex,cex.lab=cex)

    # add the errorbars as lines
    arrows(errLims1,h,errLims2,h,code=3,length=0)

  }

} # addPPlinePlot



