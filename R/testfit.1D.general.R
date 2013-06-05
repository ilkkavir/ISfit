testfit.1D.general <- function(h=200,fradar=933e6,plotTest=F,diffLimit=1e-2,absLimit=5,nacf=50,nlag=50,maxLambda=1e30){
#
# Test the 1D plasma parameter fit with simulated ACF data
#
#
#
#
#
#
#
#
#
#

  # parameters from iri model
  ptmp           <- iriParams()

  # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
  # the densities in outfmsis are in cm^-3
#  ioncoll        <- (2.44e-16*ptmp$outfmsis[2,h] + 4.34e-16*ptmp$outfmsis[3,h] + 4.28e-16*ptmp$outfmsis[4,h])*1e6
  ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp[,h])['NO+',] )

  # an approximation for electron-neutral collision frequency
  elecoll        <- ioncoll*.35714

  # electron parameters
#  ele            <- c(ptmp$outf[c(1,4,4),h],elecoll,0,0,0)
  ele            <- c(ptmp[c('e-','Te','Te'),h],elecoll,0,0,0)

  # ion parameters
# ion             <- list(
#                         c(30.5,(ptmp$outf[8,h]+ptmp$outf[9,h]),rep(ptmp$outf[3,h],2),ioncoll,0,0,0),
#                         c(16.0,ptmp$outf[5,h],rep(ptmp$outf[3,h],2),ioncoll,0,0,0),
#                         c(1.0,ptmp$outf[6,h],rep(ptmp$outf[3,h],2),ioncoll,0,0,0)
#                         )

 ion             <- list(
                         c(30.5,(ptmp['O2+',h]+ptmp['NO+',h])/ele[1],ptmp[c('Ti','Ti'),h],ioncoll,0,0,0),
                         c(16.0,ptmp[c('O+'),h]/ele[1],ptmp[c('Ti','Ti'),h],ioncoll,0,0,0),
                         c(1.0,ptmp[c('H+'),h]/ele[1],ptmp[c('Ti','Ti'),h],ioncoll,0,0,0)
                         )

  nIon           <- length(ion)

  # simulated ACF
  ndata          <- nacf*nlag
  freq           <- seq(-10000,10000)*10
  lags           <- seq(nlag)*10e-6
  simudata       <- rep(simuACF(ele=ele,ion=ion,fradar=fradar,freq=freq,lags=lags),nacf) + (rnorm(ndata) + 1i*rnorm(ndata))*1e-20/sqrt(2)
  plot(Re(simuACF(ele=ele,ion=ion,fradar=fradar,freq=freq,lags=lags)))

  # variance
  simuvar        <- rep(1e-40,ndata)

  # parameters in the form used by ISparamfit
  par            <- ISparamList2Vec.general(ele,ion,c(1))

  # parameter scaling factors
  parScales      <- ISparamScales.general(par,nIon)

  # parameter value limits
  parLimits      <- ISparamLimits.general(nIon,1)

  # copy of the original parameters
  parcopy        <- par


  # scale the initial parameter values
  initParam      <- scaleParams( par , parScales , inverse=F)

  # scale the parameter limits
  limitParam     <- parLimits
  limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
  limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)
  
  # apriori information
  apriori        <- ISapriori.1D.general( initParam , nIon )

  initParam[1] <- max(abs(simudata),na.rm=TRUE)*1e18
#  initParam[9] <- 2
  initParam[2] <- 3

  print(system.time(fitpar <- ISparamfit.general( acf=simudata,
                                var=simuvar,
                                nData=length(simudata),
                                iSite=1,
                                fSite=fradar,
                                aSite=180,
                                initParam=initParam,
                                invAprioriCovar=apriori$invAprioriCovar,
                                aprioriTheory=apriori$aprioriTheory,
                                aprioriMeas=apriori$aprioriMeas,
                                nIon=3,
                                paramLimits=limitParam,
                                directTheory=ISdirectTheory.general,
                                absLimit=absLimit,
                                diffLimit=diffLimit,
                                maxIter=100,
                                scaleFun=scaleParams,
                                scale=parScales,
                                lags=lags,
                                plotTest=plotTest,
                                B=c(1,0,0),
                                kSite=list(c(1,0,0)),
                                maxLambda=maxLambda
                                )))



  parfit         <- scaleParams(fitpar$param,scale=parScales,inverse=T)
  fitstd         <- scaleParams(sqrt(diag(fitpar$covar)),scale=parScales,inverse=T)

  for (k in seq(length(parcopy))){
    cat(sprintf("%25.10f %25.10f %25.10f %25.10f \n",parcopy[k],parfit[k],fitstd[k],parfit[k]/parcopy[k]))
  }
  invisible(fitpar)
} # testfit.1D.general
