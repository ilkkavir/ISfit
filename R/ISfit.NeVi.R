ISfit.NeVi <- function(ddir='.' , odir='.' ,  heightLimits.km=NA , timeRes.s=1 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=10 , plotFit=FALSE , recursive=TRUE , calScale=5e-20, nCores=1 ){
      #
      # A simple fit of Ne and Vi on monostatic E-region data. This function calls ISfit.3D, which is the more general multistatic
      # fitting tool.
      #
      # INPUT:
      #   ddir            Data directory
      #   odir            Output directory
      #   heightLimits.km analysis height-gate limits. If NA, range gates of the reference site are used
      #   timeRes.s       time resolution (integration time)
      #   beginTime       c(year,month,day,hour,minute,seconds) analysis start time
      #   endTime         c(year,month,day,hour,minute,seconds) analysis end time
      #   absLimit        limit for absolute value of the residual.
      #                   The iteration will not be stopped (unles maxIter is reached) before the residual is below absLimit
      #   diffLimit       Upper limit for fractional change in residual in an iteration step.
      #   maxLambda       maximum Lambda value in Levenberg-Marquardt iteration
      #   maxIter         maximum number of iterations
      #   plotFit         plot measured and fitted ACF for each fit
      #   recursive       logical, should the data directories be searched recursively
      #   calScale        additional scaling factor from ionosonde calibration applied to ALL ACF samples
      #   nCores          number of parallel processes
      #
      # OUTPUT:
      #   None, the results are written to files in odir.
      #
      #
      # IV 2017
      #

    # call ISfit.3D with sufficient parameters. This produces
    ISfit.3D( ddirs=ddir , odir=odir ,  heightLimits.km=heightLimits.km , timeRes.s=timeRes.s , beginTime=beginTime , endTime=endTime , fitFun=leastSquare.lvmrq , absLimit=absLimit , diffLimit=diffLimit , maxLambda=maxLambda , maxIter=maxIter , plotTest=FALSE , plotFit=plotFit , absCalib=FALSE , TiIsotropic=TRUE , TeIsotropic=TRUE , recursive=recursive , aprioriFunction=ISaprioriH , scaleFun=acfscales , siteScales=NULL, calScale=calScale , nCores=nCores , hVi=80 , hTeTi=Inf , hTi=Inf )



}
