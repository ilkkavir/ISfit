ISfit.3D.BAFIM <- function( ddirs='.' , odir='.' ,  heightLimits.km=NA , timeRes.s=60 , beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , fitFun=leastSquare.lvmrq , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=10 , absCalib=FALSE , TiIsotropic=TRUE , TeIsotropic=TRUE , recursive=TRUE , aprioriFunction=ISaprioriH , scaleFun=acfscales , siteScales=NULL, calScale=1, MCMCsettings=list( niter=10000 , updatecov=100 , burninlength=5000 , outputlength=5000 ) , maxdev=2 , trueHessian=FALSE , nCores=1 , reverseTime=FALSE , burnin.s=1800 , ... )
  {

      # 3D incoherent scatter plasma parameter fit with Baeysian filtering and smoothing using LPI output files in ddirs
      #


      odirFw <- paste(odir,'Forward',sep='')
      odirRev <- paste(odir,'Reverse',sep='')
      odirSmooth <- paste(odir,'smoother',sep='')

      # forward filter
      ISfit.3D(
          ddirs=ddirs,
          odir=odirFw,
          heightLimits.km=heightLimits.km,
          timeRes.s=timeRes.s,
          beginTime=beginTime,
          endTime=endTime,
          fitFun=fitFun,
          absLimit=absLimit,
          diffLimit=diffLimit,
          maxLambda=maxLambda,
          maxIter=maxIter,
          absCalib=absCalib,
          TiIsotropic=TiIsotropic,
          TeIsotropic=TeIsotropic,
          recursive=recursive,
          aprioriFunction=ISaprioriBAFIM, ###
          scaleFun=scaleFun,
          siteScales=siteScales,
          calScale=calScale,
          MCMCsettings=MCMCsettings,
          maxdev=maxdev,
          trueHessian=trueHessian,
          nCores=nCores,
          reverseTime=FALSE,
          burnin.s=burnin.s,
          ...
      )


      # reverse in time filter
      ISfit.3D(
          ddirs=ddirs,
          odir=odirRev,
          heightLimits.km=heightLimits.km,
          timeRes.s=timeRes.s,
          beginTime=beginTime,
          endTime=endTime,
          fitFun=fitFun,
          absLimit=absLimit,
          diffLimit=diffLimit,
          maxLambda=maxLambda,
          maxIter=maxIter,
          absCalib=absCalib,
          TiIsotropic=TiIsotropic,
          TeIsotropic=TeIsotropic,
          recursive=recursive,
          aprioriFunction=ISaprioriBAFIM, ###
          scaleFun=scaleFun,
          siteScales=siteScales,
          calScale=calScale,
          MCMCsettings=MCMCsettings,
          maxdev=maxdev,
          trueHessian=trueHessian,
          nCores=nCores,
          reverseTime=TRUE,
          burnin.s=burnin.s,
          ...
      )

      # smoother
      ISfit.3D(
          ddirs=ddirs,
          odir=odirSmooth,
          resdirForward=odirFw,
          resdirReverse=odirRev,
          heightLimits.km=heightLimits.km,
          timeRes.s=timeRes.s,
          beginTime=beginTime,
          endTime=endTime,
          fitFun=fitFun,
          absLimit=absLimit,
          diffLimit=diffLimit,
          maxLambda=maxLambda,
          maxIter=maxIter,
          absCalib=absCalib,
          TiIsotropic=TiIsotropic,
          TeIsotropic=TeIsotropic,
          recursive=recursive,
          aprioriFunction=ISaprioriBAFIMsmoother, ###
          scaleFun=scaleFun,
          siteScales=siteScales,
          calScale=calScale,
          MCMCsettings=MCMCsettings,
          maxdev=maxdev,
          trueHessian=trueHessian,
          nCores=nCores,
          reverseTime=FALSE,
          burnin.s=0,
          ...
      )


      

  }
