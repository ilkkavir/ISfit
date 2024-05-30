ISfit.3D.BAFIM.multinode <- function( ddirs='.' , odir='.' ,  heightLimits.km=NA , timeRes.s=60 , timeResFirst.s=timeRes.s , mlatLimits.deg=c(-90,90),mlonLimits.deg=c(-360,360), beginTime=c(1970,1,1,0,0,0) , endTime=c(2100,1,1,0,0,0) , fitFun=leastSquare.lvmrq , absLimit=5 , diffLimit=1e-2 , maxLambda=1e30 , maxIter=10 , absCalib=FALSE , TiIsotropic=TRUE , TeIsotropic=TRUE , recursive=TRUE , aprioriFunction=ISaprioriH , scaleFun=acfscales , siteScales=NULL, calScale=1, MCMCsettings=list( niter=10000 , updatecov=100 , burninlength=5000 , outputlength=5000 ) , maxdev=2 , trueHessian=FALSE , nCores=1 , reverseTime=FALSE , burnin.s=1800 , cl=NULL , ... ){




    ## find the data directories for each mlat, mlon bin
    iddirs <- ddirSelect(mlatLimits.deg,mlonLimits.deg,heightLimits.km,ddirs)

    nbinmlat <- length(mlatLimits.deg) - 1
    nbinmlon <- length(mlonLimits.deg) - 1

    nbinlatlon <- nbinmlat * nbinmlon

    if(is.null(cl)){
        for(ibin in seq(nbinlatlon)){
            ISfit:::ISfit.3D.BAFIM(
                        ddirs=ddirs,
                        odir=odir,
                        heightLimits.km=heightLimits.km,
                        timeRes.s=timeRes.s,
                        timeResFirst.s=timeResFirst.s,
                        mlatLimits.deg=mlatLimits.deg,
                        mlonLimits.deg=mlonLimits.deg,
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
                        scaleFun=scaleFun,
                        siteScales=siteScales,
                        calScale=calScale,
                        MCMCsettings=MCMCsettings,
                        maxdev=maxdev,
                        trueHessian=trueHessian,
                        nCores=nCores,
                        reverseTime=FALSE,
                        burnin.s=burnin.s,
                        iddirs = iddirs,
                        imlatlon=ibin,
                        ...
                    )
        }
    }else{
        # this has not been tested yet!!
        snow::clusterApplyLB( cl ,
                           seq(nbinlatlon) ,
                           fun=ISfit.3D.BAFIM ,
                           ddirs=ddirs,
                           odir=odir,
                           heightLimits.km=heightLimits.km,
                           timeRes.s=timeRes.s,
                           timeResFirst.s=timeResFirst.s,
                           mlatLimits.deg=mlatLimits.deg,
                           mlonLimits.deg=mlonLimits.deg,
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
                           scaleFun=scaleFun,
                           siteScales=siteScales,
                           calScale=calScale,
                           MCMCsettings=MCMCsettings,
                           maxdev=maxdev,
                           trueHessian=trueHessian,
                           nCores=nCores,
                           reverseTime=FALSE,
                           burnin.s=burnin.s,
                           iddirs = iddirs,
                           ...
                           )
    }
}
