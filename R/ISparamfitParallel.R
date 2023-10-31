ISparamfitParallel <- function(h,acf,var,lags,nData,fSite,aSite,kSite,iSite,B,apriori,directTheory,absLimit,diffLimit,scaleFun,maxLambda,maxIter,fitFun,MCMCsettings,trueHessian,heights,latitude,longitude,fitGate,...){



    fitpar <- NA

    if(fitGate[h]){
        fitpar <- ISparamfit(
            acf             = unlist(acf[[h]]),
            var             = unlist(var[[h]]),
            lags            = unlist(lags[[h]]),
            nData           = sum(nData[[h]]),
            fSite           = fSite,
            aSite           = aSite[[h]],
            kSite           = kSite[[h]],
            iSite           = unlist(iSite[[h]]),
            B               = B[h,],
            initParam       = apriori[[h]]$aprioriParam,
            invAprioriCovar = apriori[[h]]$invAprioriCovar,
            aprioriTheory   = apriori[[h]]$aprioriTheory,
            aprioriMeas     = apriori[[h]]$aprioriMeas,
            mIon            = apriori[[h]]$mIon,
            nIon            = apriori[[h]]$nIon,
            paramLimits     = apriori[[h]]$limitParam,
            directTheory    = directTheory,
            absLimit        = absLimit,
            diffLimit       = diffLimit,
            scaleFun        = scaleFun,
            scale           = apriori[[h]]$parScales,
            plotTest        = FALSE,
            plotFit         = FALSE,
            maxLambda       = maxLambda,
            maxIter         = maxIter,
            fitFun          = fitFun,
            MCMCsettings    = MCMCsettings,
            trueHessian     = trueHessian,
            h               = heights[h],
            lat             = latitude[h],
            lon             = longitude[h],
            flipchem        = apriori[[h]]$flipchem,
            flipchemStd     = apriori[[h]]$flipchemStd,
            aprioriUpdateFunction = apriori[[h]]$aprioriUpdateFunction
        )
    }

    return(fitpar)




    }
