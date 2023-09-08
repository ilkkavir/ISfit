ISparamfitParallel <- function(h,acf,var,lags,nData,fSite,aSite,kSite,iSite,B,initParam,apriori,mIon,nIon,paramLimits,directTheory,absLimit,diffLimit,scaleFun,scale,maxLambda,maxIter,fitFun,MCMCsettings,trueHessian,heights,fitGate){



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
            initParam       = initParam[[h]],
            invAprioriCovar = apriori[[h]]$invAprioriCovar,
            aprioriTheory   = apriori[[h]]$aprioriTheory,
            aprioriMeas     = apriori[[h]]$aprioriMeas,
            mIon            = mIon,
            nIon            = nIon,
            paramLimits     = paramLimits[[h]],
            directTheory    = directTheory,
            absLimit        = absLimit,
            diffLimit       = diffLimit,
            scaleFun        = scaleFun,
            scale           = scale[[h]],
            plotTest        = FALSE,
            plotFit         = FALSE,
            maxLambda       = maxLambda,
            maxIter         = maxIter,
            fitFun          = fitFun,
            MCMCsettings    = MCMCsettings,
            trueHessian     = trueHessian,
            h               = heights[h]
        )
    }
    
    

    return(fitpar)




    }
