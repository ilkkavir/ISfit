ISmcmc <- function( measData , measVar , initParam , aprioriTheory , aprioriMeas , invAprioriCovar , paramLimits , directTheory , MCMCsettings=list() , ... )
    {
        initres <- leastSquare.lvmrq(measData=measData , measVar=measVar , initParam=initParam , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , paramLimits=paramLimits , directTheory=directTheory , ... )

        init2 <- initres$param
        names(init2) <- names(initParam)

        mcmcargs <- c( list( f=SS , p=init2 , measData=measData , measVar=measVar , directTheory=directTheory , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , lower=paramLimits[1,] , upper=paramLimits[2,] , jump=initres$covar*2.4**2/length(init2) ) , MCMCsettings , list( ... ))

        resmcmc <- do.call( modMCMC , mcmcargs  )

#        npar <- length(initParam)
#        res <- list(
#            param=resmcmc[["bestpar"]],
#            covar=matrix(0,ncol=npar,nrow=npar),
#            chisqr=resmcmc[["bestfunp"]]/length(measData),
#            nIter=NULL,
#            fitStatus=0
#            )
#        for(k in seq(npar)){
#            for(l in seq(k,npar)){
#                res[["covar"]][k,l] <- mean( (resmcmc[["pars"]][,k] - resmcmc[["bestpar"]][k] ) * (resmcmc[["pars"]][,l] - resmcmc[["bestpar"]][l] ))
#            }
#        }
#
#        res[["covar"]] <- res[["covar"]] + t(res[["covar"]])
#        diag(res[["covar"]]) <- diag(res[["covar"]])/2
#
#        cat(sprintf("%.2f",init2),'\n')
#        cat(sprintf("%.2f",res$param),'\n')
#
        return(list(param=initres$param,covar=initres$covar,chisqr=initres$chisqr,fitStatus=initres$fitStatus,MCMC=resmcmc))

    }


