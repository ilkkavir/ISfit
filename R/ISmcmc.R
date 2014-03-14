ISmcmc <- function( measData , measVar , initParam , aprioriTheory , aprioriMeas , invAprioriCovar , paramLimits , directTheory , MCMCsettings=list() , ... )
    {
        initres <- leastSquare.lvmrq(measData=measData , measVar=measVar , initParam=initParam , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , paramLimits=paramLimits , directTheory=directTheory , ... )
        
        # use the iteration result as strating point if the iteration was successful
        # otherwise use the prior model but make sure that the parameter variances
        # are not too large, otherwise even the adaptive MCMC will be in trouble
        if( initres[["fitStatus"]]){
            init2 <- initParam
            initcovar <- diag(pmin(.1,diag(solve(t(aprioriTheory) %*% invAprioriCovar %*% aprioriTheory + diag(rep(1e-16,length(init2)))))))
        }else{
            init2 <- initres[["param"]]
            initcovar <- initres[["covar"]]
        }
        
        names(init2) <- names(initParam)
        
        
        mcmcargs <- c( list( f=SS , p=init2 , measData=measData , measVar=measVar , directTheory=directTheory , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , lower=paramLimits[1,] , upper=paramLimits[2,] , jump=initcovar*2.4^2/length(init2) ) , MCMCsettings , list( ... ))
        
        resmcmc <- do.call( modMCMC , mcmcargs  )
        
        # return the iteration output list padded with the MCMC output
        initres$MCMC <- resmcmc
        return(initres)
        
    }

