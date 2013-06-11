ISmcmc <- function( measData , measVar , initParam , aprioriTheory , aprioriMeas , invAprioriCovar , paramLimits , directTheory , ... )
    {
        initres <- leastSquare.lvmrq(measData=measData , measVar=measVar , initParam=initParam , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , paramLimits=paramLimits , directTheory=directTheory , ... )

        init2 <- initres$param
        names(init2) <- names(initParam)
print('mcmc')
        res <- modMCMC( f=SS , p=init2 , measData=measData , measVar=measVar , directTheory=directTheory , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , lower=paramLimits[1,] , upper=paramLimits[2,] , jump=initres$covar ,  ... )

        return(list(iter=initres,mcmc=res))
    }
