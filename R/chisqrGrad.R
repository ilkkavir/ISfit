chisqrGrad <- function( measData , measVar , param , directTheory , aprioriTheory , aprioriMeas , invAprioriCovar ,  ... )
    {

        #
        # 
        # Gradient of the chi-square with finite differences
        # 
        #
        # I. Virtanen 2014, 2023
        #
        #

        diff <- .0001
        
        npar <- length(param)

        chisqrGrad <- rep(0,npar)
        
        # chi-squared at the current point
        chisqr0 <-  chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=param , ... ) ,
                              param=param , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )

        # The gradients
        for( k in seq( npar ) ){
            # chi-squared at points param[k] + diff
            par2 <- param
            par2[k] <- par2[k] + diff
            chisqrGrad[k] <- chisqrGrad[k] + chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                         param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
            # subtract chi-squared at points param[k] - diff
            par2 <- param
            par2[k] <- par2[k] - diff
            chisqrGrad[k] <- chisqrGrad[k] - chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                         param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
            # divide by 2'diff
            chisqrGrad[k] <- chisqrGrad[k] / diff / 2
        }

        
        ## # The gradients
        ## for( k in seq( npar ) ){
        ##     # chi-squared at points param[k] + diff
        ##     par2 <- param
        ##     par2[k] <- par2[k] + diff
        ##     chisqrGrad[k] <- ( chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) , param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... ) - chisqr0 ) / diff
        ## }

        return(chisqrGrad)
        
    }
