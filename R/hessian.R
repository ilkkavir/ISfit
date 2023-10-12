hessian <- function( measData , measVar , param , directTheory , aprioriTheory , aprioriMeas , invAprioriCovar , ... )
    {

        #
        # 
        # Hessian matrix from actual derivatives of chi-squared by means of finite differences
        # 
        #
        # I. Virtanen 2014, 2023
        #
        #

        diff <- .0001
        
        npar <- length(param)

        H <- matrix(0,ncol=npar,nrow=npar)
        

        # lower triangular part without the main diagonal
        for( k in seq( npar - 1 ) ){
            for( l in seq( (k+1) , npar ) ){
                # chi-squared at points param[c(k,l)]+diff
                par2 <- param
                par2[c(k,l)] <- par2[c(k,l)] + diff
                H[k,l] <-  H[k,l] + chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                              param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
                # chi-squared at points param[c(k,l)]-diff
                par2 <- param
                par2[c(k,l)] <- par2[c(k,l)] - diff
                H[k,l] <- H[k,l] + chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                             param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
                # chi-squared at points param[c(k,l)]+c(diff,-diff)
                par2 <- param
                par2[c(k,l)] <- par2[c(k,l)] + c( diff , -diff )
                H[k,l] <- H[k,l] - chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                             param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
                # chi-squared at points param[c(k,l)]+c(-diff,diff)
                par2 <- param
                par2[c(k,l)] <- par2[c(k,l)] + c( -diff , diff )
                H[k,l] <- H[k,l] - chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                             param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
                # normalization
                H[k,l] <- H[k,l] / ( 4 * diff**2 )
            }
        }

        # Fill the upper triangular part
        H <- H + t(H)

        # chi-squared at the current point
        chisqr0 <-  chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=param , ... ) ,
                              param=param , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )

        # The main diagonal
        for( k in seq( npar ) ){
            # chi-squared at points param[k] + diff
            par2 <- param
            par2[k] <- par2[k] + diff
            H[k,k] <- H[k,k] + chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                         param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
            # chi-squared at points param[k] - diff
            par2 <- param
            par2[k] <- par2[k] - diff
            H[k,k] <- H[k,k] + chiSquare( measData=measData , measVar=measVar , dirtheData=directTheory( param=par2 , ... ) ,
                                         param=par2 , aprioriTheory=aprioriTheory , aprioriMeas=aprioriMeas , invAprioriCovar=invAprioriCovar , ... )
            # subtract twice the chi-squared at the central point
            H[k,k] <- H[k,k] - 2*chisqr0

            # normalization
            H[k,k] <- H[k,k] / diff**2
        }

        return(H)
        
    }
