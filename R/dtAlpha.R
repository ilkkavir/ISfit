dtAlpha <- function( dirtheGrad , measVar , aprioriTheory , invAprioriCovar)
    {
        #
        # Real part of one half of complex Hessian matrix at the current point.
        # A sometimes better approximation of the Hessian, by means of finite
        # differences of chi-squared estimates, can be calculated with the
        # function "hessian".
        #
        # In incoherent scatter plasma parameter estimation the Hessian *should*
        # be almost real because the plasma parameters are real.
        #  
        # INPUT:
        #  dirtheGrad      a matrix of numerical differences from dtGradient
        #  measVar         variances of the measured data points
        #  aprioriTheory   theory matrix of the linear apriori theory
        #  invAprioriCovar inverse of apriori covariance matrix
        #  

#  return( Re( dirtheGrad %*% diag(1/measVar) %*% t(Conj(dirtheGrad))) + t(aprioriTheory)%*%invAprioriCovar%*%aprioriTheory )

        return( Re( dirtheGrad %*% (t(Conj(dirtheGrad))/measVar) ) + t(aprioriTheory)%*%invAprioriCovar%*%aprioriTheory )
        
    }
