leastSquare.lvmrq <- function(measData , measVar , initParam , aprioriTheory, aprioriMeas , invAprioriCovar , paramLimits , directTheory , diffLimit=1e-2 , absLimit=5 , maxIter=50 , maxLambda=1e20 , plotTest=F , plotFit=F , trueHessian=F , aprioriUpdateFunction=NULL , ... ){
# 
# Non-linear least-squares fit using the Levenberg-Marquardt algorithm
# 
# 
# INPUT:
#  measData       A vector of measured data values
#  measVar        A vector of measurement variances
#  initParam      c(p1,p2,p3, ... ,pN ) a vector of initial parameter values
#  aprioriTheory  M x N theory matrix for linear apriori theory, M>0. Use a row of zeros to disable apriori information.
#  aprioriMeas    M vector of the imaginary "measurements" correspoding the rows of aprioriTheory
#  invAprioriCovar   M x M inverse of the covariance matrix of aprioriMeas
#  paramLimits    2 x N matrix of smallest (paramLimits[1,]) and largest (paramLimits[2,]) allowed value for each parameter 
#  directTheory   function( param , ... ) , where param is a vector similar to initParam and ... are the additional arguments
#                 given to leastSquare.lvmrq. Should return a vector similar to measData
#  diffLimit      If relative change is chi-squared is smaller than diffLimit AND absolute residual is smaller than absLimit,
#                 the iteration is stopped
#  absLimit       see above
#  maxIter        Maximum number of iteration steps
#  plotTest       If true, a plot of the measured data and the currently tested direct theory is generated at each iteration
#  plotFit        If true, a plot of the measured data and the fitted direct theory is generated
#  trueHessian    If true, the Hessian matrix for final error estimation is calculated from finite differences of chi-squared,
#                 otherwise the real part of complex Hessian as calculated from gradient of direct theory is used. With
#                 trueHessian=TRUE the analysis is somewhat more time-consuming.
#  ...            Additional parameters passed to directTheory
# 
# OUTPUT:
#  param          The solved least squares point
#  covar          Covariance matrix approximated from linearised theory in vicinity of the least squares solution
#  chisqr         Chi-squared of the solution
#  nIter          Number of iterations in the algorithm
#  fitStatus      Status of the fit: 0=succeeded, 1=exceeded maxIter, 2=exceeded maxLambda
# 
# I. Virtanen 2012
# 

  # length of the measurement vector
  nMeas            <- length(measData)

  # Parameter vector, the apriori values are used as a starting point
  # Current best parameters are in param[1,], currently tested values will be stored in param[2,]
  param            <- matrix( rep(initParam,2) , nrow=2 , byrow=T )

  # measurement values from direct theory
  # Current best values (corresponding param[1,]) are stored in dirtheData[1,], the currently tested values
  # (corresponding param[2,]) are stored in dirtheData[2,]
  dirtheData       <- matrix( rep(directTheory( param[1,] , ... ),2) , nrow=2 , byrow=T )

  # chi-square
  chisqr           <- rep( chiSquare( measData , measVar , dirtheData[1,] , param[1,] ,
                                     aprioriTheory, aprioriMeas , invAprioriCovar , ... ) , 2 )

  # status of the fit, initialise to 1
  fitStatus        <- 1
  
  # Initialise the iteration counter
  nIter            <- 0

  # The gradient and Hessian will need to be updated in the first iteration
  newPoint         <- T

  # start value for lambda
  lambda           <- 1e-3

  # repeat until convergence or until reaching maximum number of iterations
  repeat{

    # check the iteration counter, and break if the maximum value is reached
    if( nIter > maxIter ) break

    # calculate gradient and Hessian if necessary
    if( newPoint ){

        # Increment the iteration counter
        nIter          <- nIter + 1
        
        ## if(!is.null(aprioriUpdateFunction)){
        ##     aprioriNew <- aprioriUpdateFunction(param=param[1,],aprioriTheory=aprioriTheory,aprioriMeas=aprioriMeas,invAprioriCovar=invAprioriCovar,paramLimits=paramLimits , ... )
        ##     aprioriTheory <- aprioriNew$A
        ##     aprioriMeas <- aprioriNew$m
        ##     invAprioriCovar <- aprioriNew$invAprioriCovar
        ## }

        if( trueHessian ){

            alpha <- .5 * hessian( measData , measVar , param[1,] , directTheory , aprioriTheory , aprioriMeas , invAprioriCovar ,  ... )
            beta <- -.5 * chisqrGrad( measData , measVar , param[1,] , directTheory , aprioriTheory , aprioriMeas , invAprioriCovar , ... )


        }else{
          
            # calculate  gradients
            dirtheGrad      <- dtGradient( dirtheData[1,] , directTheory , param[1,] , ... )

            # Hessian
            alpha          <- dtAlpha( dirtheGrad , measVar , aprioriTheory , invAprioriCovar )

            # Gradient of the cost function
            beta           <- dtBeta( measData , measVar , dirtheGrad , dirtheData[1,] , param[1,] ,
                                     aprioriTheory , aprioriMeas , invAprioriCovar )
            
        }

        
    }

    # add lambda to the diagonal of alpha
    alpha2           <- alpha
    diag(alpha2)     <- diag(alpha2) + lambda

    # solve the problem alpha2 * dx = beta
    dx               <- tryCatch(solve( alpha2 , beta ),error=function(e){return(NA)})
    if(any(is.na(dx))){
      # if the problem could not be solved, break the iteration with exit code 3
      fitStatus <- 3
      break
    }

    # add the solved dx to the parameters
    param[2,]        <- param[1,] + dx

    # force the new point inside the given limits
    param[2,]        <- pmax(param[2,],paramLimits[1,])
    param[2,]        <- pmin(param[2,],paramLimits[2,])

    # direct theory with the new parameters
    dirtheData[2,]   <- directTheory( param[2,] , ... )

    # chi-squared with the current direct theory
    chisqr[2]        <- chiSquare( measData , measVar , dirtheData[2,] , param[2,] , aprioriTheory , aprioriMeas , invAprioriCovar , ... )

    # plot the current best direct theory
    if(plotTest){
      plot(Re(measData),main=chisqr[1]/nMeas,ylim=c(-1,1)*1.5*max(abs(measData)))
      points(Im(measData),col='red')
      lines(Re(dirtheData[1,]))
      lines(Im(dirtheData[1,]),col='red')
    }

    # Test if chi-squared improved
    if(chisqr[2] < chisqr[1]){

      # if yes, decrease lambda and continue to next iteration
      lambda         <- 0.1*lambda
      newPoint       <- T
      param[1,]      <- param[2,]
      dirtheData[1,] <- dirtheData[2,]
      chiChange      <- (chisqr[1]-chisqr[2]) / chisqr[2]
      chisqr[1]      <- chisqr[2]


        
      if((chisqr[1]/nMeas) < absLimit){
        # if relative change in chi-squared is below the limit, break the loop
        if(chiChange < diffLimit){
          fitStatus    <- 0
          break
        }
      }
    }else{

      # if no, try a larger lambda
      lambda         <- 10*lambda
      newPoint       <- F

      # if lambda is above limit, break with status 2
      if(lambda > maxLambda){
        fitStatus    <- 2
        break
      }
    }
    
  } # end of iteration



  if( trueHessian ){
      
      # the true Hessian for error estimation
      H <- hessian( measData , measVar , param[1,] , directTheory , aprioriTheory , aprioriMeas , invAprioriCovar ,  ... )

  }else{
      
      # real part of the complex Hessian as calculated from direct theory gradients
      dirtheGrad <- dtGradient( dirtheData[1,] , directTheory , param[1,] , ... )
      H <- 2 * dtAlpha( dirtheGrad , measVar , aprioriTheory , invAprioriCovar )

  }

  # the covariance matrix
  covar <- tryCatch( solve(.5*H) , error=function(e){ H[,]*NA })
      
  # Do not accept negative variances
  if(any(diag(covar)<0 | is.na(diag(covar)))){
      covar <- covar[,]*NA
      param <- param[,]*NA
      fitStatus <- 4
      chisqr[] <- NA
  }

  # plot the current best direct theory
  if( ( plotTest | plotFit )  ){
      plot(Re(measData),main=sprintf("Residual %6.2f, Fit status %i ",chisqr[1]/nMeas,fitStatus),ylim=c(-1,1)*1.5*max(abs(measData)))
      points(Im(measData),col='red')
      if( fitStatus < 3 ){
          lines(Re(dirtheData[1,]))
          lines(Im(dirtheData[1,]),col='red')
      }
  }
  
  # return the solution and some other information
  return( list( param=param[1,] , covar=covar , chisqr=chisqr[1]/nMeas , nIter=nIter , fitStatus=fitStatus ) )
  
}

