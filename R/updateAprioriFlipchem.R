updateAprioriFlipchem  <- function( aprioriTheory , aprioriMeas , invAprioriCovar , param , flipchem , flipchemStd ,  lat , lon, h , scaleFun ,   ... ){
#
# update the prior theory row that contains flipchem parameters
#
#
# INPUT:
#  param            current parameter values (scaled)    
#  flipchem         a flipchem function handle
#  flipchemStd      standard deviation of the flipchem composition prediction
#  lat              geodetic latitude [deg] 
#  lon              geodetic longitude [deg]
#  h                geodetic altitude [km]
#  scaleFun         parameter scaling function
#
#    
# OUTPUT:
#  a list with entries:
#      A           the theory matrix row
#      m           the prior measurement
#      var         variance of the prior measurement
#  
# I. Virtanen 2023
#


    fcPrior <- aprioriFlipchem( param=param , flipchem=flipchem , flipchemStd=flipchemStd ,  lat=lat , lon=lon, h=h , scaleFun=scaleFun ,   ... )

    nMeas <- length(aprioriMeas)
    aprioriTheory[nMeas,] <- fcPrior$A
    aprioriMeas[nMeas] <- fcPrior$m
    priorCovar <- solve(invAprioriCovar)
    priorCovar[nMeas,nMeas] <- fcPrior$var
    
    
    return(list(A=aprioriTheory,m=aprioriMeas,invAprioriCovar=solve(priorCovar)))
    
}
