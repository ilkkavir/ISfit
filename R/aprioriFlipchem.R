aprioriFlipchem <- function( param , flipchem , flipchemStd ,  lat , lon, h , scaleFun ,   ... ){
#
# Apriori theory from linear apporiximation of flipchem composition
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

    dp <- 1e-4
    
    # scale parameters
    sparam         <- scaleFun( param , ... , inverse=T)

    # ion and electron temperatures ( (Tpar+2*Tperp)/3 )
    Te <- (sparam[4] + 2*sparam[5]) / 3
    Ti <- (sparam[2] + 2*sparam[3]) / 3
    
    # call flipchem 
    fcout <- flipchem$get_point(lat,lon,h,sparam[1],Te,Ti)

    # O+ ion fraction from flicphem
    fcOp <- fcout[[4]] / sparam[1]

    # derivatives with respect to Ne, Ti, Te, Coll
    dOp  <-  rep(0,6)
    for( k in seq(6)){
        param2 <- param
        dparam <- rep(0,6)
        dparam[k] <- dparam[k] + dp
        param2[1:6] <- param2[1:6] + dparam
        sparam2 <- scaleFun( param2 , ... , inverse=T)
        Te <- (sparam2[4] + 2*sparam2[5]) / 3
        Ti <- (sparam2[2] + 2*sparam2[3]) / 3
        fcout <- flipchem$get_point(lat,lon,h,sparam2[1],Te,Ti)
        dfcOp <- fcout[[4]] / sparam[1] - fcOp
        dOp[k] <- dfcOp/dp
    }

    apriori <- list()
    apriori$A <- param*0
    apriori$A[1:6] <- -dOp
    apriori$A[11] <- 1
    apriori$m <- fcOp - sum(dOp*param[1:6])
    apriori$var <- flipchemStd**2
#    print(c(sum(dOp*param[1:6])+apriori$m,fcOp,apriori$m))

    
    
    # return the flipchem prior
    return( apriori )
    
}
