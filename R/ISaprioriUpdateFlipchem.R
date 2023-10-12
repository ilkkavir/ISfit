ISaprioriUpdateFlipchem <- function( param , paramLimits , aprioriMeas , flipchem ,  lat , lon, h , scaleFun ,   ... ){
#
# sum of residuals 
#
#
# INPUT:
#  param            current parameter values (scaled)    
#  aprioriMeas      imaginary "apriori measurements"
#  flipchem         a flipchem function handle
#  flipchemStd      standard deviation of the flipchem composition prediction
#  lat              geodetic latitude [deg] 
#  lon              geodetic longitude [deg]
#  h                geodetic altitude [km]
#  scaleFun         parameter scaling function
#
#    
# OUTPUT:
#  chi-squared with the current parameters
#  
# I. Virtanen 2012
#

    
    # scale parameters
    sparam         <- scaleFun( param , ... , inverse=T)

    # ion and electron temperatures ( (Tpar+2*Tperp)/3 )
    Te <- (sparam[4] + 2*sparam[5]) / 3
    Ti <- (sparam[2] + 2*sparam[3]) / 3
    
    # call flipchem with parallel temperatures to be sure that 
    fcout <- flipchem$get_point(lat,lon,h,sparam[1],Te,Ti)

    # O+ and molecular ion fractions from flicphem
    fcOp <- fcout[[4]] / sparam[1]
    fcMol <- (fcout[[5]]+fcout[[6]]+fcout[[7]]) / sparam[1]


    nprMeas  <- length(aprioriMeas)
    aprioriMeas[nprMeas-1] <- fcMol
    aprioriMeas[nprMeas] <- fcOp
    
    # return the sum of the cost functions
    return( aprioriMeas )
    
}
