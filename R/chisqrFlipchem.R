chisqrFlipchem <- function( measData , measVar , dirtheData , param , aprioriTheory , aprioriMeas , invAprioriCovar , flipchem , flipchemStd , lat , lon, h , scaleFun ,  mIon , nIon ,  ... ){
#
# sum of residuals 
#
#
# INPUT:
#  measData         a vector of measured data values
#  measVar          variances of the measurements
#  dirtheData       data values from direct theory
#  param            current parameter values
#  aprioriTheory    theory matrix of linear apriori information
#  aprioriMeas      imaginary "apriori measurements"
#  invAprioriCovar  inverse of apriori measurement covariance matrix (i.e. apriori Fisher information matrix)
#  flipchem         a flipchem function handle
#  flipchemStd      standard deviation of the flipchem composition prediction
#  lat              geodetic latitude [deg] 
#  lon              geodetic longitude [deg]
#  h                geodetic altitude [km]
#  scaleFun         parameter scaling function
#  mIon             ion masses
#  nIon             number of ion species
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

    # identify the molecular and O+ ions in the mIon arrray
    MolIon <- which(mIon >= 28 & mIon <= 32)
    OpIon <- which(mIon > 15.8 & mIon < 16.2)

    # extra cost for deviations from flipchem
    fcCost <- 0
    if(length(MolIon)>0){
        if(MolIon<=nIon){
            fcCost <- fcCost + (sparam[8+MolIon] - fcMol)**2 / flipchemStd**2
        }
    }
    if(length(OpIon)>0){
        if(OpIon<=nIon){
            fcCost <- fcCost + (sparam[8+OpIon] - fcOp)**2 / flipchemStd**2
        }
    }

    # cost function from the direct theory
    dtCost <- sum( abs(measData - dirtheData)**2 / measVar )

    # cost function from the linear prior
    prCost <- t( aprioriTheory%*%param - aprioriMeas ) %*% invAprioriCovar %*% ( aprioriTheory%*%param - aprioriMeas )

    
    # return the sum of the cost functions
    return( fcCost + dtCost + prCost )
    
}
