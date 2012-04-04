ISaprioricovar.default.1D <- function( nIon , nSite , ... ){

  #                Ne   dTe  dnu_en  dve    m1     T1   nu_1n   vi            mi    Ni   dTi   dnu_in  dvi                   siteScale
  return(diag( c( 1.0 , 1e4 , 1e-3 , 1e-3 , 1e-3 , 1e4 , 1e-3 , 1e4 , rep( c(1e-3 , 1 , 1e-3 , 1e-3 , 1e-3) , (nIon-1) ) , rep(1e-3,nSite) ) ) )
 
} # ISaprioricovar.default.1D
