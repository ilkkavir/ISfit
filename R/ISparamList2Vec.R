ISparamList2Vec <- function(ele,ion,cSite){
#
# Converts the plasma parameters given as vector "ele", list "ion", and vector siteScales, into a parameter vector
# used in ISparamfit
#
# INPUT:
#  ele        c(Ne,Te_par,Te_perp,nu_en,ve_x,ve_y,ve_z)
#  ion        list( c(m1,N1,T1_par,T1_perp,nu_1n,v1_x,v1_y,v1_z) , c(m2,N2,T2_par,T2_perp,nu_2n,v2_x,v2_y,v2_z) , ... )
#  cSite       a vector of ACF scaling factors for each site
#
# OUTPUT:
#  paramVec    a parameter vector used in ISparamfit
#
# I. Virtanen 2012
#

  parvec <- c( ele[1] , ion[[1]][3:4] , ele[2:3] , ion[[1]][5], ion[[1]][6:8] , sapply(ion,FUN=function(x){x[2]}) , cSite )
  names(parvec) <- c('Ne','Tipar','Tiperp','Tepar','Teperp','coll','Vx','Vy','Vz',paste('ion',seq(length(ion)),sep=''),paste('site',seq(length(cSite))))
  return(parvec)
  
} 
