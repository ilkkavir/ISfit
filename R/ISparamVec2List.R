ISparamVec2List <- function(param,mIon){
#
# Conversion from parameter vector used in ISparamfit into lists of electron and ion parameters
# 
# INPUT:
#  param        a vector of plasma parameters and ACF scales:
#                 c( Ne , Ti_par , Ti_perp , Te_par/Ti_par , Te_perp/Ti_perp,
#                    nu_in , vi_x , vi_y , vi_z , p_m2 , p_m3 , ... p_nIon ,
#                    s_1 , s_2 , ... , s_Nsite)
#  mIon         Ion masses
#
#
  # electron parameters
  ele      <- c( param[1] , param[4:5] , param[6]*0.35714 , param[7:9])

  nIon <- length(mIon)
  
  # ion parameters
  ion      <- vector(mode='list',length=nIon)
  ion[[1]] <- c( mIon[1] , param[10] , param[2:3]  , param[6:9])
  if( nIon > 1){
      for(k in seq(2,nIon)){
          ion[[k]] <- c( mIon[k] , param[9+k] , param[2:3] , param[6:9] )
      }
  }

  # ACF scales
  cSite    <- param[(10+nIon):length(param)]


  return(list(ele=ele,ion=ion,cSite=cSite))
  
}
