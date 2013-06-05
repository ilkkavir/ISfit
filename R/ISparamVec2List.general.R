ISparamVec2List.general <- function(param,nIon){
#
# Conversion from parameter vector used in ISparamfit.general into lists of electron and ion parameters
# 
# INPUT:
#  param        a vector of plasma parameters and ACF scales:
#                 c( Ne , Te_par , Te_perp , nu_en , ve_x , ve_y , ve_z ,
#               m1 , N1 , T1_par , T1_perp , nu_1n , v1_x , v1_y , v1_z ,
#               m2 , N2 , T2_par , ...
#               ...
#               mnIon , ...                                   , vnIon_z ,
#               s_1 , s_2 , ... , s_Nsite)
#  nIon         Number of ion species included in the parameter vector
#
#
  # electron parameters
  ele      <- param[1:7]

  # ion parameters
  ion      <- vector(mode='list',length=nIon)
  i        <- 8
  for(k in seq(nIon)){
    ion[[k]] <- param[i:(i+7)]
    i        <- i + 8
  }

  # ACF scales
  cSite    <- param[i:length(param)]


  return(list(ele=ele,ion=ion,cSite=cSite))
  
} # ISparamVec2List.general
