iriParamsParFun <- function(h,date,latitude,longitude,height,fitGate,okData){

    IRIpar <- NA
    if (fitGate[h]&okData[h]){
        IRIpar  <- iriParams( time=date ,latitude=latitude[h],longitude=longitude[h],heights=height[h])
    }

    return(IRIpar)

}
