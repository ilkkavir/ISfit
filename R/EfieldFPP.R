EfieldFPP <- function(PP,hmin=200,hmax=400,vipar0=FALSE,...){
    # F-region electric field with ISfit output list
    # as input

    hh <- which( (PP[["height"]]>=hmin) & (PP[["height"]] <= 400) )

    vi <- PP[["param"]][hh,c('Vix','Viy','Viz')]
    vicov <- lapply(PP[["covar"]][hh],FUN=function(x){return(x[c('Vix','Viy','Viz'),c('Vix','Viy','Viz')])})
    B <- PP[["B"]][hh,]

    time <- as.numeric(PP[["POSIXtime"]])
    tres <- diff(PP[["timeLimits.s"]])

    latitudes <- PP[["latitude"]][hh]
    longitudes <- PP[["longitude"]][hh]
    heights <- PP[["height"]][hh]

    return( c( EfieldF(vi,vicov,B,vipar0,...) , list( time=time , tres=tres ,
              latitudes=latitudes , longitudes=longitudes , heights=heights,
              B=B , date=PP[["date"]] , sites=PP[["sites"]] , POSIXtime=PP[["POSIXtime"]]
                                                 )
              )
           )
}
