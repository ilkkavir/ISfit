EfieldFPP <- function(PP,hmin=200,hmax=400,vipar0=FALSE){
    # F-region electric field with ISfit output list
    # as input

    hh <- which( (PP[["height"]]>=hmin) & (PP[["height"]] <= 400) )

    vi <- PP[["param"]][hh,c('Vix','Viy','Viz')]
    vicov <- lapply(PP[["covar"]][hh],FUN=function(x){return(x[c('Vix','Viy','Viz'),c('Vix','Viy','Viz')])})
    B <- PP[["B"]][hh,]

    return(EfieldF(vi,vicov,B,vipar0))
}
