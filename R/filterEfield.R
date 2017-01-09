filterEfield <- function( Efield , tres=600 , startTime=NULL , ...){
    #
    # Post-integrate electric field estimates from EletricFieldsF to time resolution tres.
    #
    # INPUT:
    #  Efield output list from ElectricFieldsF
    #  tres  decimated time resolution in seconds
    #  startTime c(hour,minute,second) of the start time of the first integration period.
    #        the default is to start from 00:00 UT on the day of first data point
    #
    # OUTPUT:
    #  a list of filtered electric field data, identical elements with the ElectricFieldsF output list,
    #
    # I. Virtanen 2016


    # start time of the first integration period
    if(is.null(startTime)){
        firstTime <- floor(Efield[["time"]][1]/tres) * tres
        if(firstTime == Efield[["time"]][1]) firstTime <- firstTime - tres
    }else{
        firstTime <- floor(Efield[["time"]][1]/86400)*86400 + sum( startTime * c( 3600 , 60 , 1 ) )
    }

    # number of integration periods
    niper <- ceiling( ( max(Efield[["time"]]) - firstTime ) / tres )


    # filter Efield and Ecov
    EfieldF <- array(dim=c(niper,2))
    EcovF <- array(dim=c(niper,2,2))
    timeF <- rep(0,niper)
    for(k in seq(niper)){
        tlims <- c((k-1),k)*tres + firstTime
        dinds <- which((Efield[["time"]]>tlims[1])&(Efield[["time"]]<=tlims[2]))
        Q <- matrix(0,ncol=2,nrow=2)
        Mtmp <- c(0,0)
        for(n in dinds){
            Qn <- tryCatch( solve(Efield[["Ecov"]][n,,]) , error=function(e){matrix(0,ncol=2,nrow=2)})
            if((!(any(is.na(Qn))))&(!any(is.na(Efield[["E"]][n,])))){
                Q <- Q + Qn
                Mtmp <- Mtmp + Qn%*%Efield[["E"]][n,]
            }
        }
        EcovF[k,,] <- tryCatch( solve(Q) , error=function(e){matrix(0,ncol=2,nrow=2)})
        EfieldF[k,] <- EcovF[k,,]%*%Mtmp
        timeF[k] <- Efield[["time"]][max(dinds)]
    }

    return(list(E=EfieldF,Ecov=EcovF,time=timeF,tres=tres))


}
