filterNwind <- function( nWind , tres=600 , startTime=NULL , ...){
    #
    # Post-integrate neutral wind estimates from NeutralWinds to time resolution tres.
    #
    # INPUT:
    #  nWind output list from NeutralWinds
    #  tres  decimated time resolution in seconds
    #  startTime c(hour,minute,second) of the start time of the first integration period.
    #        the default is to start from 00:00 UT on the day of first data point
    #
    # OUTPUT:
    #  a list of filtered neutral wind data, identical elements with the NeutralWinds output list,
    #  except E and Ecov stripped off
    #
    # I. Virtanen 2016


    # start time of the first integration period
    if(is.null(startTime)){
        firstTime <- floor(nWind[["time"]][1]/tres) * tres
        if(firstTime == nWind[["time"]][1]) firstTime <- firstTime - tres
    }else{
        firstTime <- floor(nWind[["time"]][1]/86400)*86400 + sum( startTime * c( 3600 , 60 , 1 ) )
    }

    # number of integration periods
    niper <- ceiling( ( max(nWind[["time"]]) - firstTime ) / tres )

    # number of height gates
    nEgates <- dim(nWind[["nWind"]])[2]


    # filter nWind and nWcov
    nWindF <- array(dim=c(niper,nEgates,3))
    nWcovF <- array(dim=c(niper,nEgates,3,3))
    nWindBF <- array(dim=c(niper,nEgates,3))
    nWcovBF <- array(dim=c(niper,nEgates,3,3))
    timeF <- rep(0,niper)
    for(k in seq(niper)){
        tlims <- c((k-1),k)*tres + firstTime
        dinds <- which((nWind[["time"]]>tlims[1])&(nWind[["time"]]<=tlims[2]))
        for(hh in seq(nEgates)){
            Q <- matrix(0,ncol=3,nrow=3)
            QB <- matrix(0,ncol=3,nrow=3)
            Mtmp <- c(0,0,0)
            MtmpB <- c(0,0,0)
            for(n in dinds){
                Qn <- tryCatch( solve(nWind[["cov"]][n,hh,,]) , error=function(e){matrix(0,ncol=3,nrow=3)})
                QnB <- tryCatch( solve(nWind[["covB"]][n,hh,,]) , error=function(e){matrix(0,ncol=3,nrow=3)})
                if((!(any(is.na(Qn))))&(!any(is.na(nWind[["nWind"]][n,hh,])))){
                    Q <- Q + Qn
                    Mtmp <- Mtmp + Qn%*%nWind[["nWind"]][n,hh,]
                }
                if((!(any(is.na(QnB))))&(!any(is.na(nWind[["nWindB"]][n,hh,])))){
                    QB <- QB + QnB
                    MtmpB <- MtmpB + QnB%*%nWind[["nWindB"]][n,hh,]
                }
            }
            nWcovF[k,hh,,] <- tryCatch( solve(Q) , error=function(e){matrix(0,ncol=3,nrow=3)})
            nWindF[k,hh,] <- nWcovF[k,hh,,]%*%Mtmp
            nWcovBF[k,hh,,] <- tryCatch( solve(QB) , error=function(e){matrix(0,ncol=3,nrow=3)})
            nWindBF[k,hh,] <- nWcovBF[k,hh,,]%*%MtmpB
            timeF[k] <- nWind[["time"]][max(dinds)]
        }
    }

    return(list(nWind=nWindF,cov=nWcovF,nWindB=nWindBF,covB=nWcovBF,time=timeF,height=nWind[["height"]]))


}
