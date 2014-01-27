averageVelocities <- function( ddir='.' , ofile=NULL , ascii=FALSE )
    {

        # average ion vector velocity components over all measured heights and return
        # the averaged vectors in both geographic and geomagnetic coordinate systems

        # read the data files, this also puts NA's to the height gates which are not measured from at
        # least three different directions
        dlist <- readPP.3D( ddir )

        # number of integration periods
        nt <- length( dlist[["time_sec"]] )

        # matrices for velocity components and their error estimates
        Vion <- matrix(ncol=6,nrow=nt)
        Verr <- Vion

        # dimnames for the matrices
        dnames <- list( sapply( dlist[["POSIXtime"]] , format , format="%Y-%m-%d %H:%M:%S" ) , c('Vix','Viy','Viz','ViBx','ViBy','ViB') )
        dimnames(Vion) <- dnames
        dimnames(Verr) <- dnames

        # variance-weighted averages of velocity vector components and their
        # standard deviations
        Verr[,'Vix'] <- 1 / colSums( 1 / dlist[["std"]][,'Vix',]**2 , na.rm=T )
        Vion[,'Vix'] <- colSums( dlist[["param"]][,'Vix',] / dlist[["std"]][,'Vix',]**2,na.rm=T) * Verr[,'Vix']
        
        Verr[,'Viy'] <- 1 / colSums( 1 / dlist[["std"]][,'Viy',]**2 , na.rm=T )
        Vion[,'Viy'] <- colSums( dlist[["param"]][,'Viy',] / dlist[["std"]][,'Viy',]**2,na.rm=T) * Verr[,'Viy']
        
        Verr[,'Viz'] <- 1 / colSums( 1 / dlist[["std"]][,'Viz',]**2 , na.rm=T )
        Vion[,'Viz'] <- colSums( dlist[["param"]][,'Viz',] / dlist[["std"]][,'Viz',]**2,na.rm=T) * Verr[,'Viz']
        
        Verr[,'ViBx'] <- 1 / colSums( 1 / dlist[["std"]][,'ViBx',]**2 , na.rm=T )
        Vion[,'ViBx'] <- colSums( dlist[["param"]][,'ViBx',] / dlist[["std"]][,'ViBx',]**2,na.rm=T) * Verr[,'ViBx']
        
        Verr[,'ViBy'] <- 1 / colSums( 1 / dlist[["std"]][,'ViBy',]**2 , na.rm=T )
        Vion[,'ViBy'] <- colSums( dlist[["param"]][,'ViBy',] / dlist[["std"]][,'ViBy',]**2,na.rm=T) * Verr[,'ViBy']
        
        Verr[,'ViB'] <- 1 / colSums( 1 / dlist[["std"]][,'ViB',]**2 , na.rm=T )
        Vion[,'ViB'] <- colSums( dlist[["param"]][,'ViB',] / dlist[["std"]][,'ViB',]**2,na.rm=T) * Verr[,'ViB']

        Verr <- sqrt( Verr )

        outlist <- list( Vion=Vion , Verr=Verr , time_sec=dlist[["time_sec"]] , POSIXtime=dlist[["POSIXtime"]] , date=dlist[["date"]] , note="Vix=geographic east, Viy=geographic north, Viz=upwards, ViBx=geomagnetic east, ViBy=geomagnetic north, ViB=field-aligned (downwards in northern hemisphere!)")

        # optionally write the data to file
        if( !is.null( ofile ) ){
            # ascii file option for portability
            if(ascii){
                ff <- file( paste( ofile , '.txt' , sep='' ) , 'wt' )
                cat("##                                                               \n" , file=ff, append=FALSE )
                cat("## Vi = ion velocity (m/s)                                       \n" , file=ff )
                cat("## std= standard deviation of ion velocity (m/s)                 \n" , file=ff )
                cat("##                                                               \n" , file=ff )
                cat("## Timestamps are integration period end times in UTC            \n" , file=ff )
                cat("##                                                               \n" , file=ff )
                cat("## x  = geographic east                                          \n" , file=ff )
                cat("## y  = geographic north                                         \n" , file=ff )
                cat("## z  = upwards                                                  \n" , file=ff )
                cat("## Bx = geomagnetic east                                         \n" , file=ff )
                cat("## By = geomagnetic west                                         \n" , file=ff )
                cat("## B  = along magnetic field (downwards in northern hemisphere!) \n" , file=ff )
                cat("##                                                               \n" , file=ff )
                cat( "         Time (UTC)        Vix        Viy        Viz       ViBx       ViBy        ViB       stdx       stdy       stdz      stdBx      stdBy       stdB\n" , file=ff )
                for( k in seq(nt) ){
                    cat( dnames[[1]][k] , sprintf( "%10.2f" , Vion[k,] ) , sprintf( "%10.2f" , Verr[k,] ) , '\n' ,file=ff )
                }
                close(ff)
            # otherwise Rdata
            }else{
                save( Vion=Vion , Verr=Verr , time_sec=dlist[["time_sec"]] , POSIXtime=dlist[["POSIXtime"]] , date=dlist[["date"]] ,
                     note="Vix=geographic east, Viy=geographic north, Viz=upwards, ViBx=geomagnetic east, ViBy=geomagnetic north, ViB=field-aligned (downwards in northern hemisphere!)" , file=paste( ofile , '.Rdata' , sep='' ) )
            }
            # return the output invisibly
            return( invisible( outlist ) )
        }else{
            # visible return
            return( outlist )
        }

    }
