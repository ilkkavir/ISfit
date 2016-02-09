EfieldF <- function(vi,vicov,B){
    #
    # Calculate F-region electric field from
    # ion velocity measurement(s).
    #
    # input:
    #   vi        an ion velocity vector in cartesian ENU
    #             coordinates [m/s] OR
    #             a matrix with one velocity vector on each line
    #   vicov     covariance matrix of a velocity vector OR
    #             a list of covariance matrices.
    #             length(vicov) == dim(vi)[1]
    #   B         magnetic field [nT] in the same coordinate system
    #             with vi. Either a vector, or a list of vectors
    #
    #  output
    #   E         Electric field estimate
    #             a 2D-vector in a coordinate system where
    #             x-axis points toward magnetic east, and
    #             y-axis toward magnetic north. Both perpendicular
    #                    to B. 
    #
    #
    # I. Virtanen 2016
    #


    # convert vectors into matrices and make sure vicov is a list
    if(is.vector(vi)) vi <- matrix(vi[1:3],ncol=3)
    if(is.vector(B)) B <- matrix(B[1:3],ncol=3)
    if(!is.list(vicov)) vicov <- list(vicov)


    # we will not need the parallel components, but I am keeping
    # them for testing purposes
    
    # a matrix for electric fields
    Efields <- matrix(nrow=dim(vi)[1],ncol=2)
    # an array for E-field covariances
    Ecov <- array(dim=c(dim(vi)[1],2,2))

    
    # fields and covariances in each gate
    # there must be a faster way to do this, but this is
    # easy to understand and fast enough for me
    for(hh in seq(1,dim(vi)[1])){
        # horizontal component of B
        Bhor <- c(B[hh,c(1,2)],0)
        names(Bhor) <- c('x','y','z')
        # geomagnetic east
        Bx <- radarPointings:::vectorProduct.cartesian(B[hh,],Bhor)
        Bx <- Bx/sqrt(sum(Bx**2))
        names(Bx) <- c('x','y','z')
        # geomagnetic north, perpendicular to B
        By <- radarPointings:::vectorProduct.cartesian(Bx,B[hh,])
        By <- By/sqrt(sum(By**2))
        names(By) <- c('x','y','z')
        # a unit vector parallel with B, but upwards
        Bz <- -B[hh,]/sqrt(sum(B[hh,]**2))
        names(Bz) <- c('x','y','z')
        # field intensity (the input  is in nT)
        Bstrength <- sqrt(sum(B[hh,]**2)) * 1e-9
        names(Bstrength) <- NULL
        # Ion velocity projection to Bx, By, Bz coordinates
        rotMat <- matrix(c(Bx,By,Bz),ncol=3,byrow=TRUE)
        ViBxyz <- rotMat%*%vi[hh,]
        names(ViBxyz) <- NULL
        # Covariance
        ViBxyzCov <- rotMat%*%vicov[[hh]]%*%t(rotMat)
        dimnames(ViBxyzCov) <- NULL
        # The E-field components and their covariances
        rotMat2 <- matrix(c(0,1,-1,0),byrow=TRUE,ncol=2)
        Efields[hh,] <- rotMat2%*%ViBxyz[c(1,2)] * Bstrength
        Ecov[hh,,] <- rotMat2%*%ViBxyzCov[1:2,1:2]%*%t(rotMat2) * Bstrength**2
        
    }

    # average over all heights. This could be part of the
    # previous loop, but let us keep things simple again..
    Mtmp <- c(0,0)
    Q <- matrix(0,ncol=2,nrow=2)
    for(hh in seq(1,dim(vi)[1])){
        Qh <- tryCatch( solve(Ecov[hh,,]), error=function(e){matrix(0,ncol=2,nrow=2)})
        Q <- Q + Qh
        Mtmp <- Mtmp + Qh%*%Efields[hh,]
    }

    Ecov <- tryCatch( solve(Q) , error=function(e){matrix(0,ncol=2,nrow=2)})
    E <- Ecov%*%Mtmp
        

    return(list(E=E,cov=Ecov))
    

    
}
    
