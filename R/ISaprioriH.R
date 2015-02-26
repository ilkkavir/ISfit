ISaprioriH <- function( aprioriParam ,  nIon , absCalib=FALSE , TiIsotropic=FALSE , TeIsotropic=FALSE, refSite=1 , siteScales=NULL , h=300 , hTeTi=120 , hTi=100 , hVi, hColl=c(0,0) , ... )
    {
        #
        #
        # Apriori function with height dependence
        #
        # INPUT:
        #  aprioriParam    apriori parameter values
        #  nIon            number of ion masses
        #  absCalib        Logical, all site scales are fixed to unity with small variance if absCalib==TRUE and
        #                  siteScales==NULL
        #  TiIsotropic     TRUE if isotropic Ti is assumed, FALSE for bimaxwellian ion velocity distribution
        #  TeIsotropic     TRUE if isotropic Te is assumed, FALSE for bimaxwellian electron velocity distribution
        #  refSite         reference site, whose scale is fixed to unity with small variance
        #  siteScales      a matrix of site scales and their variances or NULL
        #  h               height in km
        #
        #  ...             arbitrary parameters to be passed forward to other functions, mainly for compatability reasons
        #
        # OUTPUT:
        #  aprioriTheory      apriori theory matrix
        #  aprioriMeas        apriori "measurements"
        #  invAprioriCovar    inverse of apriori covariance matrix
        #
        #  I. Virtanen 2012, 2013
        #

#        # height limits, these could be optional input arguments as well
#        hTi <- 100     # Ti from model below this height
#        hTeTi <- 120   # Te=Ti below this height
            
        
        # length of the parameter vector
        nPar                         <- length(aprioriParam)

        # number of imaginary apriori "measurements"
        nApriori                     <- ( nPar + 4 )

        # apriori theory matrix
        aprioriTheory                <- matrix( 0 , nrow=nApriori , ncol=nPar )

        # apriori measurement vector
        aprioriMeas                  <- vector(mode='numeric',length=nApriori)

        # the apriori covariance matrix will be diagonal, so we begin with
        # a vector of standard deviations, which is easier.
        aprioriStd                   <- vector(mode='numeric',length=nApriori)

        # apriori parameter values
        aprioriTheory[1:nPar,1:nPar] <- diag(rep(1,nPar))

        aprioriMeas[1:nPar]          <- aprioriParam

        aprioriStd[1]                <- 1e4                                        # electron density
        aprioriStd[2]                <- ifelse(h<hTi,1e-3,1)                       # parallel ion temperature
        aprioriStd[3]                <- 1                                          # perpendicular ion temperature
        aprioriStd[4]                <- 1                                          # parallel electron temperature
        aprioriStd[5]                <- 1                                          # perpendicular electron temperature
        aprioriStd[6]                <- ifelse((h>hColl[1])&(h<hColl[2]),.1,1e-3)  # ion-neutral collision frequency
        aprioriStd[7]                <- ifelse(h<hVi,.1,10)                        # ion velocity, x-component
        aprioriStd[8]                <- ifelse(h<hVi,.1,10)                        # ion velocity, y-component
        aprioriStd[9]                <- ifelse(h<hVi,.1,10)                        # ion velocity, z-component
        aprioriStd[10:(9+nIon)]      <- 1e-3                                       # ion abundances


        # remove model information about perpendicular temperatures
        aprioriTheory[3,] <- 0
        aprioriTheory[5,] <- 0
        aprioriMeas[c(3,5)] <- 0



        if(absCalib){
            aprioriStd[(nIon+10):length(aprioriParam)] <- 1e-3 # fix all sites to the same ACF scale
        }else{
            aprioriStd[(nIon+10):length(aprioriParam)] <- 1   # allow scaling for other sites
        }
        if(!is.null(siteScales)){
            if(!is.matrix(siteScales)) siteScales <- matrix(siteScales,nrow=1)
            ssinds <- which(!is.na(rowSums(siteScales)))
            aprioriMeas[ssinds+nIon+9] <- siteScales[ssinds,1]  # user-given scaling factors
            aprioriStd[ssinds+nIon+9] <- siteScales[ssinds,2]
        }

        aprioriStd[nIon+9+refSite]     <- 1e-3                 # do not allow scaling at the reference site

        # force certain parameter differences close to zero
        curRow                         <- nPar + 1

        # electron temperature anisotropy
        aprioriTheory[curRow,c(4,5)]   <- c(1,-1)
        aprioriMeas[curRow]            <- 0
        if(TeIsotropic){
            aprioriStd[curRow]             <- 1e-3
        }else{
            aprioriStd[curRow]             <- 1
        }
        curRow                         <- curRow + 1

        # ion temperature anisotropy
        aprioriTheory[curRow,c(2,3)]   <- c(1,-1)
        aprioriMeas[curRow]            <- 0
        if(TiIsotropic){
            aprioriStd[curRow]             <- 1e-3
        }else{
            aprioriStd[curRow]             <- 1
        }
        curRow                         <- curRow + 1

        # Sum of ion abundances must be unity
        aprioriTheory[curRow,10:(nIon+9)] <- 1
        aprioriMeas[curRow] <- 1
        aprioriStd[curRow] <- 1e-3

        # Te=Ti below hTeTi
        curRow                         <- curRow + 1
        aprioriTheory[curRow,c(2,4)] <- c(1,-1)
        aprioriMeas[curRow] <- 0
        aprioriStd[curRow] <- ifelse(h<hTeTi,1e-3,10)


        return(list(aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas))

    }
