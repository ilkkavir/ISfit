ISaprioriBAFIMsmoother <- function( PP , height , resFile , resdirForward , resdirReverse , updateFile=FALSE , returnParams=TRUE, ... )
{

        # THIS FUNCTION IS UNDER DEVELOPMENT!!!


        # should read the predictions from forward & backward in time runs and then fit using these as a prior
    
        #
        #
        # Prior model for plasma parameter fits by means of Bayesian Filtering
        #
        # INPUT:
        #  PP              fit output list from the previous integration period (an empty list in the first integration period)
        #  date            measurement time as c(year,month,day,hour,minute,seconds)
        #  latitude        geodetic latitudes of the measurement volumes (deg north)
        #  longitude       geodetic longitudes of the measurement volumes (deg east)
        #  height               heights in km
        #  nIon            number of ion masses
        #  absCalib        Logical, all site scales are fixed to unity with small variance if absCalib==TRUE and
        #                  siteScales==NULL
        #  TiIsotropic     TRUE if isotropic Ti is assumed, FALSE for bimaxwellian ion velocity distribution
        #  TeIsotropic     TRUE if isotropic Te is assumed, FALSE for bimaxwellian electron velocity distribution
        #  refSite         reference site, whose scale is fixed to unity with small variance
        #  siteScales      a matrix of site scales and their variances or NULL
        #  hTeTi           Te=Ti below hTeTi [km]
        #  B               magnetic field (direction). The default is considered as missing value
        #  ViPar0          logical, force field-aligned ion velocity to zero
        #  nCores          number of cpu cores to use in forks
        #
        #  ...             arbitrary parameters to be passed forward to other functions, mainly for compatability reasons
        #
        # OUTPUT:
        #  aprioriTheory      apriori theory matrix
        #  aprioriMeas        apriori "measurements"
        #  invAprioriCovar    inverse of apriori covariance matrix
        #
        #  I. Virtanen 2012, 2013, 2023

    PPorig <- PP
    
    # call ISaprioriBAFIM with the same arguments to smooth the previous fit in altitude
    PPrcorr <- ISaprioriBAFIM( PP=PP , height=height ,  updateFile=T , returnParams=T , ... )

        # read the two filter results from files
    load(file.path(resdirForward,resFile))
    PPfw <- PP
    load(file.path(resdirReverse,resFile))
    PPrev <- PP

    apriorilist <- list()
    nh <- length(height)
    for(h in seq(nh)){
        apriorilist[[h]] <- list(
            aprioriTheory = rbind( PPfw$apriori[[h]]$aprioriTheory , PPrev$apriori[[h]]$aprioriTheory ) ,
            invAprioriCovar = cbind( rbind( PPfw$apriori[[h]]$invAprioriCovar , PPfw$apriori[[h]]$invAprioriCovar * 0) ,
                                    rbind( PPrev$apriori[[h]]$invAprioriCovar *0 , PPrev$apriori[[h]]$invAprioriCovar ) ) ,
            aprioriMeas = c( PPfw$apriori[[h]]$aprioriMeas , PPrev$apriori[[h]]$aprioriMeas ) ,
            limitParam = PPfw$apriori[[h]]$limitParam ,
            parScales = PPfw$apriori[[h]]$parScales ,
            aprioriParam = ( PPfw$apriori[[h]]$aprioriParam + PPrev$apriori[[h]]$aprioriParam ) /  2 , # this is just copied in ISfit.3D
            mIon = PPfw$apriori[[h]]$mIon ,
            nIon = PPfw$apriori[[h]]$nIon
        )
    }

    PP <- PPorig
    
    if(length(PP)>0){
            
        PP$paramSmoother <- PP$param
        PP$stdSmoother <- PP$std
        PP$covarSmoother <- PP$covar
        
        PP$param  <- PP$paramSmootherRcorr <- PPrcorr$paramRcorr
        PP$std  <- PP$paramSmootherRcorr <- PPrcorr$stdRcorr
        PP$covar <- PP$covarSmootherRcorr <- PPrcorr$covarRcorr
        
        # overwrite the output file with the updated copy
        save(PP,file=file.path(PP$resDir,PP$resFile))
    }
    
    return(apriorilist)
    
}
