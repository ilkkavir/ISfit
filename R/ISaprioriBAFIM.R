ISaprioriBAFIM <- function( PP , date , latitude , longitude , height , nSite ,  nIon , absCalib=FALSE , TiIsotropic=FALSE , TeIsotropic=FALSE, refSite=1 , siteScales=NULL , hTeTi=100 , B=c(0,0,0) , ViPar0=FALSE , nCores=1, BAFIMpar=list(Ne=c(0,Inf,0.05,2.5e11),Ti=c(80,Inf,0.1,30),Te=c(100,Inf,0.1,30),Coll=c(0,0,0,0),Vi=c(80,Inf,0.05,2.5),Op=c(150,500,0.05,0.01)), ... )
    {
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
                                        #


        # THIS function has not been completed yet!

        nh <- length(height)

        # IRI parameters (would it be enough to call this with just one lat and lon?
#        IRIlist <- mclapply(seq(nh) , FUN=iriParamsParFun , date=date,latitude=latitude,longitude=longitude,height=height,fitGate=rep(T,nh) , okData=rep(T,nh) , mc.cores=nCores)
        IRIpar <- iriParams( time=date,latitude=mean(latitude),longitude=mean(longitude),heights=height) # IS this accurate enough for low-elevation measurements?


        # parameter value limits
        parLimits      <- ISparamLimits(3,nSite)


        
        apriorilist <- list()

        aprioriIRI <- list()
        for(h in seq(nh)){



            ################# IRI parameters ###################################


            
            # parameters from iri model
            ptmp <- IRIpar[,h]
            
            # an approximation for NO+-neutral colllision frequency (Schunk & Walker, Planet. Space Sci., 1971)
            # This is approximately true for all ions, because ion density is much smaller than neutral density
            ioncoll        <- sum( ionNeutralCollisionFrequency(ptmp)['NO+',] )


            # initial plasma parameter values
            cH <- max(ptmp['H+'],0)
            cO <- max(ptmp['O+'],0)
            cM <- max(sum(ptmp[c('NO+','O2+','cluster')]),0)
            cTot <- cH + cO + cM
            if(cTot<1e7){
                if(h<150){
                    cH <- 0
                    cO <- 0
                    cM <- 1
                    cTot <- 1
                }else{
                    cM <- 0
                    cTot <- cO + cH
                }
            }
                          
            parInit <- pmax( c( ptmp['e-'] , ptmp['Ti'] , ptmp['Ti'], ptmp['Te'] , ptmp['Te'] , ioncoll , 0 , 0 , 0 , cM/cTot , cO/cTot , cH/cTot , rep(1,nSite) ) , 0 )
            
            
            parInit[1]     <- max(parInit[1],1e9)
            
            mIon <- c(30.5,16.0,1)


            
            # parameter scaling factors
            parScales      <- ISparamScales(parInit,3)
            
            # scale the initial parameter values
            aprioriParam      <- scaleParams( parInit , parScales , inverse=F)

           # scale the parameter limits
            limitParam     <- parLimits
            limitParam[1,] <- scaleParams(parLimits[1,] , parScales , inverse=F)
            limitParam[2,] <- scaleParams(parLimits[2,] , parScales , inverse=F)

            aprioriIRI[[h]] <- list(aprioriParam=aprioriParam,limitParam=limitParam,parScales=parScales)
        }







        ############### smooth the plasma parameter profiles in altitude ##########################
        # Plasma scale heights
        kB <- 1.380649e-23
        amu <- 1.66053907e-27
        IRImol <- colSums(IRIpar[c('NO+','O2+','cluster'),])
        IRItot <- IRImol + IRIpar['O+',] + IRIpar['H+',]
        H <- kB * IRIpar['Ti',] * IRIpar['Te',] / 2 / ( amu * ( IRIpar['H+',]/IRItot + 16*IRIpar['O+',]/IRItot + 30.5*IRImol/IRItot ) * 9.82 * ( 6372/(6372+height) )**2 )


        
        if (length(PP)>0){
            dt <- as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - PP$time_sec
            hsAlt <- H/1000 * sqrt(dt)
            # replace unrealistic values and failed fits with the previous predictions
            for(h in seq(nh)){
                okfit <- TRUE
                if(PP$status[h] | PP$chisqr[h]>100 | any(PP$param[h,] < parLimits[1,]) | any(PP$param[h,] > parLimits[2,] ) ){
                    okfit <- FALSE
                }

#### continue from here!! we need to replace the bad points and then do the actual smoothing in altitude

                
                if(!okfit){
#                    PP$param[h,1:(9+nIon)] <- scaleParam(PP$apriori[[h]]$aprioriParam,)
                }

            }
            
            if(nh>2){
                
                
                
                
            }
        }
















        

 
        for(h in seq(nh)){
            

            aprioriParam <- aprioriIRI[[h]]$aprioriParam
            limitParam <- aprioriIRI[[h]]$limitParam
            parScales <- aprioriIRI[[h]]$parScales

            # length of the parameter vector
            nPar                         <- length(aprioriParam)

            # number of imaginary apriori "measurements"
            nApriori                     <- ( nPar + 5 )

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
            aprioriStd[2]                <- ifelse(height[h]<BAFIMpar$Ti[1],1e-3,2)                       # parallel ion temperature
            aprioriStd[3]                <- 2                                          # perpendicular ion temperature
            aprioriStd[4]                <- 2                                          # parallel electron temperature
            aprioriStd[5]                <- 2                                          # perpendicular electron temperature
            aprioriStd[6]                <- ifelse((height[h]>BAFIMpar$Coll[1])&(height[h]<BAFIMpar$Coll[2]),1,1e-3)   # ion-neutral collision frequency
            aprioriStd[7]                <- ifelse(height[h]<BAFIMpar$Vi[1],.1,10)                        # ion velocity, x-component
            aprioriStd[8]                <- ifelse(height[h]<BAFIMpar$Vi[1],.1,10)                        # ion velocity, y-component
            aprioriStd[9]                <- ifelse(height[h]<BAFIMpar$Vi[1],.1,10)                        # ion velocity, z-component
            aprioriStd[10:(9+nIon)]      <- 1e-3                                       # ion abundances


            # remove model information about perpendicular temperatures
            aprioriTheory[3,] <- 0
            aprioriTheory[5,] <- 0
            aprioriMeas[c(3,5)] <- 0




            ####################### replace the IRI predictions with range-smoothed profiles from the previous fit where appropriate ##


            # replace the IRI value with results from the previous time step (or from the previous prior if the failed)
            # and the standard deviations with that from the previous fit + the process noise
            # 
            # NOTE: this is a development version, in the final version we will smooth the profiles first!
            if(length(PP)>0){
                dt <- as.double(ISOdate(date[1],date[2],date[3],date[4],date[5],date[6])) - PP$time_sec
                processStdScale <- scaleParams(c( BAFIMpar$Ne[4] , BAFIMpar$Ti[4] , BAFIMpar$Ti[4] ,BAFIMpar$Te[4] , BAFIMpar$Te[4] , BAFIMpar$Coll[4] , BAFIMpar$Vi[4] , BAFIMpar$Vi[4] , BAFIMpar$Vi[4] , BAFIMpar$Op[4] ),parScales[1:10],inverse=F)

                if (!PP$status[h]){
                    parFit <- scaleParams(PP$param[h,],parScales,inverse=F)
                    stdFit <- scaleParams(PP$std[h,],parScales,inverse=F)
                }else{
                    parFit <- PP$apriori[[h]]$aprioriParam
                    stdFit <- 1/sqrt(diag(PP$apriori[[h]]$invAprioriCovar))
                }
                if (height[h]>=BAFIMpar$Ne[1]){
                    aprioriMeas[1] <- aprioriParam[1] <- parFit[1]
                    aprioriStd[1] <- stdFit[1] + processStdScale[1]*sqrt(dt)
                }
                if (height[h]>=BAFIMpar$Ti[1]){
                    aprioriMeas[2] <- aprioriParam[2] <- parFit[2]
                    aprioriStd[2] <- stdFit[2] + processStdScale[2]*sqrt(dt)
                    aprioriMeas[3] <- aprioriParam[3] <- parFit[3]
                    aprioriStd[3] <- stdFit[3] + processStdScale[3]*sqrt(dt)
                }
                if (height[h]>=BAFIMpar$Te[1]){
                    aprioriMeas[4] <- aprioriParam[4] <- parFit[4]
                    aprioriStd[4] <- stdFit[4] + processStdScale[4]*sqrt(dt)
                    aprioriMeas[5] <- aprioriParam[5] <- parFit[5]
                    aprioriStd[5] <- stdFit[5] + processStdScale[5]*sqrt(dt)
                }
                if (height[h]>=BAFIMpar$Vi[1]){
                    aprioriMeas[7] <- aprioriParam[7] <- parFit[7]
                    aprioriStd[7] <- stdFit[7] + processStdScale[7]*sqrt(dt)
                    aprioriMeas[8] <- aprioriParam[8] <- parFit[8]
                    aprioriStd[8] <- stdFit[8] + processStdScale[8]*sqrt(dt)
                    aprioriMeas[9] <- aprioriParam[9] <- parFit[9]
                    aprioriStd[9] <- stdFit[9] + processStdScale[1]*sqrt(dt)
                }
            }





            
            
            if(absCalib){
                aprioriStd[(nIon+10):length(aprioriParam)] <- 1e-3 # fix all sites to the same ACF scale
            }else{
                aprioriStd[(nIon+10):length(aprioriParam)] <- 1   # allow scaling for other sites
            }
            if(!is.null(siteScales)){
                if(!is.matrix(siteScales)) siteScales <- matrix(siteScales,nrow=1)
                ssinds <- which(!is.na(rowSums(siteScales)))
                aprioriMeas[ssinds+nIon+9] <- siteScales[ssinds,1]  # user-given scaling factors
                if(absCalib){
                    aprioriStd[ssinds+nIon+9] <- siteScales[ssinds,2]
                }
            }
            
            aprioriStd[nIon+9+refSite]     <- 1e-3                 # do not allow scaling at the reference site

            # force certain parameter differences close to zero
            curRow                         <- nPar + 1



            # the temperature ansitropies somewhat diffcult this way,
            # it would perhaps be better to fit the field-aligned temperature
            # and the difference Tperp - Tpar. This will require changes in a number
            # of places but could be worth it...

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
            aprioriStd[curRow] <- ifelse(height[h]<hTeTi,1e-3,10)
            if(height[h]<hTeTi){
                aprioriTheory[4,] <- 0
                aprioriTheory[5,] <- 0
                aprioriMeas[c(4,5)] <- 0
            }
            
            # optional ViPar=0
            curRow                         <- curRow + 1
            aprioriTheory[curRow,c(7,8,9)] <- B[h,]/sum(sqrt(B[h,]^2))
            aprioriMeas[curRow] <- 0
            aprioriStd[curRow] <- ifelse(ViPar0&all(B[h,]!=0),1e-3,10)

            apriorilist[[h]] <- list(aprioriParam=aprioriParam,aprioriTheory=aprioriTheory,invAprioriCovar=diag(1/aprioriStd**2),aprioriMeas=aprioriMeas,limitParam=limitParam,parScales=parScales,mIon=mIon,nIon)
        }

        return(apriorilist)

    }
