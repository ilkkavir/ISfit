acfscales <- function( sites )
    {
        #
        # acf scaling factors and site information for different sites
        #
        #
        # Columns of the scales matrix are:
        #
        #  1. fwhm of bore sight transmitter beam
        #  2. logical, is the transmitter a phased-array?
        #  3. fwhm of bore sight receiver beam
        #  4. logical, is the receiver a phased-array?
        #  5. logical, should the lagged products be complex-conjugated?
        #  6. acf scaling factors
        #
        # I. Virtanen 2013
        #



        nsites <- dim(sites)[1]

        scales <- matrix( NA , nrow=nsites , ncol=6 )


        for( n in seq( nsites ) ){
            if(any(is.na(sites[n,]))){
            }else{
                # Tromso TX
                if( abs( sum( sites[n,3:4] - radarSites()[["TRO"]][1:2] ) ) < .1 ){
                    scales[,1] <- 299792458/sites[,2]/ifelse(sites[,2]<250e6,60,32)*180/pi
                    scales[,2] <- FALSE
                }
                # Tromso RX
                if( abs( sum( sites[n,8:9] - radarSites()[["TRO"]][1:2] ) ) < .1 ){
                    scales[,3] <- 299792458/sites[,2]/ifelse(sites[,2]<250e6,60,32)*180/pi
                    scales[,4] <- FALSE
                    scales[,5] <- FALSE
                }
                # Kiruna RX
                if( abs( sum( sites[n,8:9] - radarSites()[["KIR"]][1:2] ) ) < .1 ){
                    scales[,3] <- 299792458/sites[,2]/32*180/pi
                    scales[,4] <- FALSE
                    scales[,5] <- FALSE
                }
                # Sodankyla RX
                if( abs( sum( sites[n,8:9] - radarSites()[["SOD"]][1:2] ) ) < .1 ){
                    scales[,3] <- 299792458/sites[,2]/32*180/pi
                    scales[,4] <- FALSE
                    scales[,5] <- FALSE
                }
                # Tromso monostatic
                if( abs( sum( sites[n,c(3,4,8,9)] - rep(radarSites()[["TRO"]][1:2],2) ) ) < .1 ){
                    scales[,6] <- 1
                }
                # Tromso - Kiruna bistatic
                if( abs( sum( sites[n,c(3,4,8,9)] - c(radarSites()[["TRO"]][1:2],radarSites()[["KIR"]][1:2]) ) ) < .1 ){
                    scales[,6] <- 1
                }
                # Tromso - Sodankyla bistatic
                if( abs( sum( sites[n,c(3,4,8,9)] - c(radarSites()[["TRO"]][1:2],radarSites()[["SOD"]][1:2]) ) ) < .1 ){
                    scales[,6] <- 1
                }
            }
        }
        
        return(scales)
    }
