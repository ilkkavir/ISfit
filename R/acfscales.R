acfscales <- function( sites , ranges , caltable )
    {
        #
        # interpolate calibration constants for measurements from given ranges at given sites using caltable
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



        ndata <- dim(sites)[1]

        scales <- matrix( NA , nrow=ndata , ncol=6 )

        # a bubble gum fix to make tri-static EISCAT analysis faster..
        if(is.null(caltable)){
            # this will give way too wide beam for the VHF, should add an if somewhere to tune the scale...
#            scales[,1] <- 299792458/sites[,2]/32*180/pi
            scales[,1] <- 299792458/sites[,2]/ifelse(sites[,2]<250e6,60,32)*180/pi
            scales[,2] <- FALSE
            scales[,3] <- 299792458/sites[,2]/ifelse(sites[,2]<250e6,60,32)*180/pi
#            scales[,3] <- 299792458/sites[,2]/32*180/pi
            scales[,4] <- FALSE
            scales[,5] <- FALSE
            scales[,6] <- 1
            return(scales)
        }

        nsites <- dim(caltable)[1]
        
        diffs <- matrix( ncol=nsites , nrow=9)

        for( d in seq(ndata)){

            # TXfreq, TXlat, TXlon, TXaz, TXel, RXlat,RXlon,RXaz,RXel
            sd <- sites[ d , c(2:4,6:9,11:12) ]

            diffs[,] <-  apply( caltable , FUN=function(x,y){ x[c(7,1,2,5,6,8,9,12,13)] - y }, y=sd , MARGIN=1 )

            diffs[4,] <- diffs[4,]%%360
            diffs <- abs(diffs)
            diffs[is.na(diffs)] <- 0
            # with 90 degree elevation the azimuth angle may be anything
            if(sd[5]==90) diffs[4,] <- 0
            if(sd[9]==90) diffs[8,] <- 0
            
            # pic the correct site from table, take the first one if there are several acceptable ones
            # this will be replaced with the closest on in future
            site <- which(apply( diffs , FUN=function(x,y){all(x<y)} , MARGIN=2 , y=c( 1e7, .1, .1, .1, .1, .1, .1, .1, .1) ))[1]

            # if there is value for only one range, simply pick it
            if (caltable[site,15]==1){
                scales[d,] <- caltable[site,c(3,4,10,11,14,17)]
            }else{ # otherwise we will interpolate
                r <- ranges[site]
                # all ranges from the relevant caltable row
                rtable <- caltable[site,16:(15+caltable[site,15])]

                rtableup <- rtable
                rtablup[rtable>r] <- NA
                rtabledown <- rtable
                rtabledown[rtable<=r] <- NA

                indup <- which.min(abs(rtableup-r))
                inddown <- which.min(abs(rtabledown-r))

                if(length(indup)==0) indup <- inddown
                if(length(inddown)==0) inddown <- indup

                if(indup==inddown){
                    scales[d,] <- caltable[site,c(3,4,10,11,14,15+caltable[site,15]+indup)]
                }else{
                    rup   <- rtable[indup]
                    rdown <- rtable[inddown]
                    sup   <- caltable[site,indup+15+caltable[site,15]]
                    sdown <- caltable[site,inddown+15+caltable[site,15]]
                    scales[d,] <- c( caltable[site,c( 3 , 4 , 10 , 11 , 15)], sdown + (r-rdown)/(rup-rdown)*(sup-sdown) )
                }
            }
        }

        return(scales)
    }
