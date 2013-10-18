acfscales <- function( sites , ranges , caltable )
    {
        #
        # interpolate calibration constants for measurements from given ranges at given sites using caltable
        #
        # I. Virtanen 2013
        #



        ndata <- dim(sites)[1]

        scales <- matrix( NA , nrow=ndata , ncol=6 )

        nsites <- dim(caltable)[1]
        
        diffs <- matrix( nrow=nsites , ncol=9)

        for( d in seq(ndata)){

            # TXfreq, TXlat, TXlon, TXaz, TXel, RXlat,RXlon,RXaz,RXel
            sd <- sites[ d , c(2:4,6:9,11:12) ]

            for( s in seq(nsites)){
                diffs[s,] <- caltable[s,c(7,1,2,5,6,8,9,12,13)] - sd
            }
            diffs[,4] <- diffs[,4]%%360
            diffs <- abs(diffs)
            diffs[is.na(diffs)] <- 0
            # with 90 degree elevation the azimuth angle may be anything
            if(sd[5]==90) diffs[4] <- 0
            if(sd[9]==90) diffs[8] <- 0
            
            # pict the correct site from table, take the first one if there are several acceptable ones
            # this will be replaced with the closest on in future
            site <- which(apply( diffs , FUN=function(x,y){all(x<y)} , MARGIN=1 , y=c( 1e7, .1, .1, .1, .1, .1, .1, .1, .1) ))[1]

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
