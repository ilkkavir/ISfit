siteCalib <- function( ddir='.' , recursive=FALSE )

    {


        # list data files
        dfiles <- dir( ddir , pattern="[[:digit:]]PP.Rdata" , full.names=T , recursive=recursive)

        # check that the dfiles vector is not empty
        if( length(dfiles)==0) stop("No data files")

        # all files should have identical height gates, load the first one 
        load(dfiles[1])

        # number of height gates
        nh <- length(PP[["height"]])
    
        # initialize with 1000 rows, expand later if necessary
        sites <- matrix(ncol=(12+nh*2),nrow=1000)

        # total number of sites found
        nstot <- 0
    
        for( f in dfiles){
        
            # load a new file
            load(f)

            # number of sites in this files
            ns <- dim(PP[["sites"]])[1]

            # round site locations to 0.1 degree  and pointings to 1 degree accuracy
            sscale <- rep(c(1,1,10,10,1,1,1,10,10,1,1,1),each=ns)
            sitesf <- round(PP[["sites"]]*sscale)/sscale
            sitesf[sitesf[,7]==90,6]<-90
            sitesf[sitesf[,12]==90,11]<-90

            for( s in seq(ns) ){
                # try to find this site from the existing site matrix
                sind <- which( apply( sites[,2:12] , FUN=function(x,y){all(x==y)} , MARGIN=1 , y=sitesf[s,2:12] ) )

                # if this is a new site
                if(length(sind)==0){

                    # increment the sites counter
                    nstot <- nstot + 1

                    # enlarge the sites matrix if necessary
                    if(nstot>dim(sites)[1]){
                        sitestmp <- sites
                        sites <- matrix(ncol=(12+nh),nrow=(dim(sitestmp)[1]+1000))
                        sites[1:dim(sitestmp)[1],] <- sitestmp
                    }

                    # initialize the new row
                    sites[nstot,] <- 0
                    sites[nstot,1:12] <- sitesf[s,]
                    sind <- nstot
                }

                # accept only those points that actually have contribution from this site
                hinds <- rep(0,nh)
                for(h in seq(length(PP[["contribSites"]]))){
                    hinds[h] <- any( PP[["contribSites"]][[h]]==s)
                }
            
                # ACF scales at all heights
                srow <- c(PP[["param"]][,paste('Site',sitesf[s,1],sep='')]/PP[["std"]][,paste('Site',sitesf[s,1],sep='')]**2,1/PP[["std"]][,paste('Site',sitesf[s,1],sep='')]**2) * hinds

                # replace missing values with zeros, this is ok since we are
                # summing information
                srow[is.na(srow)] <- 0
                srow[is.infinite(srow)] <- 0

                # add srow to the correct site
                sites[sind,13:(12+2*nh)] <- sites[sind,13:(12+2*nh)] + srow
            }
            cat('                                         \r',f)
        }
        cat('\n')

        # final variance
        sites[,(13+nh):(12+2*nh)] <- 1/sites[,(13+nh):(12+2*nh)]

        # final scales
        sites[,13:(12+nh)]        <- sites[,13:(12+nh)] * sites[,(13+nh):(12+2*nh)]

        # generate new site indices
        sites[1:nstot,1] <- seq(nstot)

        # return the scales for all sites
        return(sites[1:nstot,])
        
    }
