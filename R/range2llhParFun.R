range2llhParFun <- function(dind,ran,sites,sinds){



    llh <- NA
    
    if(!is.na(ran[dind])){
        llh <- range2llh( r=ran[dind] , llhT=sites[sinds[dind],3:5] , llhR=sites[sinds[dind],8:10] , azelT=sites[sinds[dind],6:7])
    }

    return(llh)
    
}
