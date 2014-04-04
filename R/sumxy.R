sumxy <- function( xxdir , yydir , odir , yyscale=1)
    {
        #
        # sum x and y polarization ACFs in KAIRA analysis
        #

        # list files in xxdir
        xxfiles <- dir(xxdir,full.names=FALSE,pattern='LP.Rdata')
        # list files in yydir
        yyfiles <- dir(yydir,full.names=FALSE,pattern='LP.Rdata')

        # common files
        xyfiles <- intersect( xxfiles , yyfiles )

        # stop if no common files
        if(length(xyfiles)==0) stop(paste("No matching files names in",xxdir,"and",yydir))

        dir.create(odir,recursive=TRUE,showWarnings=FALSE)
        for( fn in xyfiles){
            load(file.path(yydir,fn))
            yyACF <- ACF
            load(file.path(xxdir,fn))
            ACF$ACF <- ACF$ACF + yyACF$ACF * yyscale
            ACF$var <- ACF$var + yyACF$var * yyscale**2
            save(ACF,file=file.path(odir,fn))
            cat('\r',fn)
        }

        invisible(length(xyfiles))
    }
