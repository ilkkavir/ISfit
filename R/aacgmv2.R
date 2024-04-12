aacgmv2 <- function(inlat,inlon,height,date,methcode='G2A'){
    
    library(reticulate)
    Sys.setenv(RETICULATE_AUTOCREATE_PACKAGE_VENV="no")
    aacgmfile <- system.file('python','runaacgmv2.py',package='ISfit')
    source_python(aacgmfile)
    idate <- as.integer(date)
    coord <- runaacgmv2(inlat,inlon,height,idate[1],idate[2],idate[3],idate[4],idate[5],idate[6],methcode)

    return(coord)
}
