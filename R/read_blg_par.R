#' Reading Bilog Parameter File
#'
#' @param infile Bilog .par file name
#'
#' @export
#'

read_blg_par <- function( infile ){
 # read bilog .par file
 # Shin-ichi Mayekawa
 # 20230721
 # 
 
 res=read.table( infile, skip=4, fill=TRUE, stringsAsFactors=0
                , col.names=c("name","subtest","intercept","intercept_se"
                              , "a", "a_se","b","b_se","disp","disp_se"
                              , "c", "c_se", "DRIFT", "DRIFT_se"
                              , "location"))
 return( res )
} # end of read_blg_par

