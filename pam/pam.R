## Application of k-medoids to example CyCIF data
##
## by Artem Sokolov

library( readr )
suppressMessages(library( dplyr ))
library( stringr )
library( jsonlite )

## Clamps qq^th quantile to -1 and (1-qq)^th quantile to 1
quantnorm <- function( v, qq = 0.001 )
{
    stopifnot( qq < .5 )
    lo <- quantile( v, qq )
    hi <- quantile( v, 1 - qq )
    2 * (v - lo) / (hi - lo) - 1
}

main <- function()
{
    if( !dir.exists( "output" ) )
        stop( "No output/ directory exists. Please map the directory using docker -v" )
    
    ## Load the settings
    cat( "Reading settings from input/settings.json\n" )
    sts <- read_json( "input/settings.json", simplifyVector=TRUE )
    
    ## Load the data and reduce to the requested set of markers
    fnIn <- str_c("input/", sts$data)
    cat( "Loading CyCIF data from", fnIn, "\n" )
    Xraw <- read_csv( fnIn, col_types=cols() )
    cat( "  Parsed", nrow(Xraw), "points across", ncol(Xraw), "columns\n" )

    ## Check that all the markers are present
    vMissing <- setdiff( sts$markers, colnames(Xraw) )
    if( length(vMissing) > 0 )
    {
        msg <- str_c( "The following markers are missing: ",
                     str_c( vMissing, collapse=", ") )
        stop( msg )
    }

    ## Reduce to the requested set of markers
    Xsub <- Xraw %>% select( one_of( sts$markers ) )
    cat( "  Reduced data to the", ncol(Xsub), "requested markers\n" )

    ## Preprocess the data
    ##   a) log10-transform
    ##   b) quantile normalization of each channel
    cat( "Applying log transform and quantile normalization...\n" )
    X <- Xsub %>% mutate_all( ~log10(. + 1) ) %>% mutate_all( ~quantnorm(.x) )
    
    ## Apply k-medoids
    cat( "Applying k-medoids with k =", sts$k, "...\n" )
    cls <- cluster::pam( X, sts$k )

    ## Compute silhouette scores
    sl <- cluster::silhouette( cls )
    SL <- data_frame( rowIndex = as.integer(rownames(sl)),
                     Cluster = sl[,1], Neighbor = sl[,2], Silhouette = sl[,3] )

    ## Compose the final result
    Res <- X %>% mutate( rowIndex = 1:n() ) %>% inner_join( SL, by="rowIndex" ) %>%
        select(-rowIndex) %>% select( Cluster:Silhouette, everything() )
    cat( "Average Silhouette Score:", mean(Res$Silhouette), "\n" )

    ## Compose the output filename and save the output
    fnOut <- tools::file_path_sans_ext(sts$data) %>% basename %>%
        str_c( "output/", ., "-pam", sts$k, ".csv" )
    cat( "Storing results to", fnOut, "\n" )
    Res %>% write_csv( fnOut )
}

main()
