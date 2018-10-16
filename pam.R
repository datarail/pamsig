## Application of k-medoids to example CyCIF data
##
## by Artem Sokolov

library( readr )
suppressMessages(library( dplyr ))
library( stringr )
library( yaml )
library( purrr )
library( pheatmap )

## Clamps qq^th quantile to -1 and (1-qq)^th quantile to 1
quantnorm <- function( v, qq = 0.001 )
{
    stopifnot( qq < .5 )
    lo <- quantile( v, qq )
    hi <- quantile( v, 1 - qq )
    2 * (v - lo) / (hi - lo) - 1
}

## Given a data frame and cluster index computes t.test statistic for each
##   channel for that cluster vs. rest
tTestCluster <- function( iC, DF )
{
    S <- DF %>% split( .$Cluster == iC ) %>% map( ~select(.x, -Cluster) )
    map2_df( S[["TRUE"]], S[["FALSE"]], ~t.test(.x,.y)$statistic ) %>%
        mutate( Cluster = iC ) %>% select( Cluster, everything() )
}

main <- function()
{
    if( !dir.exists( "output" ) )
        stop( "No output/ directory exists. Please map the directory using docker -v" )
    
    ## Load the settings
    cat( "Reading settings from config/pamsig.yml\n" )
    sts <- read_yaml( "config/pamsig.yml" )
    
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

    ## cluster::pam() has a limit on the number of points it can handle
    ## Randomly subsample the data if above this limit
    pamLimit <- 65536
    if( nrow(X) > pamLimit )
    {
        cat( "cluster::pam() has a limit on the number of points it can process\n" )
        cat( "  Randomly subsampling to", pamLimit,"points to be within this limit\n" )
        i <- sample(1:nrow(X), pamLimit)
        X <- X[i,]
    }
    
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

    ## Compose the output filenames and save the outputs
    pfx <- tools::file_path_sans_ext(sts$data) %>% basename %>%
        str_c( "output/", ., "-pam", sts$k )
    fnCSV <- str_c( pfx, ".csv" )
    cat( "Storing results to", pfx, "files\n" )
    Res %>% write_csv( fnCSV )

    ## Write out a silhouette plot to a .pdf
    fnPDF <- str_c( pfx, ".pdf" )
    pdf( fnPDF )
    plot( sl )
    invisible(dev.off())

    ## Computes t.test statistic for each channel / cluster pair
    nC <- max( Res$Cluster )
    R <- Res %>% select( -Neighbor, -Silhouette ) %>% map( 1:nC, tTestCluster, . ) %>% bind_rows

    ## Store the results to a .csv file
    fnSigCSV <- str_c( pfx, "-sig.csv" )
    R %>% write_csv( fnSigCSV )

    ## Create a quick-and-dirty heatmap showing the signature
    M <- R %>% mutate( Cluster = str_c("Cluster", Cluster) ) %>% as.data.frame %>%
        tibble::column_to_rownames("Cluster") %>% as.matrix
    fnSigPDF <- str_c( pfx, "-sig.pdf" )
    pheatmap( M, cluster_rows=FALSE, filename=fnSigPDF, width=7, height=nC*0.25+2, silent=FALSE )
}

main()
