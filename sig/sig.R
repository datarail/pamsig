## Given a clustering solution, computes signatures that define each cluster
##
## by Artem Sokolov

library( readr )
suppressMessages(library( dplyr ))
library( stringr )

library( purrr )
library( pheatmap )

## Given a data frame and cluster index computes t.test statistic for each
##   channel for that cluster vs. rest
tTestCluster <- function( DF, iC )
{
    cat( "Computing signature for cluster", iC, "...\n" )
    S <- DF %>% split( .$Cluster == iC ) %>% map( ~select(.x, -Cluster) )
    map2_df( S[["TRUE"]], S[["FALSE"]], ~t.test(.x,.y)$statistic ) %>%
        mutate( Cluster = iC ) %>% select( Cluster, everything() )
}

main <- function()
{
    if( !dir.exists( "output" ) )
        stop( "No output/ directory exists. Please map the directory using docker -v" )
    
    ## Scan for the clustering file
    fns <- list.files( "input", "*.csv" )
    if( length(fns) != 1 )
        stop( "Expecting a single .csv file in input/ directory" )
    
    ## Load the provided file
    fn <- str_c( "input/", fns )
    cat( "Loading clustering solution from", fn, "...\n" )
    X <- read_csv( fn, col_types=cols() )

    ## Drop Neighbor and Silhouette columns, if present
    if( "Silhouette" %in% colnames(X) )
        X <- X %>% select( -Neighbor, -Silhouette )

    ## Display basic statistics about the data
    nC <- max( X$Cluster )
    cat( "Found assignments for", nC, "clusters\n" )
    
    ## Computes t.test statistic for each channel / cluster pair
    R <- map( 1:nC, ~tTestCluster( X, .x ) ) %>% bind_rows

    ## Store the results to a .csv file
    pfx <- tools::file_path_sans_ext(fn) %>% basename %>% str_c("output/", ., "-sig")
    fnCSV <- str_c( pfx, ".csv" )
    R %>% write_csv( fnCSV )

    ## Create a quick-and-dirty heatmap showing the signature
    M <- R %>% mutate( Cluster = str_c("Cluster", Cluster) ) %>% as.data.frame %>%
        tibble::column_to_rownames("Cluster") %>% as.matrix
    fnPDF <- str_c( pfx, ".pdf" )
    pheatmap( M, cluster_rows=FALSE, filename=fnPDF, width=7, height=nC*0.25+2, silent=FALSE )
}

main()
