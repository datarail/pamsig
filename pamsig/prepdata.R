## Prepares example dataset and settings files
##
## by Artem Sokolov

library( tidyverse )
library( jsonlite )

## Prepares the example dataset
## Original data taken from Ziming_2018 Benchmark
prepdata <- function()
{
    ## Load the original data
    X <- read_csv( "Ziming_LUNG3195.csv", col_types=cols() )

    ## Randomly subsample 5k points
    set.seed(100)
    X %>% slice( sample( 1:n(), 5000 ) ) %>% write_csv( "CyCIF1.csv" )

    ## Prepare default settings
    vRemove <- c( str_c("Hoechst", 0:9), "A488_0", "A488_1", "A555_0",
                 "A647_0", "Area", "Circ", "X", "Y", "frame" )
    vMarkers <- colnames(X) %>% setdiff(vRemove)
    sts <- list( k = 8, data = "input/CyCIF1.csv", markers = vMarkers )
    write_json( sts, "example-settings.json" )
}

