FROM r-base

## Install relevant packages
RUN Rscript -e "install.packages(c('readr','dplyr','stringr','jsonlite'))"

## Install additional packages
RUN Rscript -e "install.packages(c('purrr','pheatmap'))"

## Copy the R script into the image
COPY sig.R /sig.R

## Execute the R script
CMD ["Rscript", "sig.R"]
