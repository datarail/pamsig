FROM r-base

## Install relevant packages
RUN Rscript -e "install.packages(c('readr','dplyr','stringr','jsonlite'))"

## Copy the R script into the image
COPY pam.R /pam.R

## Execute the R script
CMD ["Rscript", "pam.R"]
