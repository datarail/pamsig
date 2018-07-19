# pamsig: A simple k-medoids clustering method for CyCIF data

The repository hosts source code for Docker containers used for making cell state calls on segmented CyCIF data:

  * pamsig - Clusters segmented CyCIF data using Partitioning Around Medoids (PAM) and produces signatures that define each cluster. The config file expects the following fields:
    * k - The desired number of clusters
    * data - Filename of segmented CyCIF data that will appears in `input/` directory of the container
    * markers - List of markers that should be used for clustering
    
