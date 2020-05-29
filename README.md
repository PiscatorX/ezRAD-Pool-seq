# ezRAD-Pool-seq
## An ezRAD pool-seq data analysis pipeline ##

An implementation of the ezRAD pipeline described in Phair et. al. 2019.


** Downloading the _Z. marina_ genome from NCBI  **

* The [jq](https://stedolan.github.io/jq/) and [datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/command-line-start/) commandline tools were used to find and download the Z. marina genome
  ```
  datasets assembly-descriptors tax-name 'zostera marina' --limit ALL  | jq '.datasets[].assembly_accession' -r
  ```
* Check available types of annotation available for the assemble.
  ```
  

  ```



## Citation ##

Phair NL, Toonen RJ, Knapp I, von der Heyden S. 2019. Shared genomic outliers across two divergent population clusters of a highly threatened seagrass. PeerJ 7:e6806 http://doi.org/10.7717/peerj.6806
