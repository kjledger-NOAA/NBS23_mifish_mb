# eDNA metabarcoding of 2023 NBS survey samples 

- MiFish 12S amplicon sequencing 
- demultiplexed sequences filtered using Dadasnake (config.NBS23.mifish.yaml)

1_taxonomic_assignment_blastn.Rmd
- performed taxonomic assignment of ASVs using NCBI nt database (query parameters: 96% seq identity and 98% query coverage)
- filtered taxonomic hits to retain only fish found in the NE Pacific (according to FishBase)
- based taxonomic assignment of ASVs on the top 0.5% of hits if the best match was >= 98% or on the top 1% of matches if the best match was < 98%

2_decontamination.Rmd
- performed additional "decontamination" of asv table 
  1) tag-jumping - subtracted the proportion of reads that jumped into the positive control samples from all samples
  2) removed ASVs w/o in-range fish taxonomic assignment
  3) removed ASVs that only showed up in controls and not in any environmental samples 
  4) removed replicates with very low total read count (<250)

3_summary.Rmd
- created taxon table by merging decontaminated ASV table and taxonomic assignments 

4_basis_catch.Rmd
- downloaded catch data on NBS 2023 survey from AKFIN database 

5_covariates.Rmd
- extracted some environmental covariates from eDNA sampling locations 
- data sources = afsc-gap-products/coldpool and alaska-groundfish-efh  


description of output files: 
* nbs23_taxon_reads.csv - includes eDNA metabarcoding read counts and sample metadata in long format (i.e. every unique PCR replicate/sample_ID and taxonomic assignment is a row)
extraction_ID: unique ID for eDNA extraction at ABL 
alternative_ID: original ID given to 1-L water sample (ex. NW2301_2_10_1: NW20301 = survey ID; 2 = station ID; 10 = depth; 1 = replicate ID)
collection_year: year of collection
collection_month: month of collection
collection_day: day of collection
location1: station (matches station/STATIONNUMBER in BASIS catch)
depth: depth of sample collection in meters
longitude: longitude of collection
latitude: latitude of collection
sample_ID: unique ID for each metabarcoding PCR replicate (each DNA extraction has three PCR replicates) 
replicate: replicate ID of each sample 
sample_type: designation of field sample = "sample"
taxon: taxonomic assignment 
taxonomic_level: rank or classification of taxonomic assignment (ie. species, geneus, etc.)
species to kingdom: full taxonomic classification of taxonomic assignment 
reads: number of eDNA metabarcoding reads  

* nbs23_basis_catch.csv
STATIONNUMBER is the same as location1 in "nbs23_taxon_reads.csv"

* NBS23_covariates.csv 
site: station_depth (ie.e 2_10 = station 2 @10m)
bdepth through sponge are environmental covariates used in the EFH walleye pollock SDMs: https://repository.library.noaa.gov/view/noaa/48659; see Fig 154 and Table 43 for full covariate names - covariates were extracted from single maps based on combined data from the full time series
bottomtemp_2023 and surfacetemp_2023:  temperatue data used to derive cold pool index -https://github.com/afsc-gap-products/coldpool 


