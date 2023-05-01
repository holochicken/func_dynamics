# Metabolic capacity of the gut microbiota associates with performance of broiler chickens

Development of functional attributes of chicken-associated microbial communities, and their link to chicken body weight

## Repository structure
data:  
  - metadata
  - metagenomics
  - metatranscriptomics

src: R code
  - hmsc_mg_and_mt: steps for Hierarchical Modelling of Species Communities (Hmsc). These scripts need higher memory and are not recommended to run in local computers. 
    - mg: for metagenomic data
    - mt: for metatranscriptomic data

The script is designed to run mg and then mt. S1 scripts from both directories will create models, model_fit and panels directories necessary to save intermediate files.

The rest of the code in src folder is organised according to the outline of the article. Some of these scripts use Hmsc results. Tables and figures created with the code will be saved in results directory. The S1 script will create the required directories: results/tables and results/figures.
