# STAAR_workflow: Rare variant analysis methods for WGS data
Maintainer: Sheila Gaynor
Version: 1.1

## Description:
Workflow to perform aggregate rare variant tests for sequencing studies and genetic data. Implements the variant-Set Test for Association using Annotation infoRmation (STAAR) procedure, as well as SKAT, Burden, and ACAT tests for both continuous and dichotomous traits. The STAAR method incorporates qualitative functional categories and quantitative complementary functional annotations (Li and Li et al, 2020). The workflow accounts for population structure and relatedness, and scales for large whole genome sequencing studies.

## Functionality:
The workflow contains three tasks, including two key analysis steps. The workflow fits a null model for testing, incorporating the outcome, covariates, and kinship (optional). The workflow then uses the null model, genotypes, and aggregation units (optional) to run rare variant association analyses. Lastly, results are compiled.

## Funcional inputs:
### Null model R/WDL inputs:
- **pheno_file**: [file] file containing the outcome, covariates for the null model (.csv)
- **null_file_name**: [string] string containing prefix for .Rds output from null model fitting via STAAR (string)
- **sample_name**: [string] column name in pheno_file for observation IDs (string)
- **outcome_name**: [string] column name in pheno_file for outcome (string)
- **outcome_type**: [string] type of variable of outcome, outcome_name in pheno_file, 'continuous' or 'dichotomous' (string)
- **covariate_names**: [string] optional, column names in pheno_file of covariate variables, as comma-separated string, to be treated as covariates (string)
- **kinship_file**: [file] optional, file containing the kinship matrix for null model with relatedness, row names are sample_names (.Rds, .Rdata, .csv)
- **het_var_name**: [string] optional, column name in pheno_file of variable for grouping heteroscedastic errors (string)
### Null model WDL inputs:
- **null_memory**: [int] requested memory in GB (numeric)
- **null_disk**: [int] requested disk size (numeric)

### Association test R/WDL inputs:
- **null_file**: [file] file containing output from null model fitting via STAAR (.Rds)
- **geno_file**: [file] file containing genotypes for all individuals from null model, optionally containing the given annotation channels (.gds)
- **annot_file**: [file] file containing annotations as input with columns 'chr', 'pos', 'ref', 'alt' (.Rds, .Rdata, .csv)
- **results_file**: [string] string of name of results file output (string)
- **agds_file**: [string] string indicating whether input geno is an agds file containing the annotations, 'None' [default] (string)
- **agds_annot_channels**: [string] comma-separated names of channels in agds to be treated as annotations, 'None' [default] (string)
- **agg_file**: [file] file containing the aggregation units in character strings for set-based analysis with columns 'chr', 'pos', 'ref', 'alt', 'group_id' (.Rds, .Rdata, .csv)
- **cond_file**: [file] file containing the variants to be conditioned upon with columns 'chr', 'pos', 'ref', 'alt' (.Rds, .Rdata, .csv)
- **cond_geno_files**: [file] file containing genotypes for all individuals from null model for conditional analysis; often same as geno_file (.gds)
- **cand_file**: [file] file containing units (agg_file required)/windows for candidate sets of interest with columns 'group_id' or 'chr', 'start', 'end' (.Rds, .Rdata, .csv)
- **maf_thres**: [int] AF threshold below which variants will be considered in rare variant analysis, 0.01 [default] (numeric)
- **mac_thres**: [int] AC threshold above which variants will be considered in rare variant analysis, 1 [default] (numeric)
- **window_length**: [int] length of window for region-based analysis, 2000 [default] (numeric)
- **step_length**: [int] length of overlap for region-based analysis, 1000 [default] (numeric)
- **num_cores**: [int] number of cores to be used in parallelized analysis, 3 [default] (numeric)
- **num_iterations**: [int] number of iterations to run in parallel loop, i.e. how many chunks to split sets into, 20 [default] (numeric)
### Association test WDL inputs:
- **null_memory**: [int] requested memory in GB (numeric)
- **null_disk**: [int] requested disk size (numeric)


## Resulting output:
The workflow produces a file containing the null model (.Rds) and all results of the aggregation test in a txt file (.txt).


