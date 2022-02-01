# STAAR_workflow: Rare variant analysis methods for WGS data
Maintainer: Sheila Gaynor
Version: 1.2

## Description:
Workflow to perform aggregate rare variant tests for sequencing studies and genetic data. Implements the variant-Set Test for Association using Annotation infoRmation (STAAR) procedure, as well as SKAT, Burden, and ACAT tests for both continuous and dichotomous traits. The STAAR method incorporates qualitative functional categories and quantitative complementary functional annotations (Li and Li et al, 2020). The workflow accounts for population structure and relatedness, and scales for large whole genome sequencing studies.

## Functionality:
The workflow contains three tasks, including two key analysis steps. The workflow fits a null model for testing, incorporating the outcome, covariates, and kinship (optional). The workflow then uses the null model, genotypes, and aggregation units (optional) to run rare variant association analyses. Lastly, results are compiled.

## Functional inputs:
### Null model R/WDL inputs - all user-provided with optional arguments denoted; note, this entire step will be skipped when providing a pre-computed null model (null_file_precompute argument):
- **pheno_file**: [file] file containing covariates and the outcome for the null model, with rows of samples and named columns of features (.csv)
- **null_file_name**: [string] string containing prefix for workflow-generated .Rds output from null model fitting via STAAR (string)
- **sample_name**: [string] column name in pheno_file for observation IDs (string)
- **outcome_name**: [string] column name in pheno_file for outcome (string)
- **outcome_type**: [string] type of variable of outcome in pheno_file, provided as 'continuous' [default] or 'dichotomous' (string)
- **covariate_names**: [string] optional, column names in pheno_file of covariate variables, provided as comma-separated string,  (string)
- **kinship_file**: [file] optional, file containing the kinship matrix for a null model with relatedness, where row names are observation IDs (.Rds, .Rdata, .csv)
- **het_var_name**: [string] optional, column name in pheno_file of variable for grouping heteroscedastic errors (string)
### Null model WDL inputs:
- **null_memory**: [int] optional, requested memory in GB (numeric)
- **null_disk**: [int] optional, requested disk size (numeric)

### Association test R/WDL inputs - all user-provided with optional arguments denoted and annotations provided in **either** an external file (annot_file) or within the gds file as annotation channels (agds_file_type, agds_annot_channels); note, the null model generated in the previous task is provided directly to the association task when the optional null_file_precompute argument is not provided:
- **null_file_precompute**: [file] optional, user-provided or work-flow generated, file containing output from null model fitting via STAAR (.Rds)
- **geno_file**: [file] file(s) containing genotypes for all individuals from null model, optionally containing the given annotation channels if in agds format (.gds)
- **annot_file**: [file] optional, file containing annotations as input with columns 'chr', 'pos', 'ref', 'alt' and named annotation columns for rows of variants (.Rds, .Rdata, .csv)
- **results_file_name**: [string] string of name of results file output (string)
- **agds_file_input**: [string] optional, string indicating whether input geno_file is an agds file containing the annotations, 'None' [default] or any other value (string)
- **agds_annot_channels**: [string] optional, comma-separated names of channels in agds to be treated as annotations, 'None' [default] or comma-separated names (string)
- **agg_file**: [file] optional, file containing the aggregation units in character strings for set-based analysis with columns 'chr', 'pos', 'ref', 'alt', 'group_id' for rows of variants (.Rds, .Rdata, .csv)
- **cond_file**: [file] optional, file containing the variants to be conditioned upon with columns 'chr', 'pos', 'ref', 'alt' for rows of variants (.Rds, .Rdata, .csv)
- **cond_geno_files**: [file] optional, file containing genotypes for all individuals from null model for conditional analysis; often same as geno_file (.gds)
- **cand_file**: [file] optional, file containing units (agg_file required) or windows for candidate sets of interest with columns 'group_id' or 'chr', 'start', 'end' (.Rds, .Rdata, .csv)
- **maf_thres**: [int] optional, AF threshold below which variants will be considered in rare variant analysis, 0.01 [default] (numeric)
- **mac_thres**: [int] optional, AC threshold above which variants will be considered in rare variant analysis, 1 [default] (numeric)
- **window_length**: [int] optional, length of window for region-based analysis, 2000 [default] (numeric)
- **step_length**: [int] optional, length of overlap for region-based analysis, 1000 [default] (numeric)
- **num_cores**: [int] optional, number of cores to be used in parallelized analysis, 3 [default] (numeric)
- **num_iterations**: [int] optional, number of iterations to run in parallel loop, i.e. how many chunks to split sets into, 20 [default] (numeric)
### Association test WDL inputs:
- **test_memory**: [int] optional, requested memory in GB (numeric)
- **test_disk**: [int] optional, requested disk size (numeric)


## Resulting output:
The workflow produces a file containing the null model (.Rds) and all results of the aggregation test in a txt file (.txt).


## Key notes:
- Example input files are available in this respository's testfiles/ folder
- The workflow will use the intersection of available data, for both aligning the kinship and phenotype file and aligning the annotation and genotype files, and requires the correct and available string names to be provided
- The workflow is designed to test multiple sets
- The workflow runs optimally when the annotations are provided in an integrated aGDS file, which takes precedence for provided annotations
- Memory use is related to both the number of cores (how many parallel processes run; smaller values require less memory with slower run time) and number of iterations (how many chunks of tests considered within a loop; smaller values require more memory with faster run time)
- Functional annotation data can be retrieved from many resources, including the [Functional Annotation of Variant - Online Resource (FAVOR)](http://favor.genohub.org) site
- A CWL version of the workflow has been prepared by Lea Ackovic at Seven Bridges and is included in the CWL subdirectory


## Citation:
For more details and citation, please visit https://doi.org/10.1101/2021.09.07.456116.
