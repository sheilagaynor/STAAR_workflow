# Description: Generate a null model using the STAAR package which provides a wrapper of the GMMAT package.
# Inputs:
# pheno_file : file containing the outcome, covariates for the null model (.csv)
# null_file_name : string containing prefix for .Rds output from null model fitting via STAAR (string)
# sample_name : column name in pheno_file for observation IDs (string)
# outcome_name : column name in pheno_file for outcome (string)
# outcome_type : type of variable of outcome, outcome_name in pheno_file, 'continuous' or 'dichotomous' (string)
# covariate_names : column names in pheno_file of covariate variables, as comma-separated string, to be treated as covariates (string)
# kinship_file : file containing the kinship matrix for null model with relatedness, row names are sample_names (.Rds, .Rdata, .csv)
# het_var_name : column name in pheno_file of variable for grouping heteroscedastic errors (string)
# Outputs:
# Null model STAAR object saved as null_file.Rds

## Parse arguments
args <- commandArgs(T)

## Required arguments
# File inputs
pheno_file <- args[1]
null_file <- args[2]
# Analysis inputs
sample_name <- args[3]
outcome_name <- args[4]
outcome_type <- args[5] #default 'continuous'
# Optional inputs
covariate_names <- args[6] #default 'NA'
kinship_file <- args[7] #default 'NA'
het_var_name <- args[8] #default 'NA'

#####################
# Functions for input processing
# Adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_nullmodel.R 
get_family <- function(outcome_type_from_args) {
  if (toupper(outcome_type_from_args) == "CONTINUOUS"){
    family <- gaussian()
  } else if (toupper(outcome_type_from_args) == "DICHOTOMOUS"){
    family <- binomial()
  } else {
    stop_msg <- paste("Invalid outcome type provided:", outcome_type_from_args)
    stop(stop_msg)  }
  return(family)
}
get_kinship <- function(kinship_file_from_args){
  cat('Loading Kinship Matrix: ',kinship_file_from_args,'\n')
  if (grepl('Rda$',kinship_file_from_args,ignore.case=TRUE) | grepl('Rdata$',kinship_file_from_args,ignore.case=TRUE)){
    kins <- get(load(kinship_file_from_args))
  } else if (grepl('Rds$',kinship_file_from_args,ignore.case=TRUE)){
    kins <- readRDS(kinship_file_from_args) 
  } else if (grepl('csv$',kinship_file_from_args,ignore.case=TRUE)){
    kins <- as.matrix(read.csv(kinship_file_from_args,as.is=T,check.names=F,row.names=1))
  } else {
    stop_msg <- paste("Invalid kinship file provided:", kinship_file_from_args)
    stop(stop_msg)  }
  cat('Loaded Kinship: no. rows:',nrow(kins),' no. cols:',ncol(kins),'\n')
  kins
}

# Load packages
suppressMessages(library(STAAR))

#####################
# Read phenotypes, covariates
if (!grepl('csv$',pheno_file,ignore.case=TRUE)){
  stop_msg <- paste("Invalid phenotype file provided:", pheno_file)
  stop(stop_msg)
}
pheno <- read.csv(pheno_file, header=TRUE, as.is=TRUE)
pheno <- pheno[order(pheno[,sample_name]),]
cat('Loaded Phenotypes: no. rows:',nrow(pheno),' no. cols:',ncol(pheno),'\n')
# Subset to complete cases for phenotype file, create null model formula
if ( covariate_names=='NA' & het_var_name=='NA' ){
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name)]), c(sample_name,outcome_name)]
  null_model_char <- paste0(outcome_name, "~", 1)
} else if ( covariate_names=='NA' & het_var_name!='NA' ){
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,het_var_name)]), c(sample_name,outcome_name,het_var_name)]
  null_model_char <- paste0(outcome_name, "~", 1)
} else if ( covariate_names!='NA' & het_var_name=='NA' ){
  covar_split <- unlist(strsplit(covariate_names, split=","))
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,covar_split)]), c(sample_name,outcome_name,covar_split)]
  null_model_char <- paste0(outcome_name, "~", paste(covar_split, collapse="+"))
} else {
  covar_split <- unlist(strsplit(covariate_names, split=","))
  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,het_var_name,covar_split)]), 
                           c(sample_name,outcome_name,het_var_name,covar_split)]
  null_model_char <- paste0(outcome_name, "~", paste(covar_split, collapse="+"))
}
cat('Complete Phenotypes: no. rows:',nrow(phenotype_input),' no. cols:',ncol(phenotype_input),'\n')

#####################
# Get kinship
if ( kinship_file!='NA' ){
  kinship_input <- get_kinship(kinship_file)
  shared_obs <- intersect(row.names(kinship_input), phenotype_input[,sample_name])
  if (length(shared_obs)==0){
    stop('No shared observations between phenotypes and kinship')
  } else {
    shared_samples <- intersect(rownames(kinship_input), phenotype_input[,sample_name])
    kinship_analysis <- kinship_input[match(shared_samples,rownames(kinship_input)),match(shared_samples,rownames(kinship_input))]
    phenotype_analysis <- phenotype_input[match(shared_samples,phenotype_input[,sample_name]),] 
    rm(kinship_input); rm(phenotype_input)
    cat('Matched Phenotypes with Kinship: no. rows:',nrow(phenotype_analysis),' no. cols:',ncol(phenotype_analysis),'\n')
  }
} else {
  phenotype_analysis <- phenotype_input
  rm(phenotype_input)
}

#####################
# Fit, save null model
if ( kinship_file=='NA' & het_var_name=='NA' ){
  cat('Fitting null model for unrelated samples, homogeneous variance')
  null_model <- STAAR::fit_null_glm(as.formula(null_model_char), 
                data = phenotype_analysis, family = get_family(outcome_type))
} else if ( kinship_file=='NA' & het_var_name!='NA' ){
  cat('Fitting null model for unrelated samples, heterogeneous variance')
  null_model <- STAAR::fit_null_glm(as.formula(null_model_char), 
                groups = het_var_name, data = phenotype_analysis, family = get_family(outcome_type))
} else if ( kinship_file!='NA' & het_var_name=='NA' ){
  cat('Fitting null model for related samples, homogeneous variance')
  null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_name, 
                data = phenotype_analysis, kins = kinship_analysis, family = get_family(outcome_type))
} else {
  cat('Fitting null model for related samples, heterogeneous variance')
  null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_name, 
                groups = het_var_name, data = phenotype_analysis, kins = kinship_analysis, family = get_family(outcome_type))
}
saveRDS(null_model, file=paste0(null_file,".Rds"))
