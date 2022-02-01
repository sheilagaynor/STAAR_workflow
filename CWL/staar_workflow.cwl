{
  "class": "Workflow",
  "cwlVersion": "v1.2",
  "doc": "**STAAR Rare Variant Workflow** performs rare variant association testing for sequencing studies and genetic data. The workflow Implements the variant-Set Test for Association using Annotation infoRmation (STAAR) procedure, as well as SKAT, Burden, and ACAT tests for both continuous and dichotomous traits [1]. The STAAR method incorporates qualitative functional categories and quantitative complementary functional annotations [2]. The workflow accounts for population structure and relatedness, and scales for large whole genome sequencing studies [1].\n\n**STAAR Rare Variant Workflow** is designed to first fit the null model for testing, incorporating the outcome, covariates, and kinship (optional). The workflow then uses the null model, genotypes, and aggregation units (optional) to run rare variant association analyses [1]. At the end, all the results from association testing are compiled together.\n\n**STAAR Rare Variant Workflow** fits a regression or generalized linear mixed model under the null hypothesis for related or unrelated samples. When provided with **Kinship file**, the tool will fit generalized linear mixed model for related samples. Depending on the nature of the outcome, the null model can be fitted for continuous (Gaussian) and dichotomous (binomial) **Outcome type**. When samples come from different groups (e.g., study or ancestry group), it is common to observe different variances by group for quantitative traits. In that case, it is recommended to fit null model using the parameter **Heteroscedastic variable name**.\n\nAfter fitting the null model, **STAAR Rare Variant Workflow** executes association testing. The workflow can be run in two association testing modes: ***Genomewide*** and ***Candidate***. In ***Genomewide*** mode the association testing is performed on a whole genome, while in ***Candidate*** mode the testing is restricted to the variants subset provided in **Candidate file**. ***Candidate*** mode is  triggered when **Candidate file** is provided on the input, otherwise the workflow will be run in ***Genomewide*** mode. When provided with **Aggregate file** the workflow will do the aggregate testing on the variants present in the file. The workflow also offers the conditional association testing. The conditional analysis is performed when the **Conditional file** is provided on the input.\n\nTo run the workflow following input files are required: **Phenotype file** and **Genotype file**. **STAAR Rare Variant Workflow** outputs fitted **Null model** in Rds format, association testing **STAAR results** in compressed gz format for each provided **Genotype file** and **Compiled STAAR results** with all the **STAAR results** compiled into one file. \n\nThe workflow utilises **STAAR_workflow** [1] implementation of **STAAR 0.9.5** [2].\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n### Common use cases\n\n**STAAR Rare Variant Workflow** can be used for performing both ***genomewide*** or ***candidate*** rare variant association testing. By default ***genomewide*** testing is used, and if **Candidate file** is provided it will run testing in ***candidate*** mode. Also, it is possible to run conditional association testing, which is triggered if **Conditional file** is provided on the input. The workflow will fit the null model first, using **Phenotype file** and **Kinship file** (optional), and then run specified association analysis using fitted null model, **Genotype file** and other optional inputs (**Annotation file**, **Aggregate file** and **Candidate  file**). Lastly,  the workflow will compile all the results from association testing. The workflow outputs fitted **Null model** in Rds format, association testing **STAAR results** in compressed gz format and **Compiled STAAR results** in txt format.\n\n### Changes Introduced by Seven Bridges\n- This workflow has the same implementation and set up as in the [STAAR_workflow WDL code](https://github.com/sheilagaynor/STAAR_workflow/blob/main/STAAR_analysis_workflow.wdl) with some small modifications introduced. \n- Beside the **Compiled STAAR results**, the **STAAR results** from single/parallel executions are also provided on the output.\n- In the original WDL implementation **STAAR Association Testing** step is scattered both by **Genotype file** and **Annotation file** input. In SBG implementation this step is only scattered by **Genotype file** and all parallel jobs use the same **Annotation file** as the script will take only annotations for intersected variants present both in **Genotype file** and **Annotation file**. In this case, it is required that **Annotation file** contains annotations for all provided **Genotype files**. This way the user doesn't have to pay attention to the order in which **Genotype file** and **Annotation file** are provided. \n- BASH script for results compilation was slightly adapted for platform execution, but it creates the same output file as in the WDL implementation.\n\n### Common issues and important notes\n- **Kinship file** must be provided in Rds, RData or csv format.\n- **Outcome type** must reflect the nature of outcome: it can be ***continuous*** (Gaussian) or ***dichotomous*** (binomial). If not defined, the outcome type will be set to ***continuous***.\n- Beside the **Outcome type**, the other required inputs are: **Null file name** - the file name to be assigned to the null model, **Column name of observation/id variable** - name of a column with observations or ids, **Column name of outcome variable** and **Results filename**.\n- When **Kinship file** and **Phenotype file** are provided on the input, **STAAR Null Model** processes data in following way: the tool will load both files, find the intersection between them based on sample ids and create a null model for intersected samples using defined input parameters for outcome variable, covariates and outcome type.\n- **Genotype file** input is required and has to be in gds format. If **Annotation channels** and **Annotated GDS filename** are provided, and the **Annotation file** is not present on the input, the tool will assume that annotations are present in annotated gds **Genotype file**. More information about annotated gds format can be found on the [STAAR workflow GitHub page](https://github.com/sheilagaynor/STAAR_workflow).\n- If the **Annotated GDS filename** has any other value than `None`, the tool will assume there are annotations present in the **Genotype file**.\n- If the **Annotation file** is present on the input, the tool will take those annotations and ignore the ones in the GDS file.\n- When the **Annotation file** is provided on the input, the **STAAR Association Testing** step will analyse only the variants that are present both in the **Annotation file** and **Genotype file**.\n- When multiple **Genotype files** are provided on the input together with the **Annotation file**, the **Annotation file** should contain annotations for the variants from all genotype files, otherwise it will only analyse the ones with present annotation.\n- **STAAR Association Testing** step is scattered into multiple parallel jobs when multiple **Genotype files** are provided on the input and will output resulting file for each **Genotype file**.\n- **Aggregate file** must be provided in .Rds, .Rdata or .csv format with following columns:\n\n```\n\"chr\",\"pos\",\"ref\",\"alt\",\"group_id\"\n\"1\",100000,\"T\",\"C\",\"EXAMPLE\"\n```\n- When running analysis with an **Aggregate file**, **Number of genome chunks** is the number of units to consider at a time within a parallel loop, so it should be set accordingly.\n- To run conditional analysis, **Conditional file** must be provided in .Rds, .Rdata or .csv format. Example format:\n```\n\"chr\",\"pos\",\"ref\",\"alt\"\n\"1\",100000,\"T\",\"C\"\n```\n- Other example input files can be found in the [test files directory of STAAR workflow GitHub page](https://github.com/sheilagaynor/STAAR_workflow/tree/main/testfiles).\n- **Number of cores** is a total number of cores that will be used for multi-threading. The computing instance will be chosen for the task based on the **Number of cores** value. When **Number of cores** is not provided, the tool will check if **CPUs per job** was defined. In case that both **Number of cores** and **CPUs per job** are not provided by the user, the tool will by default use the instance with 8 cores and set **Number of cores** to 8 to be used for multi-threading.\n- **Number of iterations** input defines number of iterations to run in parallel loop, i.e. how many chunks to split sets into [1]. \n- To run candidate analysis, **Candidate file** must be provided on the input and it is expected to have following columns: *group_id* or *chr* (chromosome), *start* (start position) and *end* (end position). When only *group_id* is provided in candidate file, it is required to also provide **Aggregate file** on the input.\n- **STAAR Association Testing** has quite high RAM consumption. RAM consumption is directly related to two input parameters: **Number of cores** and **Number of iterations**. **Number of cores** determines how many parallel processes will be run and **Number of iterations** determines how many chunks will be considered at a time within a parallel loop. These two inputs should be set up in accordance with the instance the task is run on. For **Number of iterations**, smaller values will require more RAM, but the analysis will finish faster, and higher values will require less RAM, but will prolong the duration. The opposite applies to **Number of cores**, higher number leads to faster execution and higher RAM requirement.\n- When running analysis in *genomewide* mode it is recommended to choose a custom instance in the execution settings tab in accordance with values for **Number of iterations** or **Number of cores**. Based on our tests, an r5.16xlarge instance proved to be most optimal (with **Number of cores** set to 32 and **Number of iterations** set to 20 for genotype file with 1 million variants, with 50000 variants analysed per iteration).\n- The tool can sometimes fail with memory error - despite the fact that instance metrics report more available memory. If the task fails with error \"core dumped\", \"segmentation fault\" or any memory related error, user should decrease **Number of cores** or increase **Number of iterations**. The other solution would be to rerun the task on an instance with more RAM.\n- For bigger files it is recommended to run *genomewide* analysis tasks on on-demand instances (i.e., disable spot instances in \"Settings\" for the project and/or \"App Settings\" for the task). This is because tasks with bigger files are long-running and take place on large/expensive instances, increasing the chance that the instance is pre-empted. At this point, the tool will by default switch to an on-demand instance and will start over the analysis (as the intermediary data is not saved).\n- As first and last tool in the workflow have considerably shorter execution times it is recommended to run the whole workflow on the one custom instance, so there are no changes of instances and file transfers during the task execution.\n\n### Performance Benchmarking\n\nPerformance benchmark was run with **Phenotype file** containing 25000 samples and Genotype files for Chr20 and Chr10 containing 22 and a half and 47 millions variants, respectively. For all the tasks **Number of cores** was set to 64 and **Number of iterations** to 220 and 470, respectively.\n      \n| Experiment type  | Duration | Cost (Instances + Attached Disks) | Instance (GCP on-demand)|\n|---------------------------|------------------------|-----------------------|--------------------------------|\n| Chr20 | 3h 41min | $11.39 ($11.21 + $0.19) | n1-standard-64 - 200GB disk |\n| Chr10 | 11h 23min | $43.71 ($43.14 + $0.57) | n1-highmem-64 - 200GB disk |\n\n\n*Cost can be significantly reduced by spot instance usage. Visit the [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### Portability\n\n**STAAR Rare Variant Workflow** was tested with cwltool version 3.1.20210816212154. The inputs `in_geno_file`, `in_pheno_file`, `in_agg_file`, `candidate_file`, `results_filename`, `agds_filename`, `agds_annot_chanels`, `null_file_name`, `sample_id_column_name`, `outcome_var_column_name` and `covariate_name` were provided in the job.yaml/job.json file and used for testing (linked [files](https://github.com/sheilagaynor/STAAR_workflow/tree/main/testfiles) and [input parameters](https://github.com/sheilagaynor/STAAR_workflow/blob/main/test.json) were used). All other input parameters are set to the same default values as in the WDL implementation inside the workflow and can't be modified by the user on the input.\n\n\n### References\n[1] [STAAR_workflow](https://github.com/sheilagaynor/STAAR_workflow)\\\n[2] [STAAR Github](https://github.com/xihaoli/STAAR)\\\n[3] [STAAR manual](https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf)",
  "label": "STAAR Rare Variant Workflow",
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "inputs": [
    {
      "id": "in_pheno_file",
      "sbg:fileTypes": "CSV",
      "type": "File",
      "label": "Phenotype file",
      "doc": "Phenotype file to be used for fitting null model.",
      "sbg:x": -555.3901977539062,
      "sbg:y": -201.5
    },
    {
      "id": "in_geno_file",
      "sbg:fileTypes": "GDS",
      "type": "File[]",
      "label": "Genotype file",
      "doc": "Genotype file, annotated GDS file containing the given annotation channels.",
      "sbg:x": -559,
      "sbg:y": 46
    },
    {
      "id": "in_kinship_file",
      "sbg:fileTypes": "RDS, RDATA, CSV",
      "type": "File?",
      "label": "Kinship file",
      "doc": "File containing the kinship matrix for related null model (.Rds, .Rdata or .csv).",
      "sbg:x": -558.3901977539062,
      "sbg:y": -75.5
    },
    {
      "id": "in_annot_file",
      "sbg:fileTypes": "RDATA, RDS, RDA, CSV",
      "type": "File?",
      "label": "Annotation file",
      "doc": "File containing annotation data.",
      "sbg:x": -546.8465576171875,
      "sbg:y": 393.25048828125
    },
    {
      "id": "in_agg_file",
      "sbg:fileTypes": "RDS, RDA, CSV, RDATA",
      "type": "File?",
      "label": "Aggregate file",
      "doc": "File variant masks, with a column for gene name and a column for unique variant identifier.",
      "sbg:x": -542.0873413085938,
      "sbg:y": 505.4912414550781
    },
    {
      "id": "candidate_file",
      "sbg:fileTypes": "RDATA, RDS, RDA, CSV",
      "type": "File?",
      "label": "Candidate file",
      "doc": "Candidate set file. This input is required for running the tool in Candidate mode.",
      "sbg:x": -538.9049072265625,
      "sbg:y": 627.4691162109375
    },
    {
      "id": "in_cond_file",
      "sbg:fileTypes": "RDA, RDS, CSV",
      "type": "File?",
      "label": "Conditional file",
      "doc": "Conditional file.",
      "sbg:x": -550.0349731445312,
      "sbg:y": 278.9359130859375
    },
    {
      "id": "in_cond_geno_files",
      "sbg:fileTypes": "GDS",
      "type": "File[]?",
      "label": "Conditional genotype files",
      "doc": "Conditional genotype files.",
      "sbg:x": -552.6834716796875,
      "sbg:y": 162.77281188964844
    },
    {
      "id": "null_file_name",
      "type": "string",
      "label": "Null file name",
      "doc": "String for naming the null model file. The name should not include file extension.",
      "sbg:exposed": true
    },
    {
      "id": "sample_id_column_name",
      "type": "string",
      "label": "Column name of observation/id variable",
      "doc": "Column name of observation/id variable.",
      "sbg:exposed": true
    },
    {
      "id": "outcome_var_column_name",
      "type": "string",
      "label": "Column name of outcome variable",
      "doc": "Column name of outcome variable.",
      "sbg:exposed": true
    },
    {
      "id": "outcome_type",
      "type": [
        "null",
        {
          "type": "enum",
          "symbols": [
            "continuous",
            "dichotomous"
          ],
          "name": "outcome_type"
        }
      ],
      "label": "Outcome type",
      "doc": "Type of outcome variable in phenofile. It can be 'continuous' or 'dichotomous'. The default is 'continuous'.",
      "sbg:toolDefaultValue": "continuous",
      "sbg:exposed": true
    },
    {
      "id": "covariate_name",
      "type": "string?",
      "label": "Covariate variables in phenotype file",
      "doc": "Comma-separated names of covariate variables in phenotype file to be treated as covariates.",
      "sbg:toolDefaultValue": "NA",
      "sbg:exposed": true
    },
    {
      "id": "het_var_name",
      "type": "string?",
      "label": "Heteroscedastic variable name",
      "doc": "Column name in phenotype file or group for heteroscedastic errors.",
      "sbg:toolDefaultValue": "NA",
      "sbg:exposed": true
    },
    {
      "id": "results_filename",
      "type": "string",
      "label": "Results filename",
      "doc": "Filename for results file.",
      "sbg:exposed": true
    },
    {
      "id": "agds_filename",
      "type": "string?",
      "label": "Annotated GDS filename",
      "doc": "Filename for annotated GDS file indicating whether input geno is an agds file containing the annotations.",
      "sbg:toolDefaultValue": "None",
      "sbg:exposed": true
    },
    {
      "id": "agds_annot_channels",
      "type": "string?",
      "label": "Annotation channels",
      "doc": "Comma-separated names of annotation channels in annotation GDS file.",
      "sbg:toolDefaultValue": "None",
      "sbg:exposed": true
    },
    {
      "id": "maf_thres",
      "type": "float?",
      "label": "Minor allele frequency threshold",
      "doc": "Minor allele frequency threshold. AF threshold below which variants will be considered in rare variant analysis. The default is 0.01.",
      "sbg:toolDefaultValue": "0.01",
      "sbg:exposed": true
    },
    {
      "id": "mac_thres",
      "type": "float?",
      "label": "Minor allele count threshold",
      "doc": "Minor allele count threshold. AC threshold above which variants will be considered in rare variant analysis. The default is 1.",
      "sbg:toolDefaultValue": "1",
      "sbg:exposed": true
    },
    {
      "id": "window_length",
      "type": "int?",
      "label": "Window length",
      "doc": "Length of window for region-based analysis. The default value is 2000.",
      "sbg:toolDefaultValue": "2000",
      "sbg:exposed": true
    },
    {
      "id": "step_length",
      "type": "int?",
      "label": "Step length",
      "doc": "Length of overlap for region-based analysis. The default value is 1000.",
      "sbg:toolDefaultValue": "1000",
      "sbg:exposed": true
    },
    {
      "id": "num_iterations",
      "type": "int?",
      "label": "Number of iterations",
      "doc": "Number of iterations to run in parallel loop, i.e. how many chunks to split sets into. Default is 20.",
      "sbg:toolDefaultValue": "20",
      "sbg:exposed": true
    },
    {
      "id": "num_cores",
      "type": "int?",
      "label": "Number of cores",
      "doc": "Total number of cores number of cores to be used in parallelized analysis. The default value is 8.",
      "sbg:toolDefaultValue": "8",
      "sbg:exposed": true
    },
    {
      "id": "mem_per_job",
      "type": "int?",
      "label": "Memory per job [MB]",
      "doc": "Memory per job [MB] to be used in task run. The default value is 16000.",
      "sbg:toolDefaultValue": "16000",
      "sbg:exposed": true
    }
  ],
  "outputs": [
    {
      "id": "null_model",
      "outputSource": [
        "staar_fit_null_model_0_9_5/null_model"
      ],
      "sbg:fileTypes": "Rds",
      "type": "File?",
      "label": "Null model",
      "doc": "Fitted null model.",
      "sbg:x": 450.1925354003906,
      "sbg:y": -113.1630630493164
    },
    {
      "id": "assoc_test_results",
      "outputSource": [
        "staar_association_test_0_9_5/assoc_test_results"
      ],
      "sbg:fileTypes": "GZ",
      "type": "File[]",
      "label": "STAAR results",
      "doc": "STAAR association test results",
      "sbg:x": 450.3555908203125,
      "sbg:y": 76.8369369506836
    },
    {
      "id": "compiled_results",
      "outputSource": [
        "staar_result_compilation_0_9_5/compiled_results"
      ],
      "sbg:fileTypes": "TXT",
      "type": "File?",
      "label": "Compiled STAAR results",
      "doc": "Compiled STAAR results from multiple STAAR results files.",
      "sbg:x": 466.3261413574219,
      "sbg:y": 310
    }
  ],
  "steps": [
    {
      "id": "staar_fit_null_model_0_9_5",
      "in": [
        {
          "id": "in_pheno_file",
          "source": "in_pheno_file"
        },
        {
          "id": "null_file_name",
          "source": "null_file_name"
        },
        {
          "id": "sample_id_column_name",
          "source": "sample_id_column_name"
        },
        {
          "id": "outcome_var_column_name",
          "source": "outcome_var_column_name"
        },
        {
          "id": "outcome_type",
          "default": "continuous",
          "source": "outcome_type"
        },
        {
          "id": "covariate_name",
          "source": "covariate_name"
        },
        {
          "id": "in_kinship_file",
          "source": "in_kinship_file"
        },
        {
          "id": "het_var_name",
          "source": "het_var_name"
        }
      ],
      "out": [
        {
          "id": "null_model"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-fit-null-model-0-9-5/22",
        "baseCommand": [
          "Rscript",
          "./STAAR_null_model.R"
        ],
        "inputs": [
          {
            "id": "in_pheno_file",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 0
            },
            "label": "Phenotype file",
            "doc": "Phenotype file to be used for fitting null model.",
            "sbg:fileTypes": "CSV"
          },
          {
            "id": "null_file_name",
            "type": "string",
            "inputBinding": {
              "shellQuote": false,
              "position": 1
            },
            "label": "Null file name",
            "doc": "String for naming the null model file. The name should not include file extension."
          },
          {
            "id": "sample_id_column_name",
            "type": "string",
            "inputBinding": {
              "shellQuote": false,
              "position": 2
            },
            "label": "Column name of observation/id variable",
            "doc": "Column name of observation/id variable."
          },
          {
            "id": "outcome_var_column_name",
            "type": "string",
            "inputBinding": {
              "shellQuote": false,
              "position": 3
            },
            "label": "Column name of outcome variable",
            "doc": "Column name of outcome variable."
          },
          {
            "sbg:toolDefaultValue": "continuous",
            "sbg:category": "Config Inputs",
            "id": "outcome_type",
            "type": [
              "null",
              {
                "type": "enum",
                "symbols": [
                  "continuous",
                  "dichotomous"
                ],
                "name": "outcome_type"
              }
            ],
            "inputBinding": {
              "shellQuote": false,
              "position": 4
            },
            "label": "Outcome type",
            "doc": "Type of outcome variable in phenofile. It can be 'continuous' or 'dichotomous'. The default is 'continuous'.",
            "default": "continuous"
          },
          {
            "sbg:toolDefaultValue": "NA",
            "id": "covariate_name",
            "type": "string?",
            "inputBinding": {
              "shellQuote": false,
              "position": 5
            },
            "label": "Covariate variables in phenotype file",
            "doc": "Comma-separated names of covariate variables in phenotype file to be treated as covariates.",
            "default": "NA"
          },
          {
            "id": "in_kinship_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 6
            },
            "label": "Kinship file",
            "doc": "File containing the kinship matrix for null model with relatedness, row names are sample_names. (.Rds, .Rdata or .csv).",
            "sbg:fileTypes": "RDS, RDA, CSV, RDATA"
          },
          {
            "sbg:toolDefaultValue": "NA",
            "id": "het_var_name",
            "type": "string?",
            "inputBinding": {
              "shellQuote": false,
              "position": 7,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 'NA'\n    }\n}"
            },
            "label": "Heteroscedastic variable name",
            "doc": "Column name in phenotype file or group for heteroscedastic errors.",
            "default": "NA"
          },
          {
            "sbg:toolDefaultValue": "16000",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job [MB]",
            "doc": "Memory per job [MB] to be used in task execution. The default value is 16000."
          },
          {
            "sbg:toolDefaultValue": "8",
            "id": "cpu_per_job",
            "type": "int?",
            "label": "CPUs per job",
            "doc": "Number of CPUs to be used per job. The default value is 8."
          }
        ],
        "outputs": [
          {
            "id": "null_model",
            "doc": "Fitted null model.",
            "label": "Null model",
            "type": "File?",
            "outputBinding": {
              "glob": "${\n    var name=inputs.null_file_name\n    var ext=name.split('.').pop()\n    if (ext != '.Rds'){\n        return name + '.Rds'\n    }\n    else {\n        return name\n    }\n}",
              "outputEval": "$(inheritMetadata(self, inputs.in_pheno_file))"
            },
            "sbg:fileTypes": "Rds"
          }
        ],
        "doc": "**STAAR Null Model** fits a regression or generalized linear mixed model under the null hypothesis. The model can be fitted for related or unrelated samples [1]. The tool fits a null model which can be used for association testing, incorporating the outcome, covariates, and kinship (optional) [2].\n\nContinuous (Gaussian) and dichotomous (binomial) outcomes are both supported (set in parameter **Outcome type**). Covariates from the **Phenotype file** are specified in parameter **Covariate variables in phenotype file**. When **Kinship file** is provided, the tool will fit generalised linear mixed model for related samples. \n\n**Heteroscedastic variable name** is an optional categorical variable indicating the groups which should be used for fitting a heteroscedastic linear mixed model (allowing residual variances in different groups to be different) [3]. When samples come from different groups (e.g., study or ancestry group), it is common to observe different variances by group for quantitative traits. In that case, it is recommended to fit null model using the parameter **Heteroscedastic variable name**. **Phenotype file** is required to run the tool.\n\nThe tool utilises **STAAR Null Model** implementation from **STAAR_Rare_Variant_Pipeline** [2].\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n### Common use cases:\n\nFitting a null model is a required step before doing the association testing. The **STAAR Null Model** will fit and output a null model in .Rds format, which can be used as an input for **STAAR association testing**. Null model file contains adjusted outcome values, the model matrix, the estimated covariance structure, elements indicating if samples are related or unrelated and additional parameters needed for association testing. \n\n### Common issues and important notes:\n- **Kinship file** must be provided in Rds, RData or csv format.\n- **Outcome type** must reflect the nature of outcome: it can be ***continuous*** (Gaussian) or ***dichotomous*** (binomial).\n- Beside the **Outcome type**, the other required inputs are: **Null file name** - the file name to be assigned to the null model, **Column name of observation/id variable** - name of a column with observations or ids, and **Column name of outcome variable**.\n- When **Kinship file** and **Phenotype file** are provided on the input, **STAAR Null Model** processes data in following way: the tool will load both files, find the intersection between them based on sample ids and create a null model for intersected samples using defined input parameters for outcome variable, covariates and outcome type.\n\n### Performance Benchmarking\n\nFor performance benchmark **Phenotype files** with 10000, 25000 and 50000 samples were used. The tasks were run with and without the kinship file. Standard set up was run with 15 and simple set up with 6 covariates. Default set up was run with only one covariate. \n      \n\n| Experiment type  | Duration | Cost (Instances + Attached Disks) | Instance (AWS and Google on-demand)|\n|---------------------------|------------------------|-----------------------|--------------------------------|\n| 10k default with kinship | 2 min | $0.02 ($0.01 + $0.01) | c5.2xlarge - 1024GB EBS |\n| 10k default without kinship | 2 min | $0.02 ($0.01 + $0.01) | c5.2xlarge - 1024GB EBS |\n| 25k standard with kinship | 20 min | $0.14 ($0.13 + $0.01) |  n1-standard-8 - 160GB disk |\n| 25k standard without kinship | 3 min | $0.02 | n1-standard-8 - 160GB disk  |\n| 25k simple with kinship | 5 min | $0.04 | n1-standard-8 - 160GB disk  |\n| 25k simple without kinship | 3 min | $0.02 | n1-standard-8 - 160GB disk  |\n| 50k default with kinship | 20 min | $0.16 ($0.12 + $0.04) | c5.2xlarge - 1024GB EBS |\n| 50k default without kinship | 2 min | $0.02 ($0.01 + $0.01) | c5.2xlarge - 1024GB EBS |\n\n\n*Cost can be significantly reduced by spot instance usage. Visit the [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1] [STAAR GitHub](https://github.com/xihaoli/STAAR)\\\n[2] [STAAR_workflow GitHub](https://github.com/sheilagaynor/STAAR_workflow)\\\n[3] [STAAR manual](https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf)",
        "label": "STAAR Null Model",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 100,
            "valueFrom": "2> staar_errors.log"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 6,
            "valueFrom": "${\n    if (!inputs.in_kinship_file){\n        return 'NA'\n    }\n    else{\n        return ''\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_per_job) {\n        return inputs.mem_per_job\n    }\n    else {\n        return 16000\n    }\n}",
            "coresMin": "${\n    if (inputs.cpu_per_job) {\n        return inputs.cpu_per_job\n    }\n    else {\n        return 8\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/lea_lenhardt_ackovic/staar-0-9-5:1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entryname": "STAAR_null_model.R",
                "entry": "# Description: Generate a null model using the STAAR package which provides a wrapper of the GMMAT package.\n# Inputs:\n# pheno_file : file containing the outcome, covariates for the null model (.csv)\n# null_file_name : string containing prefix for .Rds output from null model fitting via STAAR (string)\n# sample_name : column name in pheno_file for observation IDs (string)\n# outcome_name : column name in pheno_file for outcome (string)\n# outcome_type : type of variable of outcome, outcome_name in pheno_file, 'continuous' or 'dichotomous' (string)\n# covariate_names : column names in pheno_file of covariate variables, as comma-separated string, to be treated as covariates (string)\n# kinship_file : file containing the kinship matrix for null model with relatedness, row names are sample_names (.Rds, .Rdata, .csv)\n# het_var_name : column name in pheno_file of variable for grouping heteroscedastic errors (string)\n# Outputs:\n# Null model STAAR object saved as null_file.Rds\n\n## Parse arguments\nargs <- commandArgs(T)\n\n## Required arguments\n# File inputs\npheno_file <- args[1]\nnull_file <- args[2]\n# Analysis inputs\nsample_name <- args[3]\noutcome_name <- args[4]\noutcome_type <- args[5] #default 'continuous'\n# Optional inputs\ncovariate_names <- args[6] #default 'NA'\nkinship_file <- args[7] #default 'NA'\nhet_var_name <- args[8] #default 'NA'\n\n#####################\n# Functions for input processing\n# Adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_nullmodel.R \nget_family <- function(outcome_type_from_args) {\n  if (toupper(outcome_type_from_args) == \"CONTINUOUS\"){\n    family <- gaussian()\n  } else if (toupper(outcome_type_from_args) == \"DICHOTOMOUS\"){\n    family <- binomial()\n  } else {\n    stop_msg <- paste(\"Invalid outcome type provided:\", outcome_type_from_args)\n    stop(stop_msg)  }\n  return(family)\n}\nget_kinship <- function(kinship_file_from_args){\n  cat('Loading Kinship Matrix: ',kinship_file_from_args,'\\n')\n  if (grepl('Rda$',kinship_file_from_args,ignore.case=TRUE) | grepl('Rdata$',kinship_file_from_args,ignore.case=TRUE)){\n    kins <- get(load(kinship_file_from_args))\n  } else if (grepl('Rds$',kinship_file_from_args,ignore.case=TRUE)){\n    kins <- readRDS(kinship_file_from_args) \n  } else if (grepl('csv$',kinship_file_from_args,ignore.case=TRUE)){\n    kins <- as.matrix(read.csv(kinship_file_from_args,as.is=T,check.names=F,row.names=1))\n  } else {\n    stop_msg <- paste(\"Invalid kinship file provided:\", kinship_file_from_args)\n    stop(stop_msg)  }\n  cat('Loaded Kinship: no. rows:',nrow(kins),' no. cols:',ncol(kins),'\\n')\n  kins\n}\n\n# Load packages\nsuppressMessages(library(STAAR))\n\n#####################\n# Read phenotypes, covariates\nif (!grepl('csv$',pheno_file,ignore.case=TRUE)){\n  stop_msg <- paste(\"Invalid phenotype file provided:\", pheno_file)\n  stop(stop_msg)\n}\npheno <- read.csv(pheno_file, header=TRUE, as.is=TRUE)\npheno <- pheno[order(pheno[,sample_name]),]\ncat('Loaded Phenotypes: no. rows:',nrow(pheno),' no. cols:',ncol(pheno),'\\n')\n# Subset to complete cases for phenotype file, create null model formula\nif ( covariate_names=='NA' & het_var_name=='NA' ){\n  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name)]), c(sample_name,outcome_name)]\n  null_model_char <- paste0(outcome_name, \"~\", 1)\n} else if ( covariate_names=='NA' & het_var_name!='NA' ){\n  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,het_var_name)]), c(sample_name,outcome_name,het_var_name)]\n  null_model_char <- paste0(outcome_name, \"~\", 1)\n} else if ( covariate_names!='NA' & het_var_name=='NA' ){\n  covar_split <- unlist(strsplit(covariate_names, split=\",\"))\n  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,covar_split)]), c(sample_name,outcome_name,covar_split)]\n  null_model_char <- paste0(outcome_name, \"~\", paste(covar_split, collapse=\"+\"))\n} else {\n  covar_split <- unlist(strsplit(covariate_names, split=\",\"))\n  phenotype_input <- pheno[complete.cases(pheno[,c(sample_name,outcome_name,het_var_name,covar_split)]), \n                           c(sample_name,outcome_name,het_var_name,covar_split)]\n  null_model_char <- paste0(outcome_name, \"~\", paste(covar_split, collapse=\"+\"))\n}\ncat('Complete Phenotypes: no. rows:',nrow(phenotype_input),' no. cols:',ncol(phenotype_input),'\\n')\n\n#####################\n# Get kinship\nif ( kinship_file!='NA' ){\n  kinship_input <- get_kinship(kinship_file)\n  shared_obs <- intersect(row.names(kinship_input), phenotype_input[,sample_name])\n  if (length(shared_obs)==0){\n    stop('No shared observations between phenotypes and kinship')\n  } else {\n    shared_samples <- intersect(rownames(kinship_input), phenotype_input[,sample_name])\n    kinship_analysis <- kinship_input[match(shared_samples,rownames(kinship_input)),match(shared_samples,rownames(kinship_input))]\n    phenotype_analysis <- phenotype_input[match(shared_samples,phenotype_input[,sample_name]),] \n    rm(kinship_input); rm(phenotype_input)\n    cat('Matched Phenotypes with Kinship: no. rows:',nrow(phenotype_analysis),' no. cols:',ncol(phenotype_analysis),'\\n')\n  }\n} else {\n  phenotype_analysis <- phenotype_input\n  rm(phenotype_input)\n}\n\n#####################\n# Fit, save null model\nif ( kinship_file=='NA' & het_var_name=='NA' ){\n  cat('Fitting null model for unrelated samples, homogeneous variance')\n  null_model <- STAAR::fit_null_glm(as.formula(null_model_char), \n                data = phenotype_analysis, family = get_family(outcome_type))\n} else if ( kinship_file=='NA' & het_var_name!='NA' ){\n  cat('Fitting null model for unrelated samples, heterogeneous variance')\n  null_model <- STAAR::fit_null_glm(as.formula(null_model_char), \n                groups = het_var_name, data = phenotype_analysis, family = get_family(outcome_type))\n} else if ( kinship_file!='NA' & het_var_name=='NA' ){\n  cat('Fitting null model for related samples, homogeneous variance')\n  null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_name, \n                data = phenotype_analysis, kins = kinship_analysis, family = get_family(outcome_type))\n} else {\n  cat('Fitting null model for related samples, heterogeneous variance')\n  null_model <- STAAR::fit_null_glmmkin(as.formula(null_model_char), id = sample_name, \n                groups = het_var_name, data = phenotype_analysis, kins = kinship_analysis, family = get_family(outcome_type))\n}\nsaveRDS(null_model, file=paste0(null_file,\".Rds\"))\n",
                "writable": false
              }
            ]
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:projectName": "STAAR 0.9.5 - Dev",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601359715,
            "sbg:revisionNotes": null
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601359766,
            "sbg:revisionNotes": "First revision"
          },
          {
            "sbg:revision": 2,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601374770,
            "sbg:revisionNotes": "Changed default value for mem_per_job"
          },
          {
            "sbg:revision": 3,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601452426,
            "sbg:revisionNotes": "Description in progress"
          },
          {
            "sbg:revision": 4,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601458202,
            "sbg:revisionNotes": "Description in progress"
          },
          {
            "sbg:revision": 5,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601459278,
            "sbg:revisionNotes": "Licence added"
          },
          {
            "sbg:revision": 6,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601551985,
            "sbg:revisionNotes": "Inherit metadata on output from pheno_file"
          },
          {
            "sbg:revision": 7,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603715602,
            "sbg:revisionNotes": "Wrapper license updated"
          },
          {
            "sbg:revision": 8,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603793175,
            "sbg:revisionNotes": "Added one more format for kinship  file"
          },
          {
            "sbg:revision": 9,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603803079,
            "sbg:revisionNotes": "Edited input format"
          },
          {
            "sbg:revision": 10,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1612960719,
            "sbg:revisionNotes": "Updated to the latest version from GitHub."
          },
          {
            "sbg:revision": 11,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1613650620,
            "sbg:revisionNotes": "Updated R script"
          },
          {
            "sbg:revision": 12,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1614072017,
            "sbg:revisionNotes": "Updated to the latest version from GutHub."
          },
          {
            "sbg:revision": 13,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1614076591,
            "sbg:revisionNotes": "Docker image updated."
          },
          {
            "sbg:revision": 14,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630305818,
            "sbg:revisionNotes": "Script updated to the latest version."
          },
          {
            "sbg:revision": 15,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630306150,
            "sbg:revisionNotes": "Author name corrected."
          },
          {
            "sbg:revision": 16,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630394472,
            "sbg:revisionNotes": "cpu and mem requirements changed."
          },
          {
            "sbg:revision": 17,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630479079,
            "sbg:revisionNotes": "Edited description."
          },
          {
            "sbg:revision": 18,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630479338,
            "sbg:revisionNotes": "Description edited."
          },
          {
            "sbg:revision": 19,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630495926,
            "sbg:revisionNotes": "Benchmark results  added."
          },
          {
            "sbg:revision": 20,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630577634,
            "sbg:revisionNotes": "Outcome type changed to not required."
          },
          {
            "sbg:revision": 21,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630674675,
            "sbg:revisionNotes": "Edited for cwltool."
          },
          {
            "sbg:revision": 22,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630675041,
            "sbg:revisionNotes": "Removed  \"null\": null due to cwltool error."
          }
        ],
        "sbg:image_url": null,
        "sbg:toolkit": "STAAR",
        "sbg:toolkitVersion": "0.9.5",
        "sbg:license": "GPLv3, MIT License",
        "sbg:links": [
          {
            "id": "https://github.com/xihaoli/STAAR",
            "label": "STAAR Homepage"
          },
          {
            "id": "https://github.com/sheilagaynor/STAAR_workflow",
            "label": "STAAR workflow Homepage"
          },
          {
            "id": "https://github.com/xihaoli/STAAR/releases/tag/v0.9.5",
            "label": "Source Code"
          },
          {
            "id": "https://github.com/xihaoli/STAAR/archive/refs/tags/v0.9.5.tar.gz",
            "label": "Download"
          },
          {
            "id": "https://www.nature.com/articles/s41588-020-0676-4",
            "label": "Publication"
          },
          {
            "id": "https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf",
            "label": "Documentation"
          },
          {
            "id": "https://github.com/sheilagaynor/STAAR_workflow/blob/main/STAAR_null_model.R",
            "label": "STAAR null model Code"
          }
        ],
        "sbg:wrapperLicense": "MIT License",
        "sbg:toolAuthor": "Xihao Li, Zilin Li, Han Chen",
        "sbg:wrapperAuthor": "Sheila Gaynor",
        "sbg:categories": [
          "GWAS",
          "CWL1.2"
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-fit-null-model-0-9-5/22",
        "sbg:revision": 22,
        "sbg:revisionNotes": "Removed  \"null\": null due to cwltool error.",
        "sbg:modifiedOn": 1630675041,
        "sbg:modifiedBy": "lea_lenhardt_ackovic",
        "sbg:createdOn": 1601359715,
        "sbg:createdBy": "lea_lenhardt_ackovic",
        "sbg:project": "lea_lenhardt_ackovic/staar-0-9-5-dev",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "lea_lenhardt_ackovic"
        ],
        "sbg:latestRevision": 22,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a3d1a695efb3e472146967297d8590db9b660b9d171933f061bbea93f9b972f58"
      },
      "label": "STAAR Null Model",
      "sbg:x": -300,
      "sbg:y": -120
    },
    {
      "id": "staar_association_test_0_9_5",
      "in": [
        {
          "id": "in_geno_file",
          "source": "in_geno_file"
        },
        {
          "id": "in_null_model",
          "source": "staar_fit_null_model_0_9_5/null_model"
        },
        {
          "id": "in_annot_file",
          "source": "in_annot_file"
        },
        {
          "id": "in_agg_file",
          "source": "in_agg_file"
        },
        {
          "id": "results_filename",
          "source": "results_filename"
        },
        {
          "id": "agds_filename",
          "source": "agds_filename"
        },
        {
          "id": "agds_annot_channels",
          "source": "agds_annot_channels"
        },
        {
          "id": "maf_thres",
          "source": "maf_thres"
        },
        {
          "id": "mac_thres",
          "source": "mac_thres"
        },
        {
          "id": "window_length",
          "source": "window_length"
        },
        {
          "id": "step_length",
          "source": "step_length"
        },
        {
          "id": "num_iterations",
          "source": "num_iterations"
        },
        {
          "id": "num_cores",
          "source": "num_cores"
        },
        {
          "id": "mem_per_job",
          "source": "mem_per_job"
        },
        {
          "id": "candidate_file",
          "source": "candidate_file"
        },
        {
          "id": "in_cond_file",
          "source": "in_cond_file"
        },
        {
          "id": "in_cond_geno_files",
          "source": [
            "in_cond_geno_files"
          ]
        }
      ],
      "out": [
        {
          "id": "assoc_test_results"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-association-test-0-9-5/34",
        "baseCommand": [],
        "inputs": [
          {
            "id": "in_geno_file",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 1
            },
            "label": "Genotype file",
            "doc": "Genotype file containing genotypes for all individuals from null model, optionally containing the given annotation channels (aGDS).",
            "sbg:fileTypes": "GDS"
          },
          {
            "id": "in_null_model",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 0
            },
            "label": "Null model",
            "doc": "File containing output from null model fitting via STAAR.",
            "sbg:fileTypes": "RDS"
          },
          {
            "id": "in_annot_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 2
            },
            "label": "Annotation file",
            "doc": "File containing annotations as input with columns 'chr', 'pos', 'ref', 'alt'.",
            "sbg:fileTypes": "RDS, RDA, CSV, RDATA"
          },
          {
            "id": "in_agg_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 6
            },
            "label": "Aggregate file",
            "doc": "File containing the aggregation units for set-based analysis with columns 'chr', 'pos', 'ref', 'alt', 'group_id'.",
            "sbg:fileTypes": "RDS, RDA, CSV, RDATA"
          },
          {
            "id": "results_filename",
            "type": "string",
            "inputBinding": {
              "shellQuote": false,
              "position": 3
            },
            "label": "Results filename",
            "doc": "Filename for results file."
          },
          {
            "sbg:toolDefaultValue": "None",
            "id": "agds_filename",
            "type": "string?",
            "inputBinding": {
              "shellQuote": false,
              "position": 4,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 'None'\n    }\n}"
            },
            "label": "Annotated GDS filename",
            "doc": "Filename for annotated GDS file indicating whether input geno is an agds file containing the annotations.",
            "default": "None"
          },
          {
            "id": "agds_annot_channels",
            "type": "string?",
            "inputBinding": {
              "shellQuote": false,
              "position": 5,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 'None'\n    }\n}"
            },
            "label": "Annotation channels",
            "doc": "Comma-separated names of annotation channels in annotation GDS file.",
            "default": "None"
          },
          {
            "sbg:toolDefaultValue": "0.01",
            "id": "maf_thres",
            "type": "float?",
            "inputBinding": {
              "shellQuote": false,
              "position": 10,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 0.01\n    }\n}"
            },
            "label": "Minor allele frequency threshold",
            "doc": "Minor allele frequency threshold. AF threshold below which variants will be considered in rare variant analysis. The default is 0.01.",
            "default": 0.01
          },
          {
            "sbg:toolDefaultValue": "1",
            "id": "mac_thres",
            "type": "float?",
            "inputBinding": {
              "shellQuote": false,
              "position": 11,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 1\n    }\n}"
            },
            "label": "Minor allele count threshold",
            "doc": "Minor allele count threshold. AC threshold above which variants will be considered in rare variant analysis. The default is 1.",
            "default": 1
          },
          {
            "sbg:toolDefaultValue": "2000",
            "id": "window_length",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 12,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 2000\n    }\n}"
            },
            "label": "Window length",
            "doc": "Length of window for region-based analysis. The default value is 2000.",
            "default": 2000
          },
          {
            "sbg:toolDefaultValue": "1000",
            "id": "step_length",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 13,
              "valueFrom": "${\n    if (self){\n        return self\n    }\n    else {\n        return 1000\n    }\n}"
            },
            "label": "Step length",
            "doc": "Length of overlap for region-based analysis. The default value is 1000.",
            "default": 1000
          },
          {
            "sbg:toolDefaultValue": "20",
            "id": "num_iterations",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 15,
              "valueFrom": "${  \n    var cmd = ''\n    if (self){\n        cmd = self\n    }\n    else {\n        cmd = 20\n    }\n    return cmd\n}"
            },
            "label": "Number of iterations",
            "doc": "Number of iterations to run in parallel loop, i.e. how many chunks to split sets into. Default is 20.",
            "default": 20
          },
          {
            "sbg:toolDefaultValue": "8",
            "id": "cpu_per_job",
            "type": "int?",
            "label": "CPUs per job",
            "doc": "Number of CPUs per job. The default is 8."
          },
          {
            "sbg:toolDefaultValue": "8",
            "id": "num_cores",
            "type": "int?",
            "inputBinding": {
              "shellQuote": false,
              "position": 14,
              "valueFrom": "${  \n    var cmd = ''\n    if (self){\n        cmd = self\n    }\n    else {\n        cmd = 8\n    }\n    return cmd\n}"
            },
            "label": "Number of cores",
            "doc": "Total number of cores number of cores to be used in parallelized analysis. The default value is 8.",
            "default": 8
          },
          {
            "sbg:toolDefaultValue": "16000",
            "id": "mem_per_job",
            "type": "int?",
            "label": "Memory per job [MB]",
            "doc": "Memory per job [MB] to be used in task run. The default value is 16000."
          },
          {
            "id": "candidate_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 9
            },
            "label": "Candidate file",
            "doc": "Candidate set file. File containing units/windows for candidate sets of interest with columns 'group_id' or 'chr', 'start', 'end'. This input is required for running the tool in Candidate mode. When provided, it will trigger candidate analysis.",
            "sbg:fileTypes": "RDS, RDA, CSV, RDATA"
          },
          {
            "id": "in_cond_file",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 7
            },
            "label": "Conditional file",
            "doc": "Conditional file. File containing the variants to be conditioned upon with columns 'chr', 'pos', 'ref', 'alt'.",
            "sbg:fileTypes": "RDA, RDS, CSV"
          },
          {
            "id": "in_cond_geno_files",
            "type": "File[]?",
            "inputBinding": {
              "shellQuote": false,
              "position": 8,
              "valueFrom": "${\n    var cmd = ''\n    if (self){\n        for (var i=0; i < self.length; i++){\n            cmd = cmd + self[i].path + ','\n        }\n        return cmd.slice(0, -1)\n    }\n    \n}"
            },
            "label": "Conditional genotype files",
            "doc": "Conditional genotype files. File containing genotypes for all individuals from null model for conditional analysis; often same as Genotype file.",
            "sbg:fileTypes": "GDS"
          }
        ],
        "outputs": [
          {
            "id": "assoc_test_results",
            "doc": "STAAR association test results",
            "label": "STAAR results",
            "type": "File?",
            "outputBinding": {
              "glob": "*.gz",
              "outputEval": "$(inheritMetadata(self, inputs.in_geno_file))"
            },
            "sbg:fileTypes": "GZ"
          }
        ],
        "doc": "**STAAR Association Testing** performs the variant-Set Test for Association using Annotation infoRmation (STAAR)[1]. It uses STAAR procedure as well as SKAT, Burden, and ACAT tests for both continuous and dichotomous traits. The STAAR method incorporates qualitative functional categories and quantitative complementary functional annotations [2]. The tool uses the null model, genotypes, and aggregation units (optional) to run rare variant association analyses.\n\nThe tool offers two analysis modes: genomewide and candidate. In genomewide mode, the association testing is performed on a whole genome while in candidate mode the testing is performed only on the variants subset provided in **Candidate file**, which will be run if **Candidate file** is provided on the input.  When provided with **Aggregate file** the tool will do the aggregate testing on the variants present in the file. The tool also offers the conditional association testing. The conditional analysis is performed when the **Conditional file** is provided on the input. The required inputs to run the tool are: **Null model** and **Genotype file**.\n\nThe tool utilises **STAAR Association Testing** implementation from **STAAR_workflow** [2].\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n### Common use cases:\n\n**STAAR Association Testing** is designed to run association testing using the annotation information and SKAT, Burden, and ACAT tests. The user can run genomewide or candidate association testing, which will be triggered if **Candidate file** is provided on the input. In both modes, it is possible to run conditional association testing, by providing **Conditional file** on the input.\n\n### Common issues and important notes:\n\n- Previously calculated null model must be provided on the **Null model** input in Rds format. **STAAR Null Model** should be used for fitting the null model.\n- **Genotype file** input is required and has to be in gds format. If **Annotation channels** and **Annotated GDS filename** are provided, and the **Annotation file** is not present on the input, the tool will assume that annotations are present in annotated gds **Genotype file**. More information about annotated gds format can be found on the [STAAR workflow GitHub page](https://github.com/sheilagaynor/STAAR_workflow).\n- If **Annotation file** is present on the input, the tool will take those annotations and ignore the ones in the GDS file.\n- When **Annotation file** is provided on the input, the tool will analyse only the variants that are present both in the **Annotation file** and **Genotype file**.\n- If the **Annotated GDS filename** has any other  value than `None`, the tool will assume there are annotations present in the **Genotype file**.\n- **Aggregate file** must be provided in .Rds, .Rdata or .csv format with following columns:\n\n```\n\"chr\",\"pos\",\"ref\",\"alt\",\"group_id\"\n\"1\",100000,\"T\",\"C\",\"EXAMPLE\"\n```\n- When running analysis with an **Aggregate file**, **Number of genome chunks** is the number of units to consider at a time within a parallel loop, so it should be set accordingly.\n- To run conditional analysis, **Conditional file** must be provided in .Rds, .Rdata or .csv format. Example format:\n```\n\"chr\",\"pos\",\"ref\",\"alt\"\n\"1\",100000,\"T\",\"C\"\n```\n- Other example input files can be found in the [test files directory of STAAR workflow GitHub page](https://github.com/sheilagaynor/STAAR_workflow/tree/main/testfiles).\n- **Number of cores** is a total number of cores that will be used for multi-threading. The computing instance will be chosen for the task based on the **Number of cores** value. When **Number of cores** is not provided, the tool will check if **CPUs per job** was defined. In case that both **Number of cores** and **CPUs per job** are not provided by the user, the tool will by default use the instance with 8 cores and set **Number of cores** to 8 to be used for multi-threading.\n- **Number of iterations** input defines number of iterations to run in parallel loop, i.e. how many chunks to split sets into [2]. \n- To run candidate analysis, **Candidate file** must be provided on the input and it is expected to have following three columns: *group_id* or *chr* (chromosome), *start* (start position) and *end* (end position).\n- **STAAR Association Testing** has quite high RAM consumption. RAM consumption is directly related to two input parameters: **Number of cores** and **Number of iterations**. **Number of cores** determines how many parallel processes will be run and **Number of iterations** determines how many chunks will be considered at a time within a parallel loop. So these two inputs should be set up according to the instance the task is run on. For **Number of iterations**, smaller values will require more RAM, but the analysis will finish faster, and higher values will require less RAM, but will prolong the duration. The opposite applies to **Number of cores**, higher number leads to faster execution and higher RAM requirement.\n- When running analysis in *genomewide* mode it is recommended to choose a custom instance in the execution settings tab in accordance with values for **Number of iterations** or **Number of cores**. Based on our tests, an r5.16xlarge instance proved to be most optimal (with **Number of cores** set to 32 and **Number of iterations** set to 20 for genotype file with 1 million variants, with 50000 variants analysed per iteration).\n- The tool can sometimes fail with memory error - despite the fact that instance metrics report more available memory. If the task fails with error \"core dumped\", \"segmentation fault\" or any memory related error, user should decrease **Number of cores** or increase **Number of iterations**. The other solution would be to rerun the task on an instance with more RAM.\n- For bigger files it is recommended to run *genomewide* analysis tasks on on-demand instances (i.e., disable spot instances in \"Settings\" for the project and/or \"App Settings\" for the task). This is because tasks with bigger files are long-running and take place on large/expensive instances, increasing the chance that the instance is pre-empted. At this point, the tool will by default switch to an on-demand instance and will start over the analysis (as the intermediary data is not saved).\n\n\n### Performance Benchmarking\n\nPerformance benchmark was run with **Null model** containing 9316 samples and **Genotype files** for Chr20, Chr10 and Chr1 containing 22 and a half millions, 47 millions and 80 and a half millions variants, respectively. For all the tasks **Number of cores** was set to 64 and **Number of genome chunks** to 100000.\n\n| Experiment type  | Duration | Cost (Instances + Attached Disks) | Instance (GCP on-demand)|\n|---------------------------|------------------------|-----------------------|--------------------------------|\n| Chr20 | 3h 48min | $11.65 ($11.55 + $0.10) | n1-standard-64 - 100GB disk |\n| Chr10 | 12h 18min | $46.89 ($46.58 + $0.31) | n1-highmem-64 - 100GB disk |\n| Chr1 | 1d 3h 40min | $105.54 ($104.84 + $0.70) |  n1-highmem-64 - 100GB disk |\n\n*Cost can be significantly reduced by spot instance usage. Visit the [knowledge center](https://docs.sevenbridges.com/docs/about-spot-instances) for more details.*\n\n### References\n[1] [STAAR GitHub](https://github.com/xihaoli/STAAR)\\\n[2] [STAAR_workflow GitHub](https://github.com/sheilagaynor/STAAR_workflow)\\\n[3] [STAAR manual](https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf)",
        "label": "STAAR Association Testing",
        "arguments": [
          {
            "prefix": "",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "Rscript ./STAAR_analysis.R"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 2,
            "valueFrom": "${\n    if (!inputs.in_annot_file){\n        return 'None'\n    }\n    else{\n        return ''\n    }\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 6,
            "valueFrom": "${\n    if (!inputs.in_agg_file){\n        return 'None'\n    }\n    else{\n        return ''\n    }\n}"
          },
          {
            "prefix": "2>",
            "shellQuote": false,
            "position": 2000,
            "valueFrom": "staar_error.log"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 7,
            "valueFrom": "${\n    if (!inputs.in_cond_file){\n        return 'None'\n    }\n    else{\n        return ''\n    }\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 8,
            "valueFrom": "${\n    if (!inputs.in_cond_geno_files){\n        return 'None'\n    }\n    else{\n        return ''\n    }\n}"
          },
          {
            "prefix": "",
            "shellQuote": false,
            "position": 9,
            "valueFrom": "${\n    if (!inputs.candidate_file){\n        return 'None'\n    }\n    else{\n        return ''\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": "${\n    if (inputs.mem_per_job) {\n        return inputs.mem_per_job\n    }\n    else {\n        return 16000\n    }\n}",
            "coresMin": "${\n    if (inputs.num_cores) {\n        return inputs.num_cores\n    }\n    else if (inputs.cpu_per_job) {\n        return inputs.cpu_per_job\n    }\n    else {\n        return 8\n    }\n}"
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/lea_lenhardt_ackovic/staar-0-9-5:1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entryname": "STAAR_analysis.R",
                "entry": "# Description: Run genomewide analysis using the STAAR package.\n# Inputs:\n# null_file : file containing output from null model fitting via STAAR (.Rds)\n# geno_file : file containing genotypes for at least all individuals from null model, optionally containing the given annotation channels (.gds)\n# annot_file : file containing annotations as input with columns 'chr', 'pos', 'ref', 'alt' (.Rds, .Rdata, .csv)\n# results_file : string of name of results file output (string)\n# agds_file : string indicating whether input geno is an agds file containing the annotations, 'None' if not provided (string)\n# agds_annot_channels : comma-separated names of channels in agds to be treated as annotations (string)\n# agg_file : file containing the aggregation units for set-based analysis with columns 'chr', 'pos', 'ref', 'alt', 'group_id' (.Rds, .Rdata, .csv)\n# cond_file : file containing the variants to be conditioned upon with columns 'chr', 'pos', 'ref', 'alt' (.Rds, .Rdata, .csv)\n# cond_geno_files : files containing genotypes for at least all individuals from null model for conditional analysis; often same as geno_file (.gds)\n# cand_file : file containing units/windows for candidate sets of interest with columns 'group_id' or 'chr', 'start', 'end' (.Rds, .Rdata, .csv)\n# maf_thres : AF threshold below which variants will be considered in rare variant analysis, 0.01 default (numeric)\n# mac_thres : AC threshold above which variants will be considered in rare variant analysis, 1 default (numeric)\n# window_length : length of window for region-based analysis, 2000 default (numeric)\n# step_length : length of overlap for region-based analysis, 1000 default (numeric)\n# num_cores : number of cores to be used in parallelized analysis (numeric)\n# num_iterations : number of iterations to run in parallel loop, i.e. how many chunks to split sets into, 20 [default] (numeric)\n\nSys.setenv(MKL_NUM_THREADS = 1)\n\n## Parse arguments\nargs <- commandArgs(T)\n\n## Required arguments\n# File inputs\nnull_file <- args[1]\ngeno_file <- args[2]\nannot_file <- args[3]\nresults_file <- args[4]\nagds_file <- args[5]\nagds_annot_channels <- args[6]\nagg_file <- args[7]\ncond_file <- args[8]\ncond_geno_files <- args[9]\ncand_file <- args[10]\n# Analysis inputs\nmaf_thres <- as.numeric(args[11])\nmac_thres <- as.numeric(args[12])\nwindow_length <- as.numeric(args[13])\nstep_length <- as.numeric(args[14])\n# Compute inputs\nnum_cores <- as.numeric(args[15])\nnum_iterations <- as.numeric(args[16])\n\n#####################\n# Function for input file processing\nget_file <- function(file_from_args,file_type){\n  # Read in file\n  cat('Loading File: ',file_from_args,'\\n')\n  if (grepl('Rda$',file_from_args,ignore.case=TRUE) | grepl('Rdata$',file_from_args,ignore.case=TRUE)){\n    file_in <- get(load(file_from_args))\n  } else if (grepl('Rds$',file_from_args,ignore.case=TRUE)){\n    file_in <- readRDS(file_from_args) \n  } else {\n    file_in <- fread(file_from_args,stringsAsFactors=F,sep=',',header=T,data.table=F)\n  }\n  cat('Loaded ', file_type,' file: no. rows:',nrow(file_in),' no. cols:',ncol(file_in),'\\n')\n  if (sum(names(file_in) %in% c('CHR','CHROM','Chr','chrom'))>0){\n    names(file_in)[names(file_in) %in%  c('CHR','CHROM','Chr','chrom')] = 'chr'\t\n    file_in$chr = sub('chr','',file_in$chr)\n  }\n  names(file_in) = tolower(names(file_in))\n  # Check that required columns are available\n  if (file_type=='Annotation' | file_type=='Conditional'){\n    if ( !(sum(names(file_in) %in% c('chr','pos','ref','alt'))==4) ){\n      stop(paste0(file_type, \" file does not provide necessary input \\n\")) }\n  }\n  if (file_type=='Aggregation'){\n    if ( !(sum(names(file_in) %in% c('chr','pos','ref','alt','group_id'))==5) ){\n      stop(paste0(file_type, \" file does not provide necessary input \\n\")) }\n  }\n  if (file_type=='Candidate'){\n    if ( !(sum(names(file_in) %in% c('chr','start','end'))==3 | any(names(file_in)=='group_id')) ){\n      stop(paste0(file_type, \" file does not provide necessary input \\n\")) }\n  }\n  file_in\n}\n#https://github.com/UW-GAC/analysis_pipeline/blob/master/TopmedPipeline/R/aggregateList.R\n#aggregateGRangesList function Credit to: smgogarten \naggregateGRangesList <- function(variants) {\n  stopifnot(all(c(\"group_id\", \"chr\", \"pos\") %in% names(variants)))\n  groups <- unique(variants$group_id)\n  cols <- setdiff(names(variants), c(\"group_id\", \"chr\", \"pos\"))\n  GRangesList(lapply(setNames(groups, groups), function(g) {\n    x <- variants[variants$group_id == g,]\n    gr <- GRanges(seqnames=x$chr, ranges=IRanges(start=x$pos, width=1))\n    mcols(gr) <- x[,cols]\n    gr\n  }))\n}\n\n# Load packages\nsuppressMessages(library(gdsfmt))\nsuppressMessages(library(SeqArray))\nsuppressMessages(library(STAAR))\nsuppressMessages(library(SeqVarTools))\nsuppressMessages(library(dplyr))\nsuppressMessages(library(doMC))\nsuppressMessages(library(GenomicRanges))\nsuppressMessages(library(R.utils))\nsuppressMessages(library(data.table))\n\n# Read in files: null model, genotypes\nnull_model <- readRDS(null_file)\ngeno_all <- seqOpen(geno_file)\ncat('Read in provided null model and genotype files \\n')\n# Read in files: optional files as provided\nif (annot_file=='None' & agds_file=='None'){\n  cat('No annotation file or aGDS provided, proceeding without annotations \\n')\n} else if (annot_file!='None' & agds_file=='None') {\n  annot_table <- get_file(annot_file,'Annotation')\n}\nif (agg_file!='None'){\n  cat('Aggregate testing based on grouping file will be used \\n')\n  agg_units <- get_file(agg_file,'Aggregation')\n} else {\n  window_length <- ifelse(window_length>0, window_length, 2000)\n  step_length <- ifelse(step_length>0, step_length, 1000)\n  cat(paste0('Sliding window testing will be done with window length: ',window_length,', step length: ',step_length, ' \\n'))\n}\nif (cand_file!='None'){\n  cand_in <- get_file(cand_file,'Candidate')\n  if (agg_file!='None' & any(names(cand_in)=='group_id')){\n    agg_units <- agg_units[agg_units$group_id %in% cand_in$group_id,]\n    cat(\"Candidate aggregate testing based on groups to be run \\n\")\n  } else {\n    cat(\"Candidate sliding window testing to be run \\n\")\n  }\n}\n\n# Set up potential multi-core analysis\nn_cores <- min(c(num_cores, parallel::detectCores(logical = TRUE)))\n\n#####################\n# Define sample and variants of interest\nif(any(class(null_model)=='glmmkin')){\n  pheno_id <- as.character(null_model$id_include)\n} else {\n  pheno_id <- as.character(null_model$data$ID)\n}\nvariant_id <- seqGetData(geno_all, \"variant.id\")\nif (agds_file!='None'){\n  filter <- seqGetData(geno_all, \"annotation/filter\")\n  AVGDP <- seqGetData(geno_all, \"annotation/info/AVGDP\")\n  SNVlist <- filter == \"PASS\" & AVGDP > 10 & isSNV(geno_all)\n  rm(filter); rm(AVGDP)\n  seqSetFilter(geno_all,sample.id=pheno_id,variant.id=variant_id[SNVlist])\n} else {\n  seqSetFilter(geno_all,sample.id=pheno_id,variant.id=variant_id)\n}\n\n#Below adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R\nvariant_info <- variantInfo(geno_all, alleles = FALSE, expanded=FALSE)\nchr <- variant_info$chr[1]\nif(length(unique(chr)) > 1) stop(\"Multiple chromosomes detected; terminating \\n\")\nchr <- chr[1] \n#Get the aggregation units\nchunk <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))\nif(agg_file!='None'){\n  agg_units <- agg_units[agg_units$chr == chr,]\n  agg_names <- unique(agg_units$group_id)\n  agg_chunks <- chunk(agg_names,min(num_iterations,length(agg_names)))\n  n_chunk <- length(agg_chunks)\n} else {\n  if (cand_file=='None'){\n    window_start <- seq(min(variant_info$pos), max(variant_info$pos), step_length)\n    windows <- GRanges(seqnames=chr, ranges=IRanges(start=window_start, width=window_length))\n    window_chunks <- chunk(1:length(windows),min(num_iterations,length(windows)))\n    n_chunk <- length(window_chunks)\n  } else {\n    #Get the genome chunks for windows from candidate input\n    grange_df <- data.frame(chr=cand_in$chr, start=cand_in$start, end=cand_in$end)\n    range_data <- makeGRangesFromDataFrame(grange_df)\n    n_chunk <- length(range_data) ; rm(grange_df)\n  }\n}\nseqClose(geno_all)\ngc()\n\n#####################\n# Define, prepare conditional set\nif (cond_file != 'None'){\n  cond_in <- get_file(cond_file,'Conditional')\n  cond_geno_vec <- unlist(strsplit(cond_geno_files,','))\n  cond_matrix <- c()\n  for (cond_geno_file in cond_geno_vec){\n    cond_geno <- seqOpen(cond_geno_file)\n    cond_geno_variant_info_all <- variantInfo(cond_geno, alleles = TRUE, expanded=FALSE)\n    cond_geno_variant_info <- cond_geno_variant_info_all[cond_geno_variant_info_all$pos %in% cond_in$pos,]\n    avail_cond_geno <- merge(cond_geno_variant_info, cond_in, by=c('chr','pos','ref','alt'))\n    cond_var_list <- cond_geno_variant_info_all$variant.id[cond_geno_variant_info_all$variant.id %in% avail_cond_geno$variant.id]\n    rm(cond_geno_variant_info_all); rm(cond_geno_variant_info); rm(avail_cond_geno)\n    seqSetFilter(cond_geno,sample.id=pheno_id,variant.id=cond_var_list)\n    cond_id.genotype.match <- match(pheno_id, seqGetData(cond_geno,\"sample.id\"))\n    cond_genotypes <- seqGetData(cond_geno, \"$dosage\")\n    cond_matrix <- cbind(cond_matrix, cond_genotypes[cond_id.genotype.match,])\n    rm(cond_var_list); rm(cond_id.genotype.match); rm(cond_genotypes)\n    seqClose(cond_geno)\n  }\n  cat('Prepared conditional variant lists: no. conditioning variants:',ncol(cond_matrix),'\\n')\n}\ngc()\n\n\n#####################\n# Function for running tests on the large chunk\ntest_chunk <- function( indx ){\n  print(paste0('Iteration: ', indx))\n  #Read in genotype data\n  geno <- seqOpen(geno_file)\n  seqSetFilter(geno,sample.id=pheno_id,verbose=F)\t\n  \n  #First for gene based/agg unit based (candidate or full chromosome)\n  #Agg unit option adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R\n  if(agg_file!='None'){\n    if (agds_file!='None'){\n      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_id[SNVlist],verbose=F)\t\n    }\n    agg_var <- aggregateGRangesList(agg_units[agg_units$group_id %in% agg_chunks[[indx]],])\n    if( length(names(agg_var)) < 1) return(data.frame())\n    seqSetFilter(geno, variant.sel = agg_var, verbose = TRUE)\n    #Subset to rare variants for efficiency\n    chunk_variant_id <- seqGetData(geno,'variant.id')\n    if (length(chunk_variant_id)>1){\n      n_samp <- length(seqGetData(geno,\"sample.id\"))\n      ct_vec <- seqAlleleCount(geno)\n      min_ct_vec <- pmin(ct_vec, (2*n_samp) - ct_vec)\n      min_freq_vec <- min_ct_vec / (2*n_samp)\n      rare_inc <- (min_freq_vec < maf_thres) & (min_ct_vec != 0) & (min_ct_vec > mac_thres)\n      seqSetFilter(geno, variant.id=chunk_variant_id[rare_inc], verbose = TRUE)\n      rm(ct_vec); rm(min_ct_vec); rm(min_freq_vec); rm(n_samp); rm(rare_inc)\n      #Subset annotations for efficiency [retains position; intersection not required here]\n      if(annot_file!='None'){\n        chunk_variant_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)\n        geno_matching <- paste(chunk_variant_id$chr, chunk_variant_id$pos, chunk_variant_id$ref, chunk_variant_id$alt, sep='_')\n        annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')\n        annot_chunk <- annot_table[annot_matching %in% geno_matching,]\n      }\n      #Get iterator for looping tests\n      if(length(seqGetData(geno, \"variant.id\"))>1){\n        iterator <- SeqVarListIterator(geno, agg_var, verbose = T)\n        var_info_iter <- list(variantInfo(iterator))\n        iter <- 2\n        while(iterateFilter(iterator)) {\n          var_info_iter[[iter]] <- variantInfo(iterator) ; iter <- iter + 1\n        }\n        agg_var <- agg_var[sapply(var_info_iter, function(x) nrow(x)>0)]\n        var_info_iter <- Filter(function(x) nrow(x) > 0, var_info_iter)\n        results <- c()\n        for ( var_set in 1:length(var_info_iter)) {\n          seqSetFilter(geno, sample.id=pheno_id, variant.id=var_info_iter[[var_set]]$variant.id, verbose = TRUE)\n          if (annot_file=='None' & agds_annot_channels=='None'){\n            #Proceed without annotations\n            #Subset to the genotypes of interest, match to phenotypes\n            sample.id.match <- match(pheno_id, seqGetData(geno,\"sample.id\"))\n            genotypes <- seqGetData(geno, \"$dosage\")\n            genotypes <- genotypes[sample.id.match,]\n            pvalues <- 0\n            if(cond_file=='None'){\n              try(pvalues <- STAAR(genotypes,null_model))\n            } else {\n              try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model))\n            }\n          } else if(annot_file=='None' & agds_annot_channels!='None') {\n            #Proceed with aGDS annotations\n            #Subset to the genotypes of interest, match to phenotypes\n            sample.id.match <- match(pheno_id, seqGetData(geno,\"sample.id\"))\n            genotypes <- seqGetData(geno, \"$dosage\")\n            genotypes <- genotypes[sample.id.match,]\n            ## Get annotations\n            annot_str_spl <- unlist(strsplit(agds_annot_channels, split=\",\"))\n            annot_tab <- c()\n            for (kk in 1:length(annot_str_spl)){\n              annot_tab <- cbind(annot_tab, seqGetData(geno, annot_str_spl[kk]))\n            }\n            annot_chunk <- data.frame(annot_tab) ; rm(annot_tab)\n            names(annot_chunk) <- annot_str_spl\n            ## Remove missing values\n            if (is.null(dim(annot_chunk))){\n              non_miss <- !is.na(annot_chunk)\n              annot_chunk <- annot_chunk[non_miss]\n            } else {\n              non_miss <- rowSums(is.na(annot_chunk)) == 0\n              annot_chunk <- annot_chunk[non_miss,]\n            }\n            pvalues <- 0\n            if (!is.null(dim(genotypes))){\n              genotypes <- genotypes[,non_miss]\n              if(cond_file=='None'){\n                try(pvalues <- STAAR(genotypes,null_model,annot_chunk))\n              } else {\n                try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model,annot_chunk))\n              }\n            }\n          } else {\n            #Proceed with annotations\n            #Subset to the genotypes, annots of interest, match to phenotypes\n            variant_info_var_set <- variantInfo(geno, alleles = TRUE, expanded=FALSE)\n            sample.id.match <- match(pheno_id, seqGetData(geno,\"sample.id\"))\n            genotypes <- seqGetData(geno, \"$dosage\")\n            genotypes <- genotypes[sample.id.match,]\n            geno_matching <- paste(variant_info_var_set$chr, variant_info_var_set$pos, variant_info_var_set$ref, variant_info_var_set$alt, sep='_')\n            annot_matching <- paste(annot_chunk$chr, annot_chunk$pos, annot_chunk$ref, annot_chunk$alt, sep='_')\n            geno_annot_var <- intersect(geno_matching, annot_matching)\n            annot_var_set <- annot_chunk[match(geno_annot_var,annot_matching),]\n            genotypes_var_set <- genotypes[,match(geno_annot_var,geno_matching)]\n            if (!is.null(dim(annot_var_set))){\n              annot_var_set <- annot_var_set[,!(names(annot_var_set) %in% c(\"variant.id\",\"chr\",\"pos\",\"ref\",\"alt\"))]\n              pvalues <- 0\n              if(cond_file=='None'){\n                try(pvalues <- STAAR(genotypes_var_set,null_model,annot_var_set))\n              } else {\n                try(pvalues <- STAAR_cond(genotypes_var_set,cond_matrix,null_model,annot_var_set))\n              }\n            }\n          }\n          if(class(pvalues)==\"list\"){\n            results_temp <- c(chr, names(agg_var[var_set]), unlist(pvalues[-2]))\n            results <- rbind(results,results_temp)\n          } \n          seqResetFilter(geno)\n        }\n      }\n    }\n  }\n\n  #Next for window based\n  if(agg_file=='None' & cand_file=='None'){\n    #Extract variants in region\n    variant_info <- variantInfo(geno, alleles = FALSE, expanded=FALSE)\n    current_windows <- windows[window_chunks[[indx]]]\n    indx_vars <- (variant_info$pos>=start(current_windows)[1]) & (variant_info$pos<=tail(end(current_windows),1))\n    if (agds_file!='None'){\n      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info$variant.id[SNVlist & indx_vars])\n    } else {\n      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info$variant.id[indx_vars])\n    }\n    #Subset to rare variants for efficiency or break out\n    chunk_variant_id <- seqGetData(geno,'variant.id')\n    if (length(chunk_variant_id)>1) {\n      n_samp <- length(seqGetData(geno,\"sample.id\"))\n      ct_vec <- seqAlleleCount(geno)\n      min_ct_vec <- pmin(ct_vec, (2*n_samp) - ct_vec)\n      min_freq_vec <- min_ct_vec / (2*n_samp)\n      rare_inc <- (min_freq_vec < maf_thres) & (min_ct_vec != 0) & (min_ct_vec > mac_thres)\n      seqSetFilter(geno, variant.id=chunk_variant_id[rare_inc], verbose = TRUE)\n      rm(ct_vec); rm(min_ct_vec); rm(min_freq_vec); rm(n_samp); rm(rare_inc)\n      # Match the genotype and phenotype ids\n      sample.id.match <- match(pheno_id, seqGetData(geno,\"sample.id\"))\n      genotypes <- seqGetData(geno, \"$dosage\")\n      genotypes <- genotypes[sample.id.match,]\n      geno_variant_rare_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)\n      # Get annotations if provided\n      if(annot_file=='None' & agds_annot_channels!='None') {\n        # Get annotations from agds\n        annot_str_spl <- unlist(strsplit(agds_annot_channels, split=\",\"))\n        annot_tab <- c()\n        for (kk in 1:length(annot_str_spl)){\n          annot_tab <- cbind(annot_tab, seqGetData(geno, annot_str_spl[kk]))\n        }\n        annot_chunk <- data.frame(annot_tab)\n        names(annot_chunk) <- annot_str_spl\n      }  else if(annot_file!='None'){ \n        geno_matching <- paste(geno_variant_rare_id$chr, geno_variant_rare_id$pos, geno_variant_rare_id$ref, geno_variant_rare_id$alt, sep='_')\n        annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')\n        geno_annot_var <- intersect(geno_matching, annot_matching)\n        annot_chunk <- annot_table[match(geno_annot_var,annot_matching),]\n        genotypes <- genotypes[,match(geno_annot_var,geno_matching)]\n        geno_variant_rare_id <- geno_variant_rare_id[match(geno_annot_var,geno_matching),]\n      }\n      # Loop through the windows\n      results <- c()\n      for ( window_indx in 1:length(current_windows)) {\n        # Select the region from the geno matrix\n        geno_region <- genotypes[,(geno_variant_rare_id$pos>= start(current_windows[window_indx])) & (geno_variant_rare_id$pos<= end(current_windows[window_indx]))]\n        if (exists(\"annot_chunk\")){\n          # Select annotations from chunk matrix\n          annot_region <- annot_chunk[(geno_variant_rare_id$pos>= start(current_windows[window_indx])) & (geno_variant_rare_id$pos<= end(current_windows[window_indx])),]\n          annot_region <- annot_region[,!(names(annot_region) %in% c(\"variant.id\",\"chr\",\"pos\",\"ref\",\"alt\"))]\n        }\n        pvalues <- 0\n        if(cond_file=='None'){\n          if (annot_file=='None' & agds_annot_channels=='None'){\n            try(pvalues <- STAAR(geno_region,null_model))\n          } else {\n            try(pvalues <- STAAR(geno_region,null_model,annot_region))\n          }\n        } else {\n          if (annot_file=='None' & agds_annot_channels=='None'){\n            try(pvalues <- STAAR_cond(geno_region,cond_matrix,null_model))\n          } else {\n            try(pvalues <- STAAR_cond(geno_region,cond_matrix,null_model,annot_region))\n          }\n        }\n        if(class(pvalues)==\"list\") {\n          results_temp <- c(chr, start(current_windows[window_indx]), end(current_windows[window_indx]), unlist(pvalues[-2]))\n          results <- rbind(results,results_temp)\n        }\n      }\n    }\n  }\n  \n  #Next for candidate window based\n  if(agg_file=='None' & cand_file!='None'){\n    #Extract variants in region\n    variant_info_chunk <- variantInfo(geno, alleles = FALSE, expanded=FALSE)\n    indx_vars <- (variant_info_chunk$pos>=start(range_data[indx]@ranges)) & (variant_info_chunk$pos<=end(range_data[indx]@ranges))\n    if (agds_file!='None'){\n      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info_chunk$variant.id[SNVlist & indx_vars])\n    } else {\n      seqSetFilter(geno,sample.id=pheno_id,variant.id=variant_info_chunk$variant.id[indx_vars])\n    }\n    #Subset to rare variants for efficiency or break out\n    chunk_variant_id <- seqGetData(geno,'variant.id')\n    if (length(chunk_variant_id)>1) {\n      n_samp <- length(seqGetData(geno,\"sample.id\"))\n      ct_vec <- seqAlleleCount(geno)\n      min_ct_vec <- pmin(ct_vec, (2*n_samp) - ct_vec)\n      min_freq_vec <- min_ct_vec / (2*n_samp)\n      rare_inc <- (min_freq_vec < maf_thres) & (min_ct_vec != 0) & (min_ct_vec > mac_thres)\n      seqSetFilter(geno, variant.id=chunk_variant_id[rare_inc], verbose = TRUE)\n      rm(ct_vec); rm(min_ct_vec); rm(min_freq_vec); rm(n_samp); rm(rare_inc)\n      # Match the genotype and phenotype ids\n      sample.id.match <- match(pheno_id, seqGetData(geno,\"sample.id\"))\n      # Get genotype as matrix\n      genotypes <- seqGetData(geno, \"$dosage\")\n      genotypes <- genotypes[sample.id.match,]\n      geno_variant_rare_id <- variantInfo(geno, alleles = TRUE, expanded=FALSE)\n      # Get annotations\n      if(annot_file=='None' & agds_annot_channels!='None') {\n        # Get annotations from agds\n        annot_str_spl <- unlist(strsplit(agds_annot_channels, split=\",\"))\n        annot_tab <- c()\n        for (kk in 1:length(annot_str_spl)){\n          annot_tab <- cbind( annot_tab, seqGetData(geno, annot_str_spl[kk]))\n        }\n        annot_chunk <- data.frame(annot_tab)\n        names(annot_chunk) <- annot_str_spl\n      }  else if(annot_file!='None'){\n        geno_matching <- paste(geno_variant_rare_id$chr, geno_variant_rare_id$pos, geno_variant_rare_id$ref, geno_variant_rare_id$alt, sep='_')\n        annot_matching <- paste(annot_table$chr, annot_table$pos, annot_table$ref, annot_table$alt, sep='_')\n        geno_annot_var <- intersect(geno_matching, annot_matching)\n        annot_chunk <- annot_table[match(geno_annot_var,annot_matching),]\n        genotypes <- genotypes[,match(geno_annot_var,geno_matching)]\n        annot_chunk <- annot_chunk[,!(names(annot_chunk) %in% c(\"variant.id\",\"chr\",\"pos\",\"ref\",\"alt\"))]\n      }\n      #Genome chunks is the window of interest\n      pvalues <- 0\n      if(cond_file=='None'){\n        if (annot_file=='None' & agds_annot_channels=='None'){\n          try(pvalues <- STAAR(genotypes,null_model))\n        } else {\n          try(pvalues <- STAAR(genotypes,null_model,annot_chunk))\n        }\n      } else {\n        if (annot_file=='None' & agds_annot_channels=='None'){\n          try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model))\n        } else {\n          try(pvalues <- STAAR_cond(genotypes,cond_matrix,null_model,annot_chunk))\n        }\n      }\n      if(class(pvalues)==\"list\") {\n        results <- c(chr, start(range_data[indx]@ranges), end(range_data[indx]@ranges), unlist(pvalues[-2]))\n      }\n    }\n  }\n  \n  #Assemble results\n  if(!exists(\"results\")){\n    results <- data.frame()\n  }\n  seqClose(geno)\n  gc()\n  results\n}\n\n#####################\n# Function for running full analysis\n# Below adapted from https://github.com/AnalysisCommons/genesis_wdl/blob/master/genesis_tests.R\nrun_analysis <- function( n_cores ){\n  cat('Running Analysis with ', n_cores,' cores of ', num_cores,' specified cores \\n')\n  print(paste('Running in', n_chunk,' analysis units'))\n  if (n_cores>1){\n    doMC::registerDoMC(cores = n_cores)\n    mc_options <- list(preschedule=FALSE, set.seed=FALSE)\n    out <- foreach(i=1:n_chunk, .combine=rbind, .inorder=FALSE, .options.multicore = mc_options) %dopar% test_chunk(i)\n  } else {\n    out <- c()\n    for (i in 1:n_chunk){\n      out_temp <- test_chunk(i)\n      out <- rbind(out, out_temp)\n    }\n  }\n  cat(\"\\nFinished analysis \\n\")\n  gc()\n  out\n}\n\n\n# Run, clean up results\nresults <- run_analysis( n_cores )\ngc()\nif(!is.null(results)){\n  colnames(results) <- colnames(results, do.NULL = FALSE, prefix = \"col\")\n  if(agg_file!='None'){\n    colnames(results)[1:2] <- c('chr','group_id')\n  } else {\n    colnames(results)[1:3] <- c('chr','start_pos','end_pos')\n  }\n}\n\n\n# Save output\nfwrite(results, file=paste0(results_file,'_chr',chr,\".csv.gz\"))",
                "writable": false
              }
            ]
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:projectName": "STAAR 0.9.5 - Dev",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601359831,
            "sbg:revisionNotes": null
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601359864,
            "sbg:revisionNotes": "First revision"
          },
          {
            "sbg:revision": 2,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601464316,
            "sbg:revisionNotes": "Description in progress"
          },
          {
            "sbg:revision": 3,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601464387,
            "sbg:revisionNotes": "File types added for candidate file"
          },
          {
            "sbg:revision": 4,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601553627,
            "sbg:revisionNotes": "Output inherits metadata from geno_file"
          },
          {
            "sbg:revision": 5,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601890269,
            "sbg:revisionNotes": "Added error log for stdout and error log"
          },
          {
            "sbg:revision": 6,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601890487,
            "sbg:revisionNotes": "Only error log in the end"
          },
          {
            "sbg:revision": 7,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601891323,
            "sbg:revisionNotes": "Edited inputs for chunk_size and num_cores"
          },
          {
            "sbg:revision": 8,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601891580,
            "sbg:revisionNotes": "Just error in error log"
          },
          {
            "sbg:revision": 9,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601891863,
            "sbg:revisionNotes": "Wrong name of R script edited in cmd line"
          },
          {
            "sbg:revision": 10,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1601894502,
            "sbg:revisionNotes": "Updated R scripts, now can run conditional analysis"
          },
          {
            "sbg:revision": 11,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603715572,
            "sbg:revisionNotes": "Wrapper license updated"
          },
          {
            "sbg:revision": 12,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603723248,
            "sbg:revisionNotes": "Candidate script edited row 226 SeqData"
          },
          {
            "sbg:revision": 13,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603793523,
            "sbg:revisionNotes": "Edited file formats for some inputs"
          },
          {
            "sbg:revision": 14,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1603803189,
            "sbg:revisionNotes": "Edited input format"
          },
          {
            "sbg:revision": 15,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1604043786,
            "sbg:revisionNotes": "Edited condition i scripts"
          },
          {
            "sbg:revision": 16,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1604301256,
            "sbg:revisionNotes": "Testing parallel comp"
          },
          {
            "sbg:revision": 17,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1605166649,
            "sbg:revisionNotes": "Corrected genomewide script for parallelization"
          },
          {
            "sbg:revision": 18,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1605177693,
            "sbg:revisionNotes": "Description edited"
          },
          {
            "sbg:revision": 19,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1605177752,
            "sbg:revisionNotes": "Description edited"
          },
          {
            "sbg:revision": 20,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1612960195,
            "sbg:revisionNotes": "Updated to the latest version of scripts on GitHub."
          },
          {
            "sbg:revision": 21,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1613650131,
            "sbg:revisionNotes": "Updated R script"
          },
          {
            "sbg:revision": 22,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1614071870,
            "sbg:revisionNotes": "Updated to the latest version from GitHub."
          },
          {
            "sbg:revision": 23,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1614076610,
            "sbg:revisionNotes": "Docker image updated."
          },
          {
            "sbg:revision": 24,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1619673927,
            "sbg:revisionNotes": "STAAR script updated"
          },
          {
            "sbg:revision": 25,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630306222,
            "sbg:revisionNotes": "R script updated to the latest version."
          },
          {
            "sbg:revision": 26,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630407385,
            "sbg:revisionNotes": "Updated description."
          },
          {
            "sbg:revision": 27,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630478422,
            "sbg:revisionNotes": "Updated description."
          },
          {
            "sbg:revision": 28,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630578337,
            "sbg:revisionNotes": "Description updated."
          },
          {
            "sbg:revision": 29,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630664014,
            "sbg:revisionNotes": "Added benchmark results to description."
          },
          {
            "sbg:revision": 30,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630674944,
            "sbg:revisionNotes": "Edited JS for cwltool."
          },
          {
            "sbg:revision": 31,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1633093769,
            "sbg:revisionNotes": "Updated to script version 1.1."
          },
          {
            "sbg:revision": 32,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1633512081,
            "sbg:revisionNotes": "Updated to the latest STAAR analysis version."
          },
          {
            "sbg:revision": 33,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1634635418,
            "sbg:revisionNotes": "Updated to the latest version of script from Git repo."
          },
          {
            "sbg:revision": 34,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1638257994,
            "sbg:revisionNotes": "Updated description."
          }
        ],
        "sbg:image_url": null,
        "sbg:toolkit": "STAAR",
        "sbg:toolkitVersion": "0.9.5",
        "sbg:license": "GPLv3, MIT License",
        "sbg:links": [
          {
            "id": "https://github.com/xihaoli/STAAR",
            "label": "STAAR Homepage"
          },
          {
            "id": "https://github.com/sheilagaynor/STAAR_workflow",
            "label": "STAAR workflow Homepage"
          },
          {
            "id": "https://github.com/xihaoli/STAAR/releases/tag/v0.9.5",
            "label": "Source Code"
          },
          {
            "id": "https://github.com/xihaoli/STAAR/archive/refs/tags/v0.9.5.tar.gz",
            "label": "Download"
          },
          {
            "id": "https://www.nature.com/articles/s41588-020-0676-4",
            "label": "Publication"
          },
          {
            "id": "https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf",
            "label": "Documentation"
          },
          {
            "id": "https://github.com/sheilagaynor/STAAR_workflow/blob/main/STAAR_analysis.R",
            "label": "STAAR analysis Code"
          }
        ],
        "sbg:wrapperLicense": "MIT License",
        "sbg:wrapperAuthor": "Sheila Gaynor",
        "sbg:toolAuthor": "Xihao Li, Zilin Li, Han Chen",
        "sbg:categories": [
          "GWAS",
          "CWL1.2"
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-association-test-0-9-5/34",
        "sbg:revision": 34,
        "sbg:revisionNotes": "Updated description.",
        "sbg:modifiedOn": 1638257994,
        "sbg:modifiedBy": "lea_lenhardt_ackovic",
        "sbg:createdOn": 1601359831,
        "sbg:createdBy": "lea_lenhardt_ackovic",
        "sbg:project": "lea_lenhardt_ackovic/staar-0-9-5-dev",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "lea_lenhardt_ackovic"
        ],
        "sbg:latestRevision": 34,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a8bf764a364b6d7dde68b1ed066f53f07164febfa2cb1b6a8d4e2457f13522be6"
      },
      "label": "STAAR Association Testing",
      "scatter": [
        "in_geno_file"
      ],
      "sbg:x": -102.57142639160156,
      "sbg:y": 74.85713958740234
    },
    {
      "id": "staar_result_compilation_0_9_5",
      "in": [
        {
          "id": "in_staar_results",
          "linkMerge": "merge_flattened",
          "source": [
            "staar_association_test_0_9_5/assoc_test_results"
          ]
        }
      ],
      "out": [
        {
          "id": "compiled_results"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.2",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-result-compilation-0-9-5/7",
        "baseCommand": [
          "bash",
          "compile.sh"
        ],
        "inputs": [
          {
            "id": "in_staar_results",
            "type": "File[]",
            "inputBinding": {
              "shellQuote": false,
              "position": 0
            },
            "label": "STAAR results",
            "doc": "STAAR results compressed files.",
            "sbg:fileTypes": "GZ"
          }
        ],
        "outputs": [
          {
            "id": "compiled_results",
            "doc": "Compiled STAAR results from multiple STAAR results files.",
            "label": "Compiled STAAR results",
            "type": "File?",
            "outputBinding": {
              "glob": "*_results.txt",
              "outputEval": "$(inheritMetadata(self, inputs.in_staar_results))"
            },
            "sbg:fileTypes": "TXT"
          }
        ],
        "doc": "**STAAR Result Compilation** compiles all the results from **STAAR Association Testing** into one file.\n\n*A list of **all inputs and parameters** with corresponding descriptions can be found at the bottom of the page.*\n\n***Please note that any cloud infrastructure costs resulting from app and pipeline executions, including the use of public apps, are the sole responsibility of you as a user. To avoid excessive costs, please read the app description carefully and set the app parameters and execution settings accordingly.***\n\n### Common use cases\n\nMerging results of **STAAR Association Testing** run on multiple genotype files.\n\n\n### Changes Introduced by Seven Bridges\n\nThe original code from [WDL implementation](https://github.com/sheilagaynor/STAAR_workflow/blob/main/STAAR_analysis_workflow.wdl) was adapted for platform execution, but it creates the same output as the original code.\n\n### Common issues and important notes\n\nThere are no known issues.\n\n### Performance Benchmarking\n\nIn the following table you can find estimates of running time and cost. \n      \n\n| Sample Count | Total sample size (GB) | Duration  | Cost - spot ($) |  Instance (AWS)  |\n|-------------------|-------------------|------------------|----------|-------------|------------|------------|\n|  |                     |                  |  |   |\n|  |         |                  |     |    |\n|  |                     |                  |      |  |\n\n### References\n[1] [STAAR GitHub](https://github.com/xihaoli/STAAR)\\\n[2] [STAAR_workflow GitHub](https://github.com/sheilagaynor/STAAR_workflow)\\\n[3] [STAAR manual](https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf)",
        "label": "STAAR Result Compilation",
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "ResourceRequirement",
            "ramMin": 16000,
            "coresMin": 8
          },
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/lea_lenhardt_ackovic/staar-0-9-5:1"
          },
          {
            "class": "InitialWorkDirRequirement",
            "listing": [
              {
                "entryname": "compile.sh",
                "entry": "${\n    return  '#!/bin/bash\\n{\\n  gunzip -c $1; shift\\n  for file in $@; do gunzip -c $file | sed \\'1d\\'; done\\n  } > compiled_results.txt'\n}",
                "writable": false
              }
            ]
          },
          {
            "class": "InlineJavascriptRequirement",
            "expressionLib": [
              "\nvar setMetadata = function(file, metadata) {\n    if (!('metadata' in file)) {\n        file['metadata'] = {}\n    }\n    for (var key in metadata) {\n        file['metadata'][key] = metadata[key];\n    }\n    return file\n};\nvar inheritMetadata = function(o1, o2) {\n    var commonMetadata = {};\n    if (!o2) {\n        return o1;\n    };\n    if (!Array.isArray(o2)) {\n        o2 = [o2]\n    }\n    for (var i = 0; i < o2.length; i++) {\n        var example = o2[i]['metadata'];\n        for (var key in example) {\n            if (i == 0)\n                commonMetadata[key] = example[key];\n            else {\n                if (!(commonMetadata[key] == example[key])) {\n                    delete commonMetadata[key]\n                }\n            }\n        }\n        for (var key in commonMetadata) {\n            if (!(key in example)) {\n                delete commonMetadata[key]\n            }\n        }\n    }\n    if (!Array.isArray(o1)) {\n        o1 = setMetadata(o1, commonMetadata)\n        if (o1.secondaryFiles) {\n            o1.secondaryFiles = inheritMetadata(o1.secondaryFiles, o2)\n        }\n    } else {\n        for (var i = 0; i < o1.length; i++) {\n            o1[i] = setMetadata(o1[i], commonMetadata)\n            if (o1[i].secondaryFiles) {\n                o1[i].secondaryFiles = inheritMetadata(o1[i].secondaryFiles, o2)\n            }\n        }\n    }\n    return o1;\n};"
            ]
          }
        ],
        "sbg:projectName": "STAAR 0.9.5 - Dev",
        "sbg:revisionsInfo": [
          {
            "sbg:revision": 0,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630318522,
            "sbg:revisionNotes": null
          },
          {
            "sbg:revision": 1,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630321011,
            "sbg:revisionNotes": "First draft."
          },
          {
            "sbg:revision": 2,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630321657,
            "sbg:revisionNotes": "Testing."
          },
          {
            "sbg:revision": 3,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630322497,
            "sbg:revisionNotes": "#!/bin/bash added as base command."
          },
          {
            "sbg:revision": 4,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630323628,
            "sbg:revisionNotes": "Still testing."
          },
          {
            "sbg:revision": 5,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630326726,
            "sbg:revisionNotes": "Changed bash code."
          },
          {
            "sbg:revision": 6,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630327498,
            "sbg:revisionNotes": "Changed file types for input and output."
          },
          {
            "sbg:revision": 7,
            "sbg:modifiedBy": "lea_lenhardt_ackovic",
            "sbg:modifiedOn": 1630481107,
            "sbg:revisionNotes": "Edited description."
          }
        ],
        "sbg:image_url": null,
        "sbg:toolkit": "STAAR",
        "sbg:toolkitVersion": "0.9.5",
        "sbg:wrapperAuthor": "Sheila Gaynor",
        "sbg:wrapperLicense": "MIT License",
        "sbg:toolAuthor": "Xihao Li, Zilin Li, Han Chen",
        "sbg:license": "GPLv3, MIT License",
        "sbg:categories": [
          "GWAS",
          "CWL1.2"
        ],
        "sbg:links": [
          {
            "id": "https://github.com/xihaoli/STAAR",
            "label": "STAAR Homepage"
          },
          {
            "id": "https://github.com/sheilagaynor/STAAR_workflow",
            "label": "STAAR workflow Homepage"
          },
          {
            "id": "https://github.com/xihaoli/STAAR/releases/tag/v0.9.5",
            "label": "Source Code"
          },
          {
            "id": "https://github.com/xihaoli/STAAR/archive/refs/tags/v0.9.5.tar.gz",
            "label": "Download"
          },
          {
            "id": "https://www.nature.com/articles/s41588-020-0676-4",
            "label": "Publication"
          },
          {
            "id": "https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf",
            "label": "Documentation"
          },
          {
            "id": "https://github.com/sheilagaynor/STAAR_workflow/blob/main/STAAR_analysis_workflow.wdl",
            "label": "Result Compilation Code"
          }
        ],
        "sbg:appVersion": [
          "v1.2"
        ],
        "sbg:id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-result-compilation-0-9-5/7",
        "sbg:revision": 7,
        "sbg:revisionNotes": "Edited description.",
        "sbg:modifiedOn": 1630481107,
        "sbg:modifiedBy": "lea_lenhardt_ackovic",
        "sbg:createdOn": 1630318522,
        "sbg:createdBy": "lea_lenhardt_ackovic",
        "sbg:project": "lea_lenhardt_ackovic/staar-0-9-5-dev",
        "sbg:sbgMaintained": false,
        "sbg:validationErrors": [],
        "sbg:contributors": [
          "lea_lenhardt_ackovic"
        ],
        "sbg:latestRevision": 7,
        "sbg:publisher": "sbg",
        "sbg:content_hash": "a7b2c9729bdd4f88b36373aaa7f567dc20194ab497f76f89572a6cdb65af107f3"
      },
      "label": "STAAR Result Compilation",
      "sbg:x": 201.16307067871094,
      "sbg:y": 311.8585510253906
    }
  ],
  "requirements": [
    {
      "class": "ScatterFeatureRequirement"
    },
    {
      "class": "InlineJavascriptRequirement"
    },
    {
      "class": "StepInputExpressionRequirement"
    }
  ],
  "sbg:projectName": "STAAR 0.9.5 - Dev",
  "sbg:revisionsInfo": [
    {
      "sbg:revision": 0,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1601533859,
      "sbg:revisionNotes": null
    },
    {
      "sbg:revision": 1,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1601537297,
      "sbg:revisionNotes": "First draft"
    },
    {
      "sbg:revision": 2,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1601544415,
      "sbg:revisionNotes": "Testing dot product Scatter"
    },
    {
      "sbg:revision": 3,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1601544479,
      "sbg:revisionNotes": "Testing scatter"
    },
    {
      "sbg:revision": 4,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1601546638,
      "sbg:revisionNotes": "Mem per job exposed"
    },
    {
      "sbg:revision": 5,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1603275827,
      "sbg:revisionNotes": ""
    },
    {
      "sbg:revision": 6,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1603276886,
      "sbg:revisionNotes": "Updated for conditional analysis"
    },
    {
      "sbg:revision": 7,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1603715520,
      "sbg:revisionNotes": "License updated"
    },
    {
      "sbg:revision": 8,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1603715659,
      "sbg:revisionNotes": "Updated licenses"
    },
    {
      "sbg:revision": 9,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1603794700,
      "sbg:revisionNotes": "Edited input formats"
    },
    {
      "sbg:revision": 10,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1603803228,
      "sbg:revisionNotes": "Edited input format"
    },
    {
      "sbg:revision": 11,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1612960781,
      "sbg:revisionNotes": "Updated tools in WF."
    },
    {
      "sbg:revision": 12,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1613650651,
      "sbg:revisionNotes": "Updated R scripts for both tools."
    },
    {
      "sbg:revision": 13,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1614072512,
      "sbg:revisionNotes": "Updated Null model and Association Testing to the latest version from GitHub."
    },
    {
      "sbg:revision": 14,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1614076643,
      "sbg:revisionNotes": "Docker images updated for both tools."
    },
    {
      "sbg:revision": 15,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1623147023,
      "sbg:revisionNotes": "STAAR script updated to the latest commit."
    },
    {
      "sbg:revision": 16,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630327449,
      "sbg:revisionNotes": "Added Result compilation and other tools updated."
    },
    {
      "sbg:revision": 17,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630327548,
      "sbg:revisionNotes": "Updated file types for result compilation step."
    },
    {
      "sbg:revision": 18,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630328804,
      "sbg:revisionNotes": "Scatter only with geno file."
    },
    {
      "sbg:revision": 19,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630394495,
      "sbg:revisionNotes": "STAAR null model changed."
    },
    {
      "sbg:revision": 20,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630483800,
      "sbg:revisionNotes": "Updated tools and description."
    },
    {
      "sbg:revision": 21,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630579352,
      "sbg:revisionNotes": "Updated description."
    },
    {
      "sbg:revision": 22,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630666186,
      "sbg:revisionNotes": "Updated benchmark info and description."
    },
    {
      "sbg:revision": 23,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1630676101,
      "sbg:revisionNotes": "Code and tools updated for cwltool execution. Portability added to description."
    },
    {
      "sbg:revision": 24,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1633094501,
      "sbg:revisionNotes": "Updated Assoc analysis step."
    },
    {
      "sbg:revision": 25,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1633512128,
      "sbg:revisionNotes": "Updated STAAR Association Testing step."
    },
    {
      "sbg:revision": 26,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1634635466,
      "sbg:revisionNotes": "STAAR Associat Testing updated to the latest revision."
    },
    {
      "sbg:revision": 27,
      "sbg:modifiedBy": "lea_lenhardt_ackovic",
      "sbg:modifiedOn": 1639727859,
      "sbg:revisionNotes": "Updated description."
    }
  ],
  "sbg:image_url": "https://platform.sb.biodatacatalyst.nhlbi.nih.gov/ns/brood/images/lea_lenhardt_ackovic/staar-0-9-5-dev/staar-rare-variant-workflow-0-9-5/27.png",
  "sbg:license": "GPLv3, MIT License",
  "sbg:links": [
    {
      "id": "https://github.com/xihaoli/STAAR",
      "label": "STAAR Homepage"
    },
    {
      "id": "https://github.com/sheilagaynor/STAAR_workflow",
      "label": "STAAR workflow Homepage"
    },
    {
      "id": "https://github.com/xihaoli/STAAR/releases/tag/v0.9.5",
      "label": "Source Code"
    },
    {
      "id": "https://github.com/xihaoli/STAAR/archive/refs/tags/v0.9.5.tar.gz",
      "label": "Download"
    },
    {
      "id": "https://www.nature.com/articles/s41588-020-0676-4",
      "label": "Publication"
    },
    {
      "id": "https://content.sph.harvard.edu/xlin/dat/STAAR-manual.pdf",
      "label": "Documentation"
    },
    {
      "id": "https://github.com/sheilagaynor/STAAR_workflow/blob/main/STAAR_analysis_workflow.wdl",
      "label": "STAAR workflow Code"
    }
  ],
  "sbg:wrapperLicense": "",
  "sbg:toolAuthor": "Sheila Gaynor",
  "sbg:categories": [
    "GWAS",
    "CWL1.2"
  ],
  "sbg:appVersion": [
    "v1.2"
  ],
  "id": "https://api.sb.biodatacatalyst.nhlbi.nih.gov/v2/apps/lea_lenhardt_ackovic/staar-0-9-5-dev/staar-rare-variant-workflow-0-9-5/27/raw/",
  "sbg:id": "lea_lenhardt_ackovic/staar-0-9-5-dev/staar-rare-variant-workflow-0-9-5/27",
  "sbg:revision": 27,
  "sbg:revisionNotes": "Updated description.",
  "sbg:modifiedOn": 1639727859,
  "sbg:modifiedBy": "lea_lenhardt_ackovic",
  "sbg:createdOn": 1601533859,
  "sbg:createdBy": "lea_lenhardt_ackovic",
  "sbg:project": "lea_lenhardt_ackovic/staar-0-9-5-dev",
  "sbg:sbgMaintained": false,
  "sbg:validationErrors": [],
  "sbg:contributors": [
    "lea_lenhardt_ackovic"
  ],
  "sbg:latestRevision": 27,
  "sbg:publisher": "sbg",
  "sbg:content_hash": "adb004cf133326d15c19d6ac25500d3826df7d94253ced5d63e8cc2f2ff634653"
}