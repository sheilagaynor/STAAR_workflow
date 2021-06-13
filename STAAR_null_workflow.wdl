workflow STAAR_null {
    File pheno_file
    String null_file
    String sample_name
    String outcome_name
    String? outcome_type = "continuous"
    String? covariate_names = "NA"
    File? kinship_file = "NA"
    String? het_var_name = "NA"
    Int? null_memory = 25
    Int? null_disk = 50

    call run_null_model {
            input:
                pheno_file = pheno_file,
                null_file = null_file,
                sample_name = sample_name,
                outcome_name = outcome_name,
                outcome_type = outcome_type,
                covariate_names = covariate_names,
                kinship_file = kinship_file,
                het_var_name = het_var_name,
                null_memory = null_memory,
                null_disk = null_disk
    }

    output {
        File null_model = run_null_model.null_model
    }
}

task run_null_model {
    File pheno_file
    String null_file
    String sample_name
    String outcome_name
    String outcome_type
    String covariate_names
    File kinship_file
    String het_var_name
    Int null_memory
    Int null_disk

    command {
        Rscript /staar_workflow/STAAR_null_model.R ${pheno_file} ${null_file} ${sample_name} ${outcome_name} ${outcome_type} ${covariate_names} ${kinship_file} ${het_var_name}
    }

    runtime {
        docker: "quay.io/sheilagaynor/staar_workflow"
        memory: "${null_memory} GB"
        disks: "local-disk ${null_disk} HDD"
    }

    output {
        File null_model = select_first(glob("${null_file}*"))
    }
}
