workflow STAAR_analysis {

    # run_null_model inputs
    File? null_file_precompute
    File? pheno_file
    String? null_file_name
    String? sample_name
    String? outcome_name
    String? outcome_type = "continuous"
    String? covariate_names = "NA"
    File? kinship_file
    String? het_var_name = "NA"
    Int? null_memory = 25
    Int? null_disk = 50

    # run_analysis inputs
    Array[File] geno_files
    Array[File]? annot_files
    String results_file
    String? agds_file = "None"
    String? agds_annot_channels = "None"
    File? agg_file
    File? cond_file
    Array[File]? cond_geno_files
    File? cand_file
    Float? maf_thres = 0.05
    Int? mac_thres = 1
    Int? window_length = 2000
    Int? step_length = 1000
    # Compute inputs
    Int? num_cores = 3
    Int? num_chunk_divisions = 100000
    Int? test_memory = 25
    Int? test_disk = 50

    if (!defined(null_file_precompute)) {
        call run_null_model {
            input:
                pheno_file = pheno_file,
                null_file_name = null_file_name,
                sample_name = sample_name,
                outcome_name = outcome_name,
                outcome_type = outcome_type,
                covariate_names = covariate_names,
                kinship_file = kinship_file,
                het_var_name = het_var_name,
                null_memory = null_memory,
                null_disk = null_disk
        }
    }

    File null_file = select_first([null_file_precompute, run_null_model.null_model])

    Array[File] annot_avail = select_first([ annot_files, []])

    if (length(annot_avail) == length(geno_files)) {
        Array[Pair[File,File]] geno_annot_pairs = zip(geno_files, annot_avail)
        scatter (geno_annot_set in geno_annot_pairs) {
        File geno_in = geno_annot_set.left
        File annot_in = geno_annot_set.right
            call run_analysis {
                input:
                    null_file = null_file,
                    geno_file = geno_in,
                    annot_file = annot_in,
                    results_file = results_file,
                    agg_file = agg_file,
                    cond_file = cond_file,
                    cond_geno_files = cond_geno_files,
                    cand_file = cand_file,
                    maf_thres = maf_thres,
                    mac_thres = mac_thres,
                    window_length = window_length,
                    step_length = step_length,
                    num_cores = num_cores,
                    num_chunk_divisions = num_chunk_divisions,
                    test_memory = test_memory,
                    test_disk = test_disk
            }
        }
    }

    if (!defined(annot_files)) {
        scatter (geno_in in geno_files) {
            call run_analysis as run_analysis_annotfree {
                input:
                  null_file = null_file,
                  geno_file = geno_in,
                  results_file = results_file,
                  agds_file = agds_file,
                  agds_annot_channels = agds_annot_channels,
                  agg_file = agg_file,
                  cond_file = cond_file,
                  cond_geno_files = cond_geno_files,
                  cand_file = cand_file,
                  maf_thres = maf_thres,
                  mac_thres = mac_thres,
                  window_length = window_length,
                  step_length = step_length,
                  num_cores = num_cores,
                  num_chunk_divisions = num_chunk_divisions,
                  test_memory = test_memory,
                  test_disk = test_disk
            }
        }
    }

    call run_compilation {
		input:
			results_array = if (length(annot_avail) == length(geno_files)) then run_analysis.results else run_analysis_annotfree.results
	  }

    output {
        File null_model = null_file
        File staar_results = run_compilation.compiled_results
    }

    meta {
                author: "Sheila Gaynor"
                email: "sheilagaynor@hsph.harvard.edu"
                description: "Run rare variant analysis incorporating functional annotations using STAAR and return compiled results."
        }
}



task run_null_model {
    File pheno_file
    String null_file_name
    String sample_name
    String outcome_name
    String outcome_type
    String covariate_names
    File? kinship_file
    String het_var_name
    Int null_memory
    Int null_disk

    command {
        Rscript /STAAR_null_model.R ${pheno_file} ${null_file_name} ${sample_name} ${outcome_name} ${outcome_type} ${covariate_names} ${default="NA" kinship_file} ${het_var_name}
    }

    runtime {
        docker: "quay.io/sheilagaynor/staar_workflow"
        memory: "${null_memory} GB"
        disks: "local-disk ${null_disk} HDD"
    }

    output {
        File null_model = select_first(glob("${null_file_name}*"))
    }
}

task run_analysis {
    File null_file
    File geno_file
    File? annot_file
    String results_file
    String agds_file
    String agds_annot_channels
    File? agg_file
    File? cond_file
    Array[File]? cond_geno_files
    File? cand_file
    String maf_thres
    String mac_thres
    Int window_length
    Int step_length
    Int num_cores
    Int num_chunk_divisions
    Int test_memory
    Int test_disk

    command {
        Rscript /STAAR_analysis.R ${null_file} ${geno_file} ${default="None" annot_file} ${results_file} ${default="None" agds_file} ${default="None" agds_annot_channels} ${default="None" agg_file} ${default="None" cond_file} ${default="None" sep="," cond_geno_files} ${default="None" cand_file} ${maf_thres} ${mac_thres} ${window_length} ${step_length} ${num_cores} ${num_chunk_divisions}
    }

    runtime {
        docker: "quay.io/sheilagaynor/staar_workflow"
        preemptible: 1
        maxRetries: 3
        memory: "${test_memory} GB"
        disks: "local-disk ${test_disk} HDD"
    }

    output {
        File results = select_first(glob("*.gz"))
    }
}

task run_compilation {
  Array[File] results_array

  command <<<
  set -- results_array
  {
    gunzip -c ${results_array[0]}; shift
    for file in ${sep=" " results_array}; do gunzip 0c $file | sed '1d'; done
    } > compiled_results.txt
  >>>

    runtime {
      docker: "ubuntu:latest"
      disks: "local-disk 10 HDD"
    }

    output {
    File compiled_results = "compiled_results.txt"
  }
}
