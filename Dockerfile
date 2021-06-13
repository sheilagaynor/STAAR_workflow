# Set base image
FROM uwgac/r-3.6.3-mkl:3.6.3@sha256:b50536df2e6b1ce1d11261e422aaf08e86f6213bc660b0fd45ed61dc1eb7b4ee

# Prepare image
RUN sudo apt-get update && sudo apt-get -y install git
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

# Install dependencies
RUN Rscript -e 'install.packages(c("BiocManager","Rcpp","Matrix","RcppArmadillo","readr","data.table","dplyr","doMC","R.utils"))'
RUN Rscript -e 'BiocManager::install(c("SeqArray","gdsfmt","SeqVarTools","foreach","GMMAT","CompQuadForm","GENESIS","TxDb.Hsapiens.UCSC.hg38.knownGene"))'

# Install STAAR v0.9.5 R package from source
COPY STAAR_0.9.5.tar.gz /STAAR_0.9.5.tar.gz
RUN Rscript -e 'install.packages("STAAR_0.9.5.tar.gz", repos=NULL, type="source")'

# Copy in R scripts
COPY STAAR_analysis.R STAAR_null_model.R /
