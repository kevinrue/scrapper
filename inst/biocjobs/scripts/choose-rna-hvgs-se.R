## BiocJobs job script: choose-rna-hvgs-se
## The declared interface lives in ../choose-rna-hvgs-se.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

## Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/choose-rna-hvgs-se.yaml"), params = list(infile = "test-data/sce.normalized.h5ad", outfile = "hvgs.txt"))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . choose-rna-hvgs-se --out wrappers/choose-rna-hvgs-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . choose-rna-hvgs-se --out wrappers/choose-rna-hvgs-se.nf

params <- BiocJobs::jobParams("scrapper", "choose-rna-hvgs-se")

## -- your analysis ---------------------------------------------------------
## Read inputs from params$<input name>, honour params$<option name>, and
## write every declared output to params$<output name>. Fail loudly with
## stop() on bad input; never prompt or open devices interactively.

suppressPackageStartupMessages({
  library(anndataR)
  library(scrapper)
})

sce <- anndataR::read_h5ad("test-data/sce.normalized.h5ad", as = "SingleCellExperiment")

sce <- chooseRnaHvgs.se(sce)

hvg_rownames <- rownames(sce)[SummarizedExperiment::rowData(sce)$hvg]

write(x = hvg_rownames, file = params$hvgs)

# TODO: add QC plots outputs

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
