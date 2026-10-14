## BiocJobs job script: run-pca-se
## The declared interface lives in ../run-pca-se.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

# Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/run-pca-se.yaml"), params = list(infile = "test-data/sce.normalized.h5ad", features = "test-data/hvgs.txt", n_pcs = 25, outfile = "sce.pca.h5ad"))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . run-pca-se --out wrappers/run-pca-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . run-pca-se --out wrappers/run-pca-se.nf

params <- BiocJobs::jobParams("scrapper", "run-pca-se")

## -- your analysis ---------------------------------------------------------
## Read inputs from params$<input name>, honour params$<option name>, and
## write every declared output to params$<output name>. Fail loudly with
## stop() on bad input; never prompt or open devices interactively.

suppressPackageStartupMessages({
  library(anndataR)
  library(scrapper)
})

sce <- anndataR::read_h5ad(params$infile, as = "SingleCellExperiment")

feature_rownames <- scan(params$features, what = "character", quiet = TRUE)

sce <- runPca.se(
  x = sce,
  features = feature_rownames,
  number = params$n_pcs
)

# anndataR::write_h5ad() cannot overwrite an existing output file
# [Galaxy]: this is necessary because Galaxy creates an empty output file during tests
if (file.exists(params$outfile)) {
  file.remove(params$outfile)
}

anndataR::write_h5ad(
  object = sce,
  compression = "gzip",
  path = params$outfile
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
