## BiocJobs job script: normalize-rna-counts-se
## The declared interface lives in ../normalize-rna-counts-se.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

# Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/normalize-rna-counts-se.yaml"), params = list(infile = "test-data/sce.rna-qc.h5ad", outfile = "sce.normalized.h5ad"))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . normalize-rna-counts-se --out wrappers/normalize-rna-counts-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . normalize-rna-counts-se --out wrappers/normalize-rna-counts-se.nf

params <- BiocJobs::jobParams("scrapper", "normalize-rna-counts-se")

## -- your analysis ---------------------------------------------------------
## Read inputs from params$<input name>, honour params$<option name>, and
## write every declared output to params$<output name>. Fail loudly with
## stop() on bad input; never prompt or open devices interactively.

suppressPackageStartupMessages({
  library(anndataR)
  library(scrapper)
})

sce <- anndataR::read_h5ad(params$infile, as = "SingleCellExperiment")

sce <- scrapper::normalizeRnaCounts.se(sce)

anndataR::write_h5ad(
  object = sce,
  compression = "gzip",
  path = params$outfile,
  mode = "w"
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
