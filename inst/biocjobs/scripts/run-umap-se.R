## BiocJobs job script: run-umap-se
## The declared interface lives in ../run-umap-se.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

# Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/run-umap-se.yaml"), params = list(infile = "test-data/sce.pca.h5ad", neighbours = 3, outfile = "sce.umap.h5ad"))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . run-umap-se --out wrappers/run-umap-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . run-umap-se --out wrappers/run-umap-se.nf

params <- BiocJobs::jobParams("scrapper", "run-umap-se")

suppressPackageStartupMessages({
  library(anndataR)
  library(scrapper)
  library(SingleCellExperiment)
})

sce <- anndataR::read_h5ad(params$infile, as = "SingleCellExperiment")

sce <- runUmap.se(sce, num.neighbors = params$neighbours)

anndataR::write_h5ad(
  object = sce,
  compression = "gzip",
  path = params$outfile,
  mode = "w"
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
