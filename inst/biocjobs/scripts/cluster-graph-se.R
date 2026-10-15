## BiocJobs job script: cluster-graph-se
## The declared interface lives in ../cluster-graph-se.yaml.

## Test command (R, run from the repository root; overwrites the test output)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/cluster-graph-se.yaml"), params = list(infile = "test-data/sce.pca.h5ad", outfile = "test-data/sce.clustered.h5ad", neighbours = 3, method = "multilevel", resolution = 1))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . cluster-graph-se --out wrappers/cluster-graph-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . cluster-graph-se --out wrappers/cluster-graph-se.nf

params <- BiocJobs::jobParams("scrapper", "cluster-graph-se")

suppressPackageStartupMessages({
  library(anndataR)
  library(scrapper)
  library(SingleCellExperiment)
})

sce <- anndataR::read_h5ad(params$infile, as = "SingleCellExperiment")

sce <- clusterGraph.se(
  x = sce,
  num.neighbors = params$neighbours,
  method = params$method,
  resolution = params$resolution
)

anndataR::write_h5ad(
  object = sce,
  compression = "gzip",
  path = params$outfile,
  mode = "w"
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
