## BiocJobs job script: rna-qc
##
## The declared interface lives in ../rna-qc.yaml.
## jobParams() parses the command line against that declaration, so by the
## time it returns, every value below is typed, validated and defaulted.

## Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/rna-qc.yaml"), params = list(infile = "test-data/sce.h5ad", mitochondrial_prefix = "MT-", num_mads = 3, outfile = "sce.rna-qc.h5ad"))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . rna-qc --out wrappers/rna-qc.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . rna-qc --out wrappers/rna-qc.nf

params <- BiocJobs::jobParams("scrapper", "rna-qc")

suppressPackageStartupMessages({
  library(anndataR)
  library(scrapper)
  library(SummarizedExperiment)
})

sce <- anndataR::read_h5ad(params$infile, as = "SingleCellExperiment")

gene_symbols <- SummarizedExperiment::rowData(sce)$Symbol

is_mito <- startsWith(
  gene_symbols,
  params$mitochondrial_prefix
)
is_mito[is.na(is_mito)] <- FALSE

if (any(is_mito)) {
  qc_subsets <- list(MT = is_mito)
} else {
  qc_subsets <- list()
}

sce <- scrapper::quickRnaQc.se(
  sce,
  subsets = qc_subsets,
  more.suggest.args = list(
    num.mads = params$num_mads
  )
)

anndataR::write_h5ad(
  object = sce,
  compression = "gzip",
  path = params$outfile
)

qc <- as.data.frame(
  SummarizedExperiment::colData(sce)
)

write.table(
  qc,
  params$qc_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
