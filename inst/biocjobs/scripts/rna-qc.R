params <- BiocJobs::jobParams("scrapper", "rna-qc")

suppressPackageStartupMessages({
  library(LoomExperiment)
})

sce <- import(
  params$infile,
  format = "loom",
  type = "SingleCellLoomExperiment"
)

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

export(
  sce,
  params$outfile,
  format = "loom"
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