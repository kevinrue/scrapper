## BiocJobs job script: run-umap-se
## The declared interface lives in ../run-umap-se.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

# Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/run-umap-se.yaml"), params = list(infile = "test-data/sce.pca.h5ad", neighbours = 3, outfile = "sce.umap.h5ad", plot_umap = "umap.pdf", plot_width = 7, plot_height = 5))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . run-umap-se --out wrappers/run-umap-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . run-umap-se --out wrappers/run-umap-se.nf

params <- BiocJobs::jobParams("scrapper", "run-umap-se")

suppressPackageStartupMessages({
  library(anndataR)
  library(ggplot2)
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

umap_data <- as.data.frame(SingleCellExperiment::reducedDim(sce, "UMAP"))
umap_plot <- ggplot2::ggplot(
  data = umap_data,
  mapping = ggplot2::aes(x = UMAP1, y = UMAP2)
) +
  ggplot2::geom_point() +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "UMAP 1", y = "UMAP 2")

ggplot2::ggsave(
  filename = params$plot_umap,
  plot = umap_plot,
  width = params$plot_width,
  height = params$plot_height,
  device = "pdf"
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
