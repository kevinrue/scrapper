## BiocJobs job script: cluster-graph-se
## The declared interface lives in ../cluster-graph-se.yaml.

## Test command (R, run from the repository root; overwrites the test output)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/cluster-graph-se.yaml"), params = list(infile = "test-data/sce.pca.h5ad", outfile = "test-data/sce.clustered.h5ad", neighbours = 3, method = "multilevel", resolution = 1, plot_clusters = "test-data/clusters.pdf", plot_reduced_dim = "PCA", plot_width = 7, plot_height = 5, point_size = 1))

## Validation command (Bash)
## Rscript -e 'BiocJobs::biocjobsCLI()' validate .

# Rscript -e 'BiocJobs::biocjobsCLI()' galaxy   . cluster-graph-se --out wrappers/cluster-graph-se.xml
# Rscript -e 'BiocJobs::biocjobsCLI()' nextflow . cluster-graph-se --out wrappers/cluster-graph-se.nf

params <- BiocJobs::jobParams("scrapper", "cluster-graph-se")

suppressPackageStartupMessages({
  library(anndataR)
  library(ggplot2)
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

cluster_data <- SingleCellExperiment::reducedDim(sce, params$plot_reduced_dim)
if (ncol(cluster_data) < 2L) {
  stop("The reducedDim entry selected for plotting must contain at least two dimensions.")
}
cluster_data <- data.frame(
  dimension1 = cluster_data[, 1],
  dimension2 = cluster_data[, 2],
  clusters = factor(SummarizedExperiment::colData(sce)$clusters)
)
cluster_plot <- ggplot2::ggplot(
  data = cluster_data,
  mapping = ggplot2::aes(x = dimension1, y = dimension2, colour = clusters)
) +
  ggplot2::geom_point(size = params$point_size) +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = paste(params$plot_reduced_dim, "1"),
    y = paste(params$plot_reduced_dim, "2"),
    colour = "Cluster"
  )

ggplot2::ggsave(
  filename = params$plot_clusters,
  plot = cluster_plot,
  width = params$plot_width,
  height = params$plot_height,
  device = "pdf"
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
