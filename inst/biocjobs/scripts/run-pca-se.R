## BiocJobs job script: run-pca-se
## The declared interface lives in ../run-pca-se.yaml; jobParams() parses the
## command line against it, so every value below is typed, validated and
## defaulted by the time it returns.

# Test command (R)
## BiocJobs::runJob(BiocJobs::readJob("inst/biocjobs/run-pca-se.yaml"), params = list(infile = "test-data/sce.normalized.h5ad", features = "test-data/hvgs.txt", n_pcs = 25, plot_var_width = 7, plot_var_height = 5, outfile = "sce.pca.h5ad", plot_var_explained = "pca_var_explained.pdf"))

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
  library(ggplot2)
  library(scrapper)
  library(SingleCellExperiment)
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

pc_names <- colnames(SingleCellExperiment::reducedDim(sce, "PCA"))
var_explained <- metadata(sce)$PCA$variance.explained

var_explained_data <- data.frame(
  PC = factor(pc_names, levels = pc_names),
  variance_explained = var_explained
)
pca_var_plot <- ggplot2::ggplot(
  data = var_explained_data,
  mapping = ggplot2::aes(
    x = PC,
    y = variance_explained
  )
) +
  ggplot2::geom_col() +
  ggplot2::theme_classic() +
  ggplot2::labs(
    x = "Principal component",
    y = "Variance explained"
  )

ggplot2::ggsave(
  filename = params$plot_var_explained,
  plot = pca_var_plot,
  width = params$plot_var_width,
  height = params$plot_var_height
)

## Provenance to the job log.
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
