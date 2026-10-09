# BiocJobs entry point for scrapper::scoreMarkers.se().
# Interface: ../score-markers-se.yaml
# Sources checked on 2026-09-16:
# https://github.com/libscran/scrapper/blob/master/R/se_scoreMarkers.R
# https://github.com/libscran/scrapper/blob/master/R/scoreMarkers.R
# https://github.com/almahmoud/BiocJobs/blob/main/docs/developer-guide.md

params <- BiocJobs::jobParams("scrapper", "score-markers-se")

# Register the relevant S4 classes before deserializing the input.
for (pkg in c("scrapper", "SummarizedExperiment", "SingleCellExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("Required package is not installed: ", pkg, call. = FALSE)
    }
}
if (!"scoreMarkers.se" %in% getNamespaceExports("scrapper")) {
    stop("The installed scrapper does not export scoreMarkers.se(); ",
         "install a version containing this function.", call. = FALSE)
}

# Do not accidentally replace the input with the output.
if (identical(normalizePath(params$se, winslash = "/", mustWork = TRUE),
              normalizePath(params$markers, winslash = "/", mustWork = FALSE))) {
    stop("--se and --markers must refer to different files.", call. = FALSE)
}

x <- tryCatch(
    readRDS(params$se),
    error = function(e) {
        stop("Cannot read --se as RDS: ", conditionMessage(e), call. = FALSE)
    }
)
if (!methods::is(x, "SummarizedExperiment")) {
    stop("--se must contain a SummarizedExperiment or a subclass, ",
         "such as SingleCellExperiment.", call. = FALSE)
}
if (nrow(x) == 0L || ncol(x) == 0L) {
    stop("--se must contain at least one gene and at least one cell.",
         call. = FALSE)
}

# Resolve a name first, then a positive one-based index.
assay_names <- SummarizedExperiment::assayNames(x)
n_assays <- length(SummarizedExperiment::assays(x, withDimnames = FALSE))
assay_type <- params$assay_type
if (assay_type %in% assay_names) {
    if (sum(assay_names == assay_type, na.rm = TRUE) != 1L) {
        stop("--assay_type matches multiple assay names; use an index.",
             call. = FALSE)
    }
} else if (grepl("^[1-9][0-9]*$", assay_type)) {
    assay_index <- suppressWarnings(as.numeric(assay_type))
    if (!is.finite(assay_index) || assay_index > n_assays) {
        stop("--assay_type index exceeds the number of assays: ", n_assays,
             call. = FALSE)
    }
    assay_type <- as.integer(assay_index)
} else {
    stop("--assay_type does not identify an assay. Available names: ",
         paste(assay_names, collapse = ", "),
         ". Alternatively, supply a one-based index.", call. = FALSE)
}

# A reduction checks every value without explicitly coercing the whole
# assay to a dense matrix. This also detects inaccessible backing files.
assay_range <- tryCatch(
    range(SummarizedExperiment::assay(x, assay_type)),
    error = function(e) {
        stop("Cannot inspect the selected assay: ", conditionMessage(e),
             ". File-backed assays require accessible backing files.",
             call. = FALSE)
    }
)
if (!is.numeric(assay_range) || length(assay_range) != 2L) {
    stop("The selected assay must contain numeric expression values.",
         call. = FALSE)
}
if (any(!is.finite(assay_range))) {
    stop("The selected assay contains NA, NaN, or infinite values; ",
         "resolve these before scoring markers.", call. = FALSE)
}

# Metadata is already aligned with cells by the SummarizedExperiment.
cell_data <- SummarizedExperiment::colData(x)
read_assignment <- function(column, option) {
    if (!nzchar(trimws(column)) || !column %in% colnames(cell_data)) {
        stop("--", option, " must name an existing colData() column: ",
             column, call. = FALSE)
    }
    if (sum(colnames(cell_data) == column, na.rm = TRUE) != 1L) {
        stop("--", option, " matches duplicated colData() column names.",
             call. = FALSE)
    }
    values <- cell_data[[column]]
    if (methods::is(values, "Rle")) {
        values <- as.vector(values)
    }
    if (!is.atomic(values) || !is.null(dim(values)) ||
        length(values) != ncol(x)) {
        stop("--", option, " must contain one atomic assignment per cell.",
             call. = FALSE)
    }
    if (anyNA(values) || any(!nzchar(trimws(as.character(values))))) {
        stop("--", option, " contains missing or blank assignments.",
             call. = FALSE)
    }
    if (is.numeric(values) && any(!is.finite(values))) {
        stop("--", option, " contains non-finite assignments.", call. = FALSE)
    }
    # Preserve an existing factor's order; discard unobserved levels.
    if (is.factor(values)) droplevels(values) else factor(values)
}

groups <- read_assignment(params$groups_column, "groups_column")
if (nlevels(groups) < 2L) {
    stop("--groups_column must contain at least two observed groups.",
         call. = FALSE)
}
block <- NULL
if (nzchar(params$block_column)) {
    block <- read_assignment(params$block_column, "block_column")
}

# Parse vector-valued options as data, never as R expressions.
split_csv <- function(value, option) {
    if (!nzchar(trimws(value))) {
        return(character())
    }
    pieces <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
    if (any(!nzchar(pieces)) || grepl(",[[:space:]]*$", value)) {
        stop("--", option, " contains an empty comma-separated entry.",
             call. = FALSE)
    }
    pieces
}

extra_columns <- split_csv(params$extra_columns, "extra_columns")
if (length(extra_columns)) {
    row_names <- colnames(SummarizedExperiment::rowData(x))
    if (anyDuplicated(extra_columns)) {
        stop("--extra_columns contains duplicate names.", call. = FALSE)
    }
    absent <- setdiff(extra_columns, row_names)
    if (length(absent)) {
        stop("--extra_columns not found in rowData(): ",
             paste(absent, collapse = ", "), call. = FALSE)
    }
    if (any(vapply(extra_columns, function(nm) {
        sum(row_names == nm, na.rm = TRUE) != 1L
    }, logical(1)))) {
        stop("--extra_columns matches duplicated rowData() names.",
             call. = FALSE)
    }
    reserved <- extra_columns %in% c("mean", "detected") |
        grepl("^(cohens\\.d|auc|delta\\.mean|delta\\.detected)\\.",
              extra_columns)
    if (any(reserved)) {
        stop("Rename annotations that conflict with marker statistics: ",
             paste(extra_columns[reserved], collapse = ", "), call. = FALSE)
    }
} else {
    extra_columns <- NULL
}

quantile_text <- split_csv(params$summary_quantiles, "summary_quantiles")
summary_quantiles <- NULL
if (length(quantile_text)) {
    summary_quantiles <- suppressWarnings(as.numeric(quantile_text))
    if (any(!is.finite(summary_quantiles)) ||
        any(summary_quantiles < 0 | summary_quantiles > 1)) {
        stop("--summary_quantiles must contain finite numbers in [0,1].",
             call. = FALSE)
    }
    summary_quantiles <- sort(unique(summary_quantiles))
}

if (!is.finite(params$threshold)) {
    stop("--threshold must be finite.", call. = FALSE)
}
if (!is.finite(params$block_quantile)) {
    stop("--block_quantile must be finite.", call. = FALSE)
}

# scoreMarkers.se() has no '...' argument. Underlying scoreMarkers()
# settings belong in more.marker.args, not at the top level of the call.
# Its formatter needs summary-mode results, so all.pairwise stays FALSE.
more_marker_args <- list(
    all.pairwise = FALSE,
    threshold = params$threshold,
    compute.cohens.d = params$compute_cohens_d,
    compute.auc = params$compute_auc,
    compute.delta.mean = params$compute_delta_mean,
    compute.delta.detected = params$compute_delta_detected
)
if (!is.null(summary_quantiles)) {
    more_marker_args$compute.summary.quantiles <- summary_quantiles
}
if (params$min_rank_limit > 0L) {
    more_marker_args$min.rank.limit <- params$min_rank_limit
}
if (!is.null(block)) {
    if (params$block_average_policy != "default") {
        more_marker_args$block.average.policy <- params$block_average_policy
    }
    if (params$block_weight_policy != "default") {
        more_marker_args$block.weight.policy <- params$block_weight_policy
    }
    more_marker_args$block.quantile <- params$block_quantile
}

order_by <- params$order_by
if (!nzchar(trimws(order_by))) {
    stop("--order_by must be auto, none, or an output column name.",
         call. = FALSE)
}
if (identical(order_by, "auto")) {
    order_by <- TRUE
} else if (identical(order_by, "none")) {
    order_by <- FALSE
}

message("scrapper ", utils::packageVersion("scrapper"), ": scoring ",
        nrow(x), " genes across ", ncol(x), " cells in ", nlevels(groups),
        " groups; threads=", params$num_threads)
if (!is.null(block)) {
    message("Blocking enabled: ", nlevels(block), " observed blocks.")
}

markers <- scrapper::scoreMarkers.se(
    x = x,
    groups = groups,
    block = block,
    num.threads = params$num_threads,
    more.marker.args = more_marker_args,
    assay.type = assay_type,
    extra.columns = extra_columns,
    order.by = order_by
)

# In the inspected upstream formatter, an unknown order.by column can
# silently yield empty tables. Reject it before writing a result file.
if (!methods::is(markers, "List") || length(markers) != nlevels(groups)) {
    stop("Unexpected scoreMarkers.se() result structure.", call. = FALSE)
}
if (is.character(order_by) && !order_by %in% colnames(markers[[1L]])) {
    stop("--order_by column was not produced: ", order_by,
         ". Available columns: ", paste(colnames(markers[[1L]]), collapse = ", "),
         ". It may refer to a disabled effect size.", call. = FALSE)
}
rows_per_group <- vapply(seq_along(markers), function(i) {
    nrow(markers[[i]])
}, integer(1))
if (any(rows_per_group != nrow(x))) {
    stop("Unexpected gene count in the result; check --order_by and input data.",
         call. = FALSE)
}

# RDS preserves the S4 result, including non-tabular rowData annotations.
saveRDS(markers, file = params$markers)
message("Saved marker statistics to ", params$markers)
message(paste(utils::capture.output(utils::sessionInfo()), collapse = "\n"))
