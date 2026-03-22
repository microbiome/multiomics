#' @name
#' getMantel
#'
#' @title
#' Mantel test between two tables
#'
#' @description
#' Performs a Mantel test to assess the association between two data matrices
#' (e.g., different omics layers) by comparing their dissimilarity matrices.
#' Mantel test can be used to test whether there is general association between
#' two omics layers.
#'
#' @details
#' The Mantel test evaluates the correlation between two distance matrices.
#' For multi-assay objects, the function:
#' \enumerate{
#'   \item extracts shared samples across selected experiments,
#'   \item computes dissimilarity matrices for each dataset,
#'   \item optionally applies stratification (group-wise permutations),
#'   \item performs the Mantel test using \code{vegan::mantel()}.
#' }
#'
#' @references
#' Borcard, D. & Legendre, P. (2012) Is the Mantel correlogram powerful enough
#' to be useful in ecological analysis? A simulation study. Ecology 93:
#' 1473-1481.
#'
#' @return
#' An object of class \code{"mantel"} as returned by
#' \code{vegan::mantel()}.
#'
#' @param x \code{MultiAssayExperiment} or \code{SingleCellExperiment}.
#'
#' @param experiments \code{Character vector} or \code{integer vector}. Names
#' or indices of experiments selected from \code{x}.
#'
#' @param assay.types \code{Character vector}. Names of assays selected from
#' \code{experiments}. Must match the length of
#' \code{experiments}.
#'
#' @param group.by \code{Character scalar}. Column name in sample metadata used
#' for stratified permutations. If \code{NULL}, non-stratified test is
#' performed. (Default: \code{NULL})
#'
#' @param dist.methods \code{Character vector}. Distance methods used to compute
#' dissimilarity matrices (one per dataset).  Must match the length of
#' \code{experiments}. Must be values from \code{vegan::vegdist()}.
#' (Default: \code{c("bray", "euclidean")})
#'
#' @param method \code{Character scalar}. Correlation method for Mantel test.
#' Options include \code{"pearson"}, \code{"spearman"}, or \code{"kendall"}.
#' (Default: \code{"kendall"})
#'
#' @param ... additional arguments:
#'  \itemize{
#'   \item \code{npermutations}: \code{Integer scalar}. Number of permutations
#'   used to assess significance. (Default: \code{999}).
#' }
#'
#' @examples
#' library(multiomics)
#' library(mia)
#'
#' data(HintikkaXOData)
#' mae <- HintikkaXOData
#' # Replace missing values with 0
#' assay(mae[[2]])[is.na(assay(mae[[2]]))] <- 0
#'
#' # Perform Mantel test
#' getMantel(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("counts", "nmr")
#' )
#'
#' @seealso
#' \code{\link[vegan:vegdist]{vegan::vegdist()}}
#'
NULL

#' @rdname getMantel
#' @export
setMethod("getMantel",
    signature = c(x = "MultiAssayExperiment"),
    function(x, experiments, assay.types, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        #
        mat_list <- .get_shared_samples_from_mae(
            x, experiments, assay.types, ...
        )
        group.by <- .get_grouping_for_mantel_mae(x, experiments, ...)
        args <- list(x = mat_list, group.by = group.by)
        args <- c(
            args, list(...)[!names(list(...)) %in% c("group.by", "strata")]
        )
        res <- do.call(getMantel, args)
        return(res)
    }
)

#' @rdname getMantel
#' @export
setMethod("getMantel",
    signature = c(x = "SingleCellExperiment"),
    function(x, experiments, assay.types, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        #
        mat_list <- .get_shared_samples_from_tse(
            x, experiments, assay.types, ...
        )
        group.by <- .get_grouping_for_mantel_tse(x, experiments, ...)
        args <- list(x = mat_list, group.by = group.by)
        args <- c(
            args, list(...)[!names(list(...)) %in% c("group.by", "strata")]
        )
        res <- do.call(getMantel, args)
        return(res)
    }
)

#' @rdname getMantel
#' @export
setMethod("getMantel",
    signature = c(x = "ANY"),
    function(x, dist.methods = c("bray", "euclidean"), group.by = NULL, ...) {
        .check_input(x, "list", length = 2L)
        .check_input(dist.methods, "character vector", length = length(x))
        .check_input(
            group.by, c("NULL", "character scalar"),
            length = length(x)
        )
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        if (length(setdiff(colnames(x[[1L]]), colnames(x[[2L]]))) > 0L) {
            stop("Sample names (coluns) must match.")
        }
        # Change orientation so that samples are in rows
        x <- lapply(x, t)
        x <- .get_dissimilarity_matrices(x, dist.methods, ...)
        args <- list(x = x, strata = group.by)
        args <- c(
            args, list(...)[!names(list(...)) %in% c("group.by", "strata")]
        )
        res <- do.call(.run_mantel_test, args)
        return(res)
    }
)


################################ HELP FUNCTIONS ################################

# Get grouping variable from colData of MAE
.get_grouping_for_mantel_mae <- function(x, experiments, group.by = strata, strata = NULL, ...) {
    if (!is.null(group.by)) {
        df <- .retrieve_sample_metadata_from_mae(x, experiments, ...)
        .check_input(
            group.by, "character scalar",
            supported_values = colnames(df)
        )
        group.by <- df[[group.by]]
    }
    return(group.by)
}

# Get grouping variable from colData of TSE
.get_grouping_for_mantel_tse <- function(x, experiments, group.by = strata, strata = NULL, ...) {
    if (!is.null(group.by)) {
        df <- .retrieve_sample_metadata_from_tse(x, experiments, ...)
        .check_input(
            group.by, "character scalar",
            supported_values = colnames(df)
        )
        group.by <- df[[group.by]]
    }
    return(group.by)
}

# Calculate dissimilarity matrices for each assay
#' @importFrom vegan vegdist
.get_dissimilarity_matrices <- function(x, dist.methods, na.rm = FALSE, ...) {
    .check_input(na.rm, "logical scalar")
    x <- lapply(x |> length() |> seq_len(), function(i) {
        vegdist(x[[i]], method = dist.methods[[i]], na.rm = na.rm)
    })
    return(x)
}

# Calculate the Mantel test statistics
#' @importFrom vegan mantel
.run_mantel_test <- function(x, method = "kendall", npermutations = 999, strata = NULL, na.rm = FALSE,
                             parallel = getOption("mc.cores"), ...) {
    .check_input(method, "character scalar")
    .check_input(npermutations, "integer scalar", limits = list(lower = 0))
    .check_input(strata, c("NULL", "character vector"))
    .check_input(na.rm, "logical scalar")
    #
    res <- mantel(
        xdis = x[[1L]],
        ydis = x[[2L]],
        method = method,
        permutations = npermutations,
        strata = strata,
        na.rm = na.rm,
        parallel = parallel
    )
    return(res)
}
