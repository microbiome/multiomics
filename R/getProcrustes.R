#' @name
#' getProcrustes
#'
#' @title
#' Procrustes test between two tables
#'
#' @description
#' Performs a Procrustes test to assess the association between two data
#' matrices (e.g., different omics layers) by comparing their ordinations.
#'
#' @details
#' For multi-assay objects, the function:
#' \enumerate{
#'   \item extracts shared samples across selected experiments,
#'   \item computes dissimilarity matrices for each dataset,
#'   \item computes principal coordinates (PCoA) from each dissimilarity
#'   matrix,
#'   \item performs the Procrustes test using \code{vegan::protest()}.
#' }
#'
#' @references
#' Peres-Neto, P.R. & Jackson, D.A. (2001) How well do multivariate data sets
#' match? The advantages of a Procrustean superimposition approach over the
#' Mantel test. Oecologia 129: 169-178.
#'
#' @return
#' An object of class \code{"protest"} as returned by
#' \code{vegan::protest()}.
#'
#' @inheritParams getMantel
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
#' # Perform Procrustes test
#' getProcrustes(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("counts", "nmr")
#' )
#'
#' @seealso
#' \code{\link[vegan:procrustes]{vegan::protest()}},
#' \code{\link[=getMantel]{getMantel()}}
#'
NULL

#' @rdname getProcrustes
#' @export
setMethod("getProcrustes",
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
        res <- getProcrustes(mat_list, ...)
        return(res)
    }
)

#' @rdname getProcrustes
#' @export
setMethod("getProcrustes",
    signature = c(x = "SingleCellExperiment"),
    function(x, experiments, assay.types, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        #
        mat_list <- .get_shared_samples_from_tse(
            x, experiments, assay.types
        )
        res <- getProcrustes(mat_list, ...)
        return(res)
    }
)

#' @rdname getProcrustes
#' @export
setMethod("getProcrustes",
    signature = c(x = "ANY"),
    function(x, dist.methods = c("bray", "euclidean"), ...) {
        .check_input(x, "list", length = 2L)
        .check_input(dist.methods, "character vector", length = length(x))
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        if (length(setdiff(colnames(x[[1L]]), colnames(x[[2L]]))) > 0L) {
            stop("Sample names (columns) must match.", call. = FALSE)
        }
        # Change orientation so that samples are in rows
        x <- lapply(x, t)
        x <- .get_dissimilarity_matrices(x, dist.methods, ...)
        res <- .run_procrustes_test(x, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Calculate PCoA from each dissimilarity matrix and run the Procrustes test
#' @importFrom stats cmdscale
#' @importFrom vegan protest
.run_procrustes_test <- function(x, npermutations = 999, ...) {
    .check_input(npermutations, "integer scalar", limits = list(lower = 0))
    #
    ord <- lapply(x, cmdscale)
    res <- protest(ord[[1L]], ord[[2L]], permutations = npermutations)
    return(res)
}
