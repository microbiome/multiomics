#' @name getCrossAssociation
#'
#' @title
#' Cross-association between two omic layers
#'
#' @description
#' Calculates pairwise correlations between features in two selected omic
#' datasets.
#'
#' @details
#' For multi-assay objects, the function selects two experiments and assays,
#' retains their shared samples, and calculates a feature-by-feature
#' correlation matrix using \code{stats::cor()}.
#'
#' This initial interface supports numeric assays and the correlation methods
#' supported by \code{stats::cor()}. Group comparisons, covariate models,
#' mixed-effects models, p-values, and table-formatted output are outside the
#' scope of this function.
#'
#' @return
#' A numeric matrix. Rows correspond to features in the first dataset and
#' columns correspond to features in the second dataset.
#'
#' @param x \code{MultiAssayExperiment}, \code{SingleCellExperiment}, or a
#' list containing two numeric matrices.
#'
#' @param experiments \code{Character vector} or \code{integer vector}. Names
#' or indices of the two experiments selected from \code{x}.
#'
#' @param assay.types \code{Character vector}. Names of the assays selected
#' from the two experiments.
#'
#' @param method \code{Character scalar}. Correlation method passed to
#' \code{stats::cor()}: \code{"pearson"}, \code{"spearman"}, or
#' \code{"kendall"}. (Default: \code{"pearson"})
#'
#' @param ... additional arguments used by the selected input method.
#'
#' @examples
#' first <- matrix(
#'     c(1, 2, 3, 4, 4, 3, 2, 1),
#'     nrow = 2L,
#'     dimnames = list(c("feature_a", "feature_b"), paste0("sample", 1:4))
#' )
#' second <- matrix(
#'     c(2, 4, 6, 8, 8, 6, 4, 2),
#'     nrow = 2L,
#'     dimnames = list(c("feature_c", "feature_d"), paste0("sample", 1:4))
#' )
#' getCrossAssociation(list(first, second))
#'
#' @seealso
#' \code{\link[stats:cor]{stats::cor()}}
#'
NULL

#' @rdname getCrossAssociation
#' @export
setMethod(
    "getCrossAssociation",
    signature = c(x = "MultiAssayExperiment"),
    function(x, experiments, assay.types, method = "pearson") {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        mat_list <- .get_shared_samples_from_mae(
            x, experiments, assay.types
        )
        getCrossAssociation(mat_list, method = method)
    }
)

#' @rdname getCrossAssociation
#' @export
setMethod(
    "getCrossAssociation",
    signature = c(x = "SingleCellExperiment"),
    function(x, experiments, assay.types, method = "pearson") {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        mat_list <- .get_shared_samples_from_tse(
            x, experiments, assay.types
        )
        getCrossAssociation(mat_list, method = method)
    }
)

#' @rdname getCrossAssociation
#' @export
setMethod(
    "getCrossAssociation",
    signature = c(x = "ANY"),
    function(x, method = "pearson") {
        .check_input(x, "list", length = 2L)
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        if (!all(vapply(x, is.numeric, logical(1L)))) {
            stop("'x' must include numeric matrices.", call. = FALSE)
        }
        if (is.null(colnames(x[[1L]])) || is.null(colnames(x[[2L]]))) {
            stop("Both matrices must have sample names.", call. = FALSE)
        }
        if (!identical(colnames(x[[1L]]), colnames(x[[2L]]))) {
            stop("Sample names and order must match.", call. = FALSE)
        }
        .check_input(
            method,
            "character scalar",
            supported_values = c("pearson", "spearman", "kendall")
        )
        stats::cor(t(x[[1L]]), t(x[[2L]]), method = method)
    }
)
