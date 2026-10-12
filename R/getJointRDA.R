#' @name
#' getJointRDA
#'
#' @title
#' Redundancy analysis between two omic layers
#'
#' @description
#' Performs redundancy analysis (RDA) to assess linear associations between two
#' datasets (e.g., microbiome composition and metabolite profiles). The method
#' evaluates how much variation in one dataset can be explained by the other.
#'
#' @details
#' The function extracts matched samples from the selected experiments and
#' assays, aligns them, and applies redundancy analysis using
#' \code{vegan::rda()}.
#'
#' RDA is a constrained ordination technique equivalent to multivariate linear
#' regression followed by principal component analysis (PCA) on fitted values.
#'
#' @references
#' Legendre, P., & Legendre, L. (2012). \emph{Numerical Ecology}.
#' Elsevier.
#'
#' @return
#' An object of class \code{rda} as returned by \code{vegan::rda()}, containing
#' ordination results, explained variance, and model statistics.
#'
#' @inheritParams getMantel
#'
#' @param ... additional arguments passed to \code{vegan::rda()}.
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
#' # Perform RDA
#' getJointRDA(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("counts", "nmr")
#' )
#'
#' @seealso
#' \code{\link[vegan:rda]{vegan::rda()}}
#'
NULL

#' @rdname getJointRDA
#' @export
setMethod("getJointRDA",
    signature = c(x = "MultiAssayExperiment"),
    function(x, experiments, assay.types, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        mat_list <- .get_shared_samples_from_mae(
            x, experiments, assay.types, ...
        )
        res <- getJointRDA(mat_list, ...)
        return(res)
    }
)

#' @rdname getJointRDA
#' @export
setMethod("getJointRDA",
    signature = c(x = "SingleCellExperiment"),
    function(x, experiments, assay.types, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        mat_list <- .get_shared_samples_from_tse(
            x, experiments, assay.types, ...
        )
        res <- getJointRDA(mat_list, ...)
        return(res)
    }
)

#' @rdname getJointRDA
#' @export
setMethod("getJointRDA",
    signature = c(x = "ANY"),
    function(x, ...) {
        .check_input(x, "list", length = 2L)
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        if (length(setdiff(colnames(x[[1L]]), colnames(x[[2L]]))) > 0L) {
            stop("Sample names (coluns) must match.", call. = FALSE)
        }
        # Change orientation so that samples are in rows
        x <- lapply(x, t)
        res <- .run_rda(x, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

#' @importFrom vegan rda
.run_rda <- function(x, ...) {
    res <- rda(X = x[[1L]], Y = x[[2L]], ...)
    class(res) <- c("JointRDA", class(res))
    return(res)
}

################################### PLOTTING ###################################

#' @rdname getJointRDA
#' @export
#' @importFrom ggvegan ordiggplot geom_ordi_axis geom_ordi_point geom_ordi_arrow
plot.JointRDA <- function(x, show.species = TRUE, show.biplot = TRUE, ...) {
    .check_input(show.species, "logical scalar")
    .check_input(show.biplot, "logical scalar")

    p <- ordiggplot(x) +
        geom_ordi_axis() +
        geom_ordi_point("sites")
    if( show.species ){
        p <- p + geom_ordi_arrow("species", colour = "red")
    }
    if( show.biplot ){
        p <- p + geom_ordi_arrow("biplot", colour = "blue")
    }

    return(p)
}
