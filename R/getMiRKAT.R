#' @name
#' getMiRKAT
#'
#' @title
#' Microbiome Regression-based Kernel Association Test (MiRKAT)
#'
#' @description
#' Performs association testing between microbiome composition and a second
#' omic layer (e.g., metabolites) using the Microbiome Regression-based Kernel
#' Association Test (MiRKAT). The method evaluates whether variation in
#' microbiome profiles is associated with the outcome by comparing kernel
#' matrices derived from dissimilarity measures.
#'
#' @details
#' MiRKAT tests the association between microbiome composition and another omic
#' layer:
#' \enumerate{
#'   \item computing distance matrices from microbiome data using specified
#'     methods,
#'   \item transforming distances into kernel matrices,
#'   \item fitting kernel-based regression models,
#'   \item combining results across multiple kernels using an omnibus test.
#' }
#'
#' For multi-assay objects, the function:
#' \enumerate{
#'   \item extracts shared samples across selected experiments,
#'   \item uses the first dataset as microbiome input,
#'   \item uses the second dataset as outcome(s) (e.g., metabolites),
#'   \item performs MiRKAT separately for each outcome variable.
#' }
#'
#' Multiple distance metrics can be used simultaneously, and results are
#' combined using the Cauchy combination test (default).
#'
#' @references
#' Zhao et al. (2015). Testing in microbiome-profiling studies with MiRKAT,
#' the microbiome regression-based kernel association test.
#' \emph{American Journal of Human Genetics}, 96(5), 797–807.
#'
#' @return
#' A \code{data.frame} where each row corresponds to an outcome variable
#' (e.g., metabolite) and columns contain:
#' \itemize{
#'   \item p-values for each distance metric,
#'   \item an omnibus p-value combining all kernels.
#' }
#'
#' @inheritParams getMantel
#'
#' @param dist.methods \code{Character vector}. Distance metrics used to compute
#' microbiome dissimilarities. Must be values from \code{vegan::vegdist()}.
#' (Default: \code{c("euclidean", "manhattan", "canberra")}).
#'
#' @param ... additional arguments to \code{MiRKAT::MiRKAT()}
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
#' # Perform MiRKAT
#' getMiRKAT(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("counts", "nmr")
#' )
#'
#' @seealso
#' \code{\link[MiRKAT:MiRKAT]{MiRKAT::MiRKAT()}}
#'
NULL

#' @rdname getMiRKAT
#' @export
setMethod("getMiRKAT",
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
        res <- getMiRKAT(mat_list, ...)
        return(res)
    }
)

#' @rdname getMiRKAT
#' @export
setMethod("getMiRKAT",
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
        res <- getMiRKAT(mat_list, ...)
        return(res)
    }
)

#' @rdname getMiRKAT
#' @export
setMethod("getMiRKAT",
    signature = c(x = "ANY"),
    function(x, dist.methods = c("euclidean", "manhattan", "canberra"), ...) {
        .check_input(x, "list", length = 2L)
        .check_input(dist.methods, "character vector")
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        if (length(setdiff(colnames(x[[1L]]), colnames(x[[2L]]))) > 0L) {
            stop("Sample names (coluns) must match.")
        }
        # Change orientation so that samples are in rows
        x <- lapply(x, t)
        x <- .get_microbiome_kernels(x, dist.methods)
        res <- .run_mirkat(x, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# This function computes microbiome kernel matrices from abundance data.
# It calculates distance matrices using specified methods and transforms
# them into kernel matrices required by MiRKAT.
#' @importFrom vegan vegdist
#' @importFrom MiRKAT D2K
.get_microbiome_kernels <- function(x, dist.methods, ...) {
    microbiome_dists <- lapply(dist.methods, function(method) {
        dist_mat <- vegdist(x[[1L]], method = method)
        dist_mat <- dist_mat |> as.matrix()
        dist_mat <- dist_mat |> D2K()
        return(dist_mat)
    })
    names(microbiome_dists) <- dist.methods
    x[[1L]] <- microbiome_dists
    return(x)
}

# This function runs MiRKAT for each metabolite using precomputed
# microbiome kernel matrices. It performs association testing between
# microbiome composition and each metabolite and returns p-values
# (including kernel-specific and omnibus tests).
#' @importFrom MiRKAT MiRKAT
#' @importFrom stats p.adjust
.run_mirkat <- function(x, method = "davies", npermutations = nperm, nperm = 999,
                        omnibus = "cauchy", p.adjust.method = "fdr", ...) {
    .check_input(method, "character scalar")
    .check_input(npermutations, "integer scalar", limits = list(lower = 0))
    .check_input(omnibus, "character scalar")
    .check_input(
        p.adjust.method, "character scalar",
        supported_values = stats::p.adjust.methods
    )
    #
    microbiome_dists <- x[[1L]]
    metabolites <- x[[2L]] |> as.data.frame()

    # Calculate association of each metabolite to microbiome composition
    res <- lapply(metabolites, function(single_metabolite) {
        temp <- MiRKAT(
            y = single_metabolite,
            Ks = microbiome_dists,
            method = method,
            omnibus = omnibus,
            out_type = "C",
            nperm = npermutations,
            returnKRV = TRUE,
            returnR2 = TRUE,
        )
        temp <- temp |> unlist()
        return(temp)
    })
    res <- do.call(rbind, res) |> as.data.frame()
    # Tidy up the names
    colnames(res) <- gsub("\\.", "_", colnames(res)) |> tolower()

    # Adjust p-values because of multiple testing
    if( "omnibus_p" %in% colnames(res) ){
        res[["omnibus_p_adj"]] <- p.adjust(
            res[["omnibus_p"]],
            method = p.adjust.method
        )
    }

    class(res) <- c("MiRKAT", class(res))

    return(res)
}

#' @rdname getMiRKAT
#' @export
plot.MiRKAT <- function(
        x,
        column = NULL,
        ...
) {
    .check_input(
        column,
        c("NULL", "character vector"),
        supported_values = colnames(x)
    )

    if( is.null(column) ){
        column <- colnames(x)
    }

    # Add feature names
    x[["feature"]] <- rownames(x)

    # Convert to long format
    df <- x |>
        tidyr::pivot_longer(
            cols = all_of(column),
            names_to = "metric",
            values_to = "value"
        )

    # Order features separately for each metric
    df <- df |>
        dplyr::group_by(metric) |>
        dplyr::mutate(
            feature = factor(
                feature,
                levels = feature[
                    order(value, decreasing = TRUE)
                ]
            )
        ) |>
        dplyr::ungroup()

    # Plot
    p <- ggplot(
        df,
        aes(
            x = value,
            y = feature
        )
    ) +
        geom_segment(
            aes(
                x = 0,
                xend = value,
                y = feature,
                yend = feature
            )
        ) +
        geom_point(
            shape = 21,
            fill = "black"
        ) +
        facet_wrap(
            ~ metric,
            scales = "free_x"
        ) +
        labs(
            x = NULL,
            y = NULL
        ) +
        theme_minimal() +
        theme(
            strip.background = element_rect(
                fill = "white",
                colour = "black"
            ),
            strip.text = element_text(
                face = "bold"
            )
        )

    return(p)
}
