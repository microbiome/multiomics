#' @name plotLASSO
#'
#' @title
#' Plot LASSO associations
#'
#' @description
#' Creates a heatmap of the coefficient matrix returned by
#' \code{\link{getLASSO}}.
#'
#' @param x A numeric coefficient matrix, such as the result returned by
#' \code{\link{getLASSO}}.
#'
#' @return
#' A \code{ggplot} object displaying LASSO coefficients as a heatmap.
#'
#' @details
#' The plot uses blue, white, and red for negative, zero, and positive
#' coefficients, respectively. The \pkg{ggplot2} package is required when
#' this function is called.
#'
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     coefficients <- matrix(
#'         c(1, NA, -2, 3, 4, NA),
#'         nrow = 2L,
#'         dimnames = list(
#'             c("taxon_a", "taxon_b"),
#'             c("metabolite_a", "metabolite_b", "metabolite_c")
#'         )
#'     )
#'     plotLASSO(coefficients)
#' }
#'
#' @seealso
#' \code{\link{getLASSO}}
#'
#' @export
plotLASSO <- function(x) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("'x' must be a numeric matrix.", call. = FALSE)
    }
    .require_package("ggplot2")

    taxa <- rownames(x)
    if (is.null(taxa)) {
        taxa <- seq_len(nrow(x))
    }
    metabolites <- colnames(x)
    if (is.null(metabolites)) {
        metabolites <- seq_len(ncol(x))
    }
    plot_data <- data.frame(
        Taxa = rep(taxa, times = ncol(x)),
        Metabolite = rep(metabolites, each = nrow(x)),
        Coefficient = as.vector(x),
        stringsAsFactors = FALSE
    )

    ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
            x = Metabolite,
            y = Taxa,
            fill = Coefficient
        )
    ) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = 0
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
}

utils::globalVariables(c("Coefficient", "Metabolite", "Taxa"))
