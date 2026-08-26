#' @name
#' getPairwiseAssociation
#'
#' @title
#' Pairwise association analysis between two feature matrices
#'
#' @description
#' Calculates pairwise associations between features in two data matrices.
#' Each feature in the first matrix is associated with each feature in the
#' second matrix. By default, associations are calculated using Kendall
#' correlation.
#'
#' @details
#' The function calculates the association between every pair of features
#' from the two input matrices.
#'
#' If \code{group.by} is specified, associations are calculated separately
#' within each group defined by the corresponding column in \code{col.data}.
#' This can be used to obtain, for example, sex-specific or treatment-specific
#' associations.
#'
#' If \code{formula} is specified, the effects of the covariates in the
#' formula are accounted for before calculating associations. For fixed-effect
#' formulas, each feature is regressed on the specified covariates and the
#' association is calculated between the resulting residuals. Thus, the
#' resulting association represents the relationship between the features
#' after accounting for the specified covariates.
#'
#' @return
#' A \code{data.frame} containing pairwise association estimates in long
#' format. \code{feature1} and \code{feature2} identify the feature pair,
#' and \code{value} contains the association estimate. If \code{group.by}
#' is specified, \code{strata} identifies the group for which the
#' association was calculated.
#'
#' @inheritParams getMantel
#'
#' @param col.data \code{data.frame}. Contains sample metadata used when
#' \code{group.by} is specified.
#'
#' @param method \code{Character scalar}. Specifies the association method.
#' One of \code{"kendall"}, \code{"pearson"}, or \code{"spearman"}.
#' (Default: \code{"kendall"}).
#'
#' @param group.by \code{Character scalar}. Column name in sample metadata
#' used to calculate associations separately within each group. If \code{NULL},
#' associations are calculated across all samples. (Default: \code{NULL})
#'
#' @param ... additional arguments:
#'  \itemize{
#'   \item \code{formula}: \code{Formula}. Specifies covariates to account for.
#'   Covariate effects are removed by fitting models and calculating
#'   associations between the resulting residuals.
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
#' # Agglomerate microbiome data
#' mae[[1]] <- agglomerateByRank(mae[[1]], rank = "Phylum")
#'
#' # Perform correlation analysis
#' res <- getPairwiseAssociation(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("counts", "nmr")
#' )
#'
#' # Visualize results
#' plotPairwiseAssociation(res)
#'
#' @seealso
#' \code{\link[getLASSO]{multiomics::getLASSO()}}
#'
NULL

#' @rdname getPairwiseAssociation
#' @export
setMethod("getPairwiseAssociation",
    signature = c(x = "MultiAssayExperiment"),
    function(x, experiments, assay.types, group.by = NULL, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        mat_list <- .get_shared_samples_from_mae(x, experiments, assay.types, ...)
        sample_data <- .retrieve_sample_metadata_from_mae(x, experiments)

        res <- getPairwiseAssociation(x = mat_list, col.data = sample_data, group.by = group.by, ...)
        return(res)
    }
)

#' @rdname getPairwiseAssociation
#' @export
setMethod("getPairwiseAssociation",
    signature = c(x = "SingleCellExperiment"),
    function(x, experiments, assay.types, group.by = NULL, ...) {
        .check_input(
            experiments, c("character vector", "integer vector"),
            length = 2L
        )
        .check_input(assay.types, "character vector", length = 2L)
        mat_list <- .get_shared_samples_from_tse(x, experiments, assay.types, ...)
        sample_data <- .retrieve_sample_metadata_from_tse(x, experiments)
        res <- getPairwiseAssociation(x = mat_list, col.data = sample_data, group.by = group.by, ...)
        return(res)
    }
)

#' @rdname getPairwiseAssociation
#' @export
setMethod("getPairwiseAssociation",
    signature = c(x = "ANY"),
    function(x, col.data = NULL, group.by = NULL, method = "kendall", ...) {
        .check_input(x, "list", length = 2L)
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        # Check sample names match between matrices
        samples_1 <- colnames(x[[1L]])
        samples_2 <- colnames(x[[2L]])
        if (!all(samples_1 %in% samples_2) || !all(samples_2 %in% samples_1)) {
            stop("Sample names (columns) must match between matrices.", call. = FALSE)
        }
        .check_input(col.data, c("NULL", "data.frame"))

        if( !is.null(group.by) ){
            if (is.null(col.data)) {
                stop("'col.data' (sample metadata) required when 'group.by' is specified.", call. = FALSE)
            }
            .check_input(group.by, "character scalar", supported_values = colnames(col.data))
        }
        .check_input(
            method,
            "character vector",
            supported_values = c("kendall", "pearson", "spearman")
        )

        # Change orientation so that samples are in rows
        x <- lapply(x, t)
        res <- .calculate_association(
            x = x[[1L]], y = x[[2L]],
            sample_data = col.data, group.by = group.by, method = method, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Main helper function for calculating correlation in different scenarios.
.calculate_association <- function(
        x, y,
        sample_data = NULL, group.by = NULL,
        formula = NULL,
        method = "kendall", ...
) {
    # If formula specified, validate metadata
    if (!is.null(formula)) {
        .check_input(formula, "formula")
        if (is.null(sample_data)) {
            stop("'col.data' (sample metadata) required when 'formula' is specified.", call. = FALSE)
        }
        if (!all(rownames(x) %in% rownames(sample_data))) {
            stop("All samples must be present in 'col.data' rownames.", call. = FALSE)
        }

        # Check all covariates in formula exist in metadata
        covariates <- .get_rhs(formula)
        missing_vars <- setdiff(covariates, colnames(sample_data))
        if (length(missing_vars) > 0L) {
            stop("Covariates not found in 'col.data': '",
                 paste0(missing_vars, collapse = "', '"), "'",
                 call. = FALSE)
        }
    }
    # Select association function
    if (is.null(formula)) {
        FUN <- .association_cor
    } else {
        rhs <- attr(terms(formula), "term.labels")
        has_random <- any(grepl("\\|", rhs))

        FUN <- if (has_random) .association_lmm else .association_lm
    }

    # Calculate associations
    if (is.null(group.by)) {

        res <- FUN(
            x = x,
            y = y,
            sample_data = sample_data,
            formula = formula,
            method = method,
            ...
        )

    } else {
        # Calculate correlation for each group separately
        strata <- sample_data[[group.by]]

        groups <- unique(strata[!is.na(strata)])

        res <- lapply(groups, function(g) {
            keep <- !is.na(strata) & strata == g

            FUN(
                x = x[keep, , drop = FALSE],
                y = y[keep, , drop = FALSE],
                sample_data = sample_data[keep, , drop = FALSE],
                formula = formula,
                method = method,
                ...
            )
        })

        names(res) <- groups
    }

    # Convert association matrix/matrices to long format
    if (is.list(res)) {
        res <- do.call(
            rbind,
            lapply(names(res), function(g) {
                out <- as.data.frame(as.table(res[[g]]))
                colnames(out) <- c("feature1", "feature2", "value")
                out[["strata"]] <- g
                return(out)
            })
        )
    } else {
        res <- as.data.frame(as.table(res))
        colnames(res) <- c("feature1", "feature2", "value")
    }

    rownames(res) <- NULL

    class(res) <- c("PairwiseAssociation", class(res))

    return(res)
}

# Calculate simple associaton between two tables
.association_cor <- function(x, y, method = "kendall", ...) {

    res <- cor(
        x,
        y,
        use = "pairwise.complete.obs",
        method = method
    )

    return(res)
}

# Correlation on residuals. First we remove the effect of covariates and then
# calculate the correlation.
.association_lm <- function(x, y, sample_data, formula, method = "kendall", ...) {

    covariates <- attr(terms(formula), "term.labels")
    model_formula <- reformulate(covariates)

    Z <- model.matrix(
        model_formula,
        data = sample_data
    )

    x_adj <- lm.fit(Z, x) |> residuals()
    y_adj <- lm.fit(Z, y) |> residuals()

    res <- .association_cor(x_adj, y_adj, method = method)

    return(res)
}

# Calculate association after adjustment using a linear mixed model.
.association_lmm <- function(x, y, sample_data, formula, method = "kendall", ...) {
    # Check lme4 availability
    .require_package("lme4")

    # Modify formula so that it is suitable for lme4
    covariates <- attr(terms(formula), "term.labels")
    formula <- covariates |> reformulate()
    formula <- update(formula, feature ~ .)

    # Extract residuals for each feature in x and y
    x_adj <- .extract_lmm_residuals(x, sample_data, formula, ...)
    y_adj <- .extract_lmm_residuals(y, sample_data, formula, ...)

    # Correlation on residuals
    res <- .association_cor(x_adj, y_adj, method = method)

    return(res)
}

# Fit mixed model for each feature and extract residuals
.extract_lmm_residuals <- function(feature_matrix, sample_data, formula, ...) {
    feature_matrix <- feature_matrix |> as.data.frame()

    residuals_list <- lapply(feature_matrix, function(col) {
        data_model <- cbind.data.frame(sample_data, feature = col)
        model <- lme4::lmer(formula, data = data_model, REML = TRUE)
        res <- model |> residuals()
        return(res)
    })

    res <- do.call(cbind, residuals_list)

    return(res)
}

################################### PLOTTING ###################################

#' @rdname getPairwiseAssociation
#' @export
plot.PairwiseAssociation <- function(x, sort = TRUE, ...) {
    .check_input(sort, "logical scalar")
    if( sort ){
        x <- .order_features(x)
    }
    p <- .plot_heatmap(x, ...)
    return(p)
}

.order_features <- function(df) {

    # Remove associations with missing/non-finite values
    df <- df[is.finite(df$value), , drop = FALSE]

    mat <- xtabs(value ~ feature1 + feature2, data = df)

    feature1_order <- rownames(mat)[hclust(dist(mat))$order]
    feature2_order <- colnames(mat)[hclust(dist(t(mat)))$order]

    df$feature1 <- factor(df$feature1, levels = feature1_order)
    df$feature2 <- factor(df$feature2, levels = feature2_order)

    return(df)
}

# Create the heatmap
#' @importFrom ggplot2 aes element_rect element_text geom_tile ggplot labs
#'     scale_fill_gradient2 theme theme_minimal facet_wrap
.plot_heatmap <- function(df, ...) {

    # Create the heatmap
    p <- ggplot(
        df,
        aes(
            x = feature1,
            y = feature2,
            fill = value
        )
    ) +
        geom_tile(
            width = 0.95,
            height = 0.95,
            colour = "white",
            linewidth = 0.2
        )

    # Apply the heatmap theme
    p <- p +
        theme_minimal() +
        theme(
            strip.background = element_rect(
                fill = "white",
                colour = "black",
                linewidth = 0.5
            ),
            strip.text = element_text(
                face = "bold",
                colour = "black"
            )
        )

    # Remove axis titles
    p <- p +
        labs(
            x = NULL,
            y = NULL
        )

    # Apply the fill scale
    p <- p +
        scale_fill_gradient2(
            low = "#2166AC",
            mid = "white",
            high = "#B2182B",
            midpoint = 0
        )

    # Facet by strata when group-wise associations were calculated
    if ("strata" %in% colnames(df)) {
        p <- p +
            facet_wrap(~ strata)
    }
    return(p)
}
