#' @name
#' getLASSO
#'
#' @title
#' Sparse association modeling between microbiome and metabolites using LASSO
#'
#' @description
#' Performs LASSO regression to identify associations between microbiome
#' features (predictors) and metabolites (responses). For each metabolite, a
#' separate LASSO model is fitted to select the most relevant taxa and estimate
#' their regression coefficients.
#'
#' @details
#' The function fits a separate LASSO model for each metabolite using the
#' microbiome data as predictors. The procedure is implemented using
#' \code{glmnet} with cross-validation to select the optimal regularization
#' parameter.
#'
#' Non-selected taxa (i.e., coefficients shrunk to zero) are represented as
#' \code{NA}.
#'
#' The resulting matrix includes only core features, i.e., features including
#' only \code{NA}s are filtered out.
#'
#' When applied to log-ratio transformed data (i.e.,
#' \code{apply.logratio=TRUE}), this approach is equivalent to
#' implementation in \code{coda4microbiome::coda_glmnet()}.
#'
#' @references
#' Calle M.L., Pujolassos, M. and Susin A. (2023).
#' coda4microbiome: compositional data analysis for microbiome cross-sectional
#' and longitudinal studies. \emph{BMC Bioinformatics}, 24, 82.
#' https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05205-3
#'
#' @return
#' A numeric matrix of class \code{matrix} with:
#' \itemize{
#'   \item Rows: taxa (predictors)
#'   \item Columns: metabolites (responses)
#'   \item Values: LASSO regression coefficients
#' }
#'
#' @inheritParams getMantel
#'
#' @param ... additional arguments passed to \code{glmnet::cv.glmnet()}.
#' \itemize{
#'   \item \code{apply.logratio}: \code{Logical scalar}. Whether to apply
#'   log-ratio transformation. (Default: \code{TRUE}).
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
#' # Perform LASSO
#' res <- getLASSO(
#'     mae,
#'     experiments = c(1, 2),
#'     assay.types = c("counts", "nmr")
#' )
#'
#' # Visualize results
#' plotLASSO(res)
#'
#' @seealso
#' \code{\link[glmnet:cv.glmnet]{glmnet::cv.glmnet()}}
#'
NULL

#' @rdname getLASSO
#' @export
setMethod("getLASSO",
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
        res <- getLASSO(mat_list, ...)
        return(res)
    }
)

#' @rdname getLASSO
#' @export
setMethod("getLASSO",
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
        res <- getLASSO(mat_list, ...)
        return(res)
    }
)

#' @rdname getLASSO
#' @export
setMethod("getLASSO",
    signature = c(x = "ANY"),
    function(x, ...) {
        if (!(is.list(x) && length(x) == 2L)) {
            stop("'x' must be a list of 2.")
        }
        if (!all(vapply(x, is.matrix, logical(1L)))) {
            stop("'x' must include matrices.", call. = FALSE)
        }
        if (length(setdiff(colnames(x[[1L]]), colnames(x[[2L]]))) > 0L) {
            stop("Sample names (coluns) must match.")
        }
        # Change orientation so that samples are in rows
        x <- lapply(x, t)
        res <- .run_lasso(x, ...)
        return(res)
    }
)

################################ HELP FUNCTIONS ################################

# Runs LASSO for each metabolite and returns taxa × metabolite coefficient
# matrix
.run_lasso <- function(x, apply.logratio = TRUE, ...) {
    .check_input(apply.logratio, "logical scalar")

    # Extract microbiome (predictors) and metabolites (responses)
    microbiome <- x[[1L]]
    metabolites <- x[[2L]] |> as.data.frame()

    # Apply log-ratio transformation if specified
    if (apply.logratio) {
        microbiome <- .apply_logratio(microbiome)
    }

    # Run LASSO separately for each metabolite
    res <- lapply(metabolites, function(single_metabolite) {
        temp <- .run_glmnet(x = microbiome, y = single_metabolite, ...)
        return(temp)
    })

    # Build taxa × metabolite coefficient matrix
    # NA = taxon not selected by LASSO
    coef_res <- lapply(res, function(temp) {
        coef <- rep(NA, ncol(microbiome))
        coef[temp[["selected_taxa"]]] <- temp[["coefficients"]]
        return(coef)
    })
    coef_res <- do.call(cbind, coef_res)
    colnames(coef_res) <- colnames(metabolites)
    rownames(coef_res) <- colnames(microbiome)

    # Remove those taxa and metabolites that do not have any coefficient, i.e.,
    # do not belong to core taxa or metabolites
    core_taxa <- rowSums(is.na(coef_res)) != ncol(coef_res)
    core_metabolites <- colSums(is.na(coef_res)) != nrow(coef_res)
    coef_res <- coef_res[core_taxa, core_metabolites, drop = FALSE]

    # Map results back to feature-level if log-ratio was applied
    if (apply.logratio) {
        coef_res <- .map_pairs_to_taxa(
            coef_res,
            attr(microbiome, "original_names")
        )
    }

    coef_res <- as.data.frame(as.table(coef_res))
    colnames(coef_res) <- c("feature1", "feature2", "value")

    class(coef_res) <- c("LASSO", class(coef_res))

    return(coef_res)
}

# This function applies log-ratio transformation and returns the transformed
# table along with mappings for feature-pairs.
.apply_logratio <- function(x) {
    .require_package("mia")
    # Convert names to indices
    nams <- x |> colnames()
    colnames(x) <- x |>
        ncol() |>
        seq_len()
    x_transf <- mia:::.apply_transformation(x, "division", MARGIN = 2L)
    if (any(x == 0)) {
        x_transf <- mia:::.apply_pseudocount(x_transf, TRUE)
    }
    x_transf <- mia:::.apply_transformation(x_transf, "log10", MARGIN = 2L)
    x_transf <- x_transf |> as.matrix()
    attr(x_transf, "original_names") <- nams
    return(x_transf)
}

# Run LASSO to get which taxa are relevant for single metabolite
# Inspired by coda4microbiome implementation
#' @importFrom glmnet cv.glmnet
.run_glmnet <- function(x, y, lambda = "lambda.1se", alpha = 0.9, nfolds = 10, ...) {
    .check_input(lambda, "character scalar", c("lambda.1se", "lambda.min"))
    .check_input(
        alpha, "numeric scalar",
        limits = list(include_lower = 0, include_upper = 1)
    )
    .check_input(nfolds, "integer scalar", limits = list(lower = 0))
    if (anyNA(x) || anyNA(y)) {
        stop("NAs are not allowed.", call. = FALSE)
    }

    # Fit LASSO with cross-validation
    # Finds optimal regularization parameter (lambda)
    fit <- cv.glmnet(
        x,
        y,
        alpha = alpha,
        type.measure = "deviance",
        nfolds = nfolds
    )

    #  Get lambda that was chose based on CV
    lambda_value <- fit[[lambda]]

    # Get coefficients for taxa
    coefs <- coef(fit, s = lambda_value) |> as.vector()
    # Remove intercept
    coefs <- coefs[-1]
    # Identify taxa selected by LASSO (non-zero coefficients)
    selected <- which(coefs != 0)

    # Predictions on training data (for model evaluation)
    preds <- predict(fit, x, s = lambda_value) |> as.numeric()
    # Calculate explained variance
    var_explained <- NA
    if (sd(preds) > 0 && sd(y) > 0) {
        var_explained <- cor(preds, y)^2
    }

    # Return a list of results
    res <- list(
        selected_taxa = selected,
        coefficients = coefs[selected],
        predictions = preds,
        var_explained = var_explained,
        lambda_1se = fit[["lambda.1se"]],
        lambda_min = fit[["lambda.min"]],
        fit = fit
    )
    return(res)
}

# Map log-ratio pairs back to taxa
.map_pairs_to_taxa <- function(x, original_names) {
    # Split pair names into two indices
    idx <- t(vapply(
        strsplit(rownames(x), "-"),
        function(z) as.integer(z),
        numeric(2L)
    ))
    i <- idx[, 1]
    j <- idx[, 2]

    # Get all indices
    row_indices <- strsplit(rownames(x), "-")
    # Get coefficients for numerator
    pos_mat <- x
    rownames(pos_mat) <- vapply(row_indices, function(x) x[[1L]], character(1L))
    # Get coefficients for denominator
    neg_mat <- -x
    rownames(neg_mat) <- vapply(row_indices, function(x) x[[2L]], character(1L))
    # Create a full matrix
    full_mat <- rbind(pos_mat, neg_mat)
    rownames(full_mat) <- original_names[rownames(full_mat) |> as.integer()]

    # Replace NA with 0 for aggregation
    full_mat[is.na(full_mat)] <- 0
    # Aggregate by rownames (taxa)
    full_mat <- rowsum(full_mat, group = rownames(full_mat), reorder = FALSE)

    return(full_mat)
}

################################### PLOTTING ###################################

#' @rdname getLASSO
#' @export
plot.LASSO <- function(x, sort = TRUE, ...) {
    .check_input(sort, "logical scalar")
    if( sort ){
        x <- .order_features(x)
    }
    p <- .plot_heatmap(x, ...)
    return(p)
}
