make_lasso_inputs <- function() {
    samples <- paste0("sample", seq_len(12L))
    predictors <- rbind(
        taxon_a = seq_len(12L),
        taxon_b = c(3, 1, 4, 2, 5, 7, 6, 8, 9, 11, 10, 12),
        taxon_c = c(12, 4, 9, 1, 8, 2, 10, 3, 7, 5, 11, 6)
    )
    colnames(predictors) <- samples
    responses <- rbind(
        metabolite_a = 2 * predictors["taxon_a", ] + rnorm(12L, sd = 0.1),
        metabolite_b = -1.5 * predictors["taxon_b", ] + rnorm(12L, sd = 0.1)
    )
    colnames(responses) <- samples
    list(predictors = predictors, responses = responses)
}

make_lasso_mae <- function(inputs) {
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            microbiome = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$predictors)
            ),
            metabolites = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$responses)
            )
        )
    )
}

test_that("getLASSO returns a named coefficient matrix", {
    set.seed(42)
    inputs <- make_lasso_inputs()

    result <- getLASSO(
        list(inputs$predictors, inputs$responses),
        apply.logratio = FALSE,
        alpha = 1,
        nfolds = 3L
    )

    expect_true(is.matrix(result))
    expect_type(result, "double")
    expect_equal(ncol(result), nrow(inputs$responses))
    expect_identical(colnames(result), rownames(inputs$responses))
    expect_true(all(rownames(result) %in% rownames(inputs$predictors)))
    expect_gt(nrow(result), 0L)
    expect_true(any(!is.na(result)))
})

test_that("getLASSO matches cv.glmnet at lambda.1se", {
    set.seed(42)
    inputs <- make_lasso_inputs()

    set.seed(42)
    result <- getLASSO(
        list(inputs$predictors, inputs$responses),
        apply.logratio = FALSE,
        alpha = 1,
        nfolds = 3L
    )

    set.seed(42)
    expected <- lapply(seq_len(nrow(inputs$responses)), function(i) {
        fit <- glmnet::cv.glmnet(
            x = t(inputs$predictors),
            y = inputs$responses[i, ],
            alpha = 1,
            type.measure = "deviance",
            nfolds = 3L
        )
        coefficients <- as.vector(coef(fit, s = fit$lambda.1se))[-1L]
        coefficients[coefficients == 0] <- NA_real_
        coefficients
    })
    expected <- do.call(cbind, expected)
    dimnames(expected) <- list(
        rownames(inputs$predictors),
        rownames(inputs$responses)
    )
    expected <- expected[
        rowSums(is.na(expected)) != ncol(expected),
        colSums(is.na(expected)) != nrow(expected),
        drop = FALSE
    ]

    expect_equal(result, expected)
})

test_that("getLASSO preserves results through a MultiAssayExperiment", {
    set.seed(42)
    inputs <- make_lasso_inputs()
    set.seed(42)
    direct_result <- getLASSO(
        list(inputs$predictors, inputs$responses),
        apply.logratio = FALSE,
        alpha = 1,
        nfolds = 3L
    )

    mae <- make_lasso_mae(inputs)
    set.seed(42)
    mae_result <- getLASSO(
        mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts"),
        apply.logratio = FALSE,
        alpha = 1,
        nfolds = 3L
    )

    expect_identical(dimnames(mae_result), dimnames(direct_result))
    expect_equal(mae_result, direct_result)
})

test_that("getLASSO rejects mismatched sample names", {
    inputs <- make_lasso_inputs()
    colnames(inputs$responses)[[1L]] <- "not-a-sample"

    expect_error(
        getLASSO(
            list(inputs$predictors, inputs$responses),
            apply.logratio = FALSE,
            alpha = 1,
            nfolds = 3L
        ),
        "Sample names"
    )
})

test_that("getLASSO validates list inputs and arguments", {
    set.seed(42)
    inputs <- make_lasso_inputs()

    expect_error(
        getLASSO(inputs$predictors, apply.logratio = FALSE),
        "list of 2"
    )
    expect_error(
        getLASSO(list(inputs$predictors), apply.logratio = FALSE),
        "list of 2"
    )
    expect_error(
        getLASSO(
            list(as.data.frame(inputs$predictors), inputs$responses),
            apply.logratio = FALSE
        ),
        "include matrices"
    )
    expect_error(
        getLASSO(
            list(inputs$predictors, inputs$responses),
            apply.logratio = FALSE,
            alpha = 1.1
        ),
        "alpha"
    )
    expect_error(
        getLASSO(
            list(inputs$predictors, inputs$responses),
            apply.logratio = FALSE,
            nfolds = 0L
        ),
        "nfolds"
    )
    expect_error(
        getLASSO(
            list(inputs$predictors, inputs$responses),
            apply.logratio = FALSE,
            nfolds = 3.5
        ),
        "nfolds"
    )

    mae <- make_lasso_mae(inputs)
    expect_error(
        getLASSO(
            mae,
            experiments = c(1L, 2L),
            assay.types = c("counts"),
            apply.logratio = FALSE
        ),
        "assay.types"
    )
    expect_error(
        getLASSO(
            mae,
            experiments = c(1L, 2L),
            assay.types = c("missing", "counts"),
            apply.logratio = FALSE
        ),
        "Cannot find assay"
    )
})
