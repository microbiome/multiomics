make_rda_inputs <- function() {
    samples <- paste0("sample", seq_len(8L))
    predictors <- rbind(
        feature_a = seq_len(8L),
        feature_b = c(2, 1, 4, 3, 6, 5, 8, 7),
        feature_c = c(8, 3, 6, 1, 7, 2, 5, 4)
    )
    responses <- rbind(
        response_a = c(3, 1, 4, 2, 5, 7, 6, 8),
        response_b = c(8, 6, 4, 2, 1, 3, 5, 7)
    )
    colnames(predictors) <- samples
    colnames(responses) <- samples
    list(predictors = predictors, responses = responses)
}

make_rda_mae <- function(inputs) {
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            predictors = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$predictors)
            ),
            responses = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$responses)
            )
        )
    )
}

test_that("getRDA matches vegan::rda for list input", {
    inputs <- make_rda_inputs()

    result <- getRDA(list(inputs$predictors, inputs$responses))
    expected <- vegan::rda(
        X = t(inputs$predictors),
        Y = t(inputs$responses)
    )

    expect_s3_class(result, "rda")
    expect_equal(result$tot.chi, expected$tot.chi)
    expect_equal(result$CCA$eig, expected$CCA$eig)
    expect_equal(result$CA$eig, expected$CA$eig)
    expect_equal(abs(result$CCA$u), abs(expected$CCA$u))
    expect_equal(abs(result$CA$u), abs(expected$CA$u))
})

test_that("getRDA preserves results through a MultiAssayExperiment", {
    inputs <- make_rda_inputs()
    mae <- make_rda_mae(inputs)

    direct_result <- getRDA(list(inputs$predictors, inputs$responses))
    mae_result <- getRDA(
        mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts")
    )

    expect_s3_class(mae_result, "rda")
    expect_equal(mae_result$tot.chi, direct_result$tot.chi)
    expect_equal(mae_result$CCA$eig, direct_result$CCA$eig)
    expect_equal(mae_result$CA$eig, direct_result$CA$eig)
    expect_equal(abs(mae_result$CCA$u), abs(direct_result$CCA$u))
    expect_equal(abs(mae_result$CA$u), abs(direct_result$CA$u))
})

test_that("getRDA validates list inputs and arguments", {
    inputs <- make_rda_inputs()

    expect_error(
        getRDA(inputs$predictors),
        "must be list"
    )
    expect_error(
        getRDA(list(inputs$predictors)),
        "length must be 2"
    )
    expect_error(
        getRDA(list(as.data.frame(inputs$predictors), inputs$responses)),
        "include matrices"
    )

    mismatched <- inputs$responses
    colnames(mismatched)[[1L]] <- "not-a-sample"
    expect_error(
        getRDA(list(inputs$predictors, mismatched)),
        "Sample names"
    )
    expect_error(
        getRDA(
            list(inputs$predictors, inputs$responses),
            scale = NA
        )
    )
})
