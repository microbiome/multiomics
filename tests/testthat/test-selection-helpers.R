make_selection_inputs <- function() {
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
    first_alt <- rbind(
        first_a = c(1, 2, 4, 8, 3, 5, 7, 11),
        first_b = c(2, 1, 5, 9, 4, 6, 10, 12)
    )
    second_alt <- rbind(
        second_a = c(2, 1, 5, 9, 4, 6, 10, 12),
        second_b = c(1, 3, 6, 10, 5, 7, 9, 13)
    )
    colnames(predictors) <- samples
    colnames(responses) <- samples
    colnames(first_alt) <- samples
    colnames(second_alt) <- samples

    make_sce <- function(main, alternative) {
        x <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = main, scaled = main + 100)
        )
        SingleCellExperiment::altExp(x, "alternative") <-
            SingleCellExperiment::SingleCellExperiment(
                assays = list(
                    counts = alternative,
                    scaled = alternative + 100
                )
            )
        x
    }

    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            first = make_sce(predictors, first_alt),
            second = make_sce(responses, second_alt)
        )
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = predictors, scaled = predictors + 100)
    )
    SingleCellExperiment::altExp(sce, "outcomes") <-
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = responses, scaled = responses + 100)
        )

    list(
        predictors = predictors,
        responses = responses,
        first_alt = first_alt,
        second_alt = second_alt,
        mae = mae,
        sce = sce
    )
}

expect_selection_rda_equivalent <- function(actual, expected) {
    expect_s3_class(actual, "rda")
    expect_equal(actual$tot.chi, expected$tot.chi)
    expect_equal(actual$CCA$eig, expected$CCA$eig)
    expect_equal(actual$CA$eig, expected$CA$eig)
    expect_equal(abs(actual$CCA$u), abs(expected$CCA$u))
    expect_equal(abs(actual$CA$u), abs(expected$CA$u))
}

test_that("MAE experiment and assay selection preserve requested order", {
    inputs <- make_selection_inputs()

    result <- multiomics::getRDA(
        inputs$mae,
        experiments = c("second", "first"),
        assay.types = c("scaled", "counts")
    )
    expected <- multiomics::getRDA(
        list(inputs$responses + 100, inputs$predictors)
    )

    expect_selection_rda_equivalent(result, expected)
})

test_that("SCE experiment and assay selection preserve requested order", {
    inputs <- make_selection_inputs()

    named_result <- multiomics::getRDA(
        inputs$sce,
        experiments = c("main", "outcomes"),
        assay.types = c("counts", "counts")
    )
    named_expected <- multiomics::getRDA(
        list(inputs$predictors, inputs$responses)
    )
    indexed_result <- multiomics::getRDA(
        inputs$sce,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts")
    )
    indexed_expected <- multiomics::getRDA(
        list(inputs$responses, inputs$predictors)
    )

    expect_selection_rda_equivalent(named_result, named_expected)
    expect_selection_rda_equivalent(indexed_result, indexed_expected)
})

test_that("numeric and character altExp selection preserve results", {
    inputs <- make_selection_inputs()

    expected <- multiomics::getRDA(
        list(inputs$first_alt, inputs$second_alt)
    )
    numeric_result <- multiomics::getRDA(
        inputs$mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts"),
        altexps = c(1L, 1L)
    )
    character_result <- multiomics::getRDA(
        inputs$mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts"),
        altexps = c("alternative", "alternative")
    )

    expect_selection_rda_equivalent(numeric_result, expected)
    expect_selection_rda_equivalent(character_result, expected)
})

test_that("invalid experiment and assay selections are rejected", {
    inputs <- make_selection_inputs()

    expect_error(
        multiomics::getRDA(
            inputs$mae,
            experiments = c("missing", "first"),
            assay.types = c("counts", "counts")
        ),
        "must specify names or index"
    )
    expect_error(
        multiomics::getRDA(
            inputs$sce,
            experiments = c(0L, 1L),
            assay.types = c("counts", "counts")
        ),
        "must specify names or index"
    )
    expect_error(
        multiomics::getRDA(
            inputs$mae,
            experiments = c("first", "second"),
            assay.types = c("missing", "counts")
        ),
        "Cannot find assay"
    )
})
