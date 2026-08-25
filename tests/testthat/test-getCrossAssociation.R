make_cross_association_inputs <- function() {
    samples <- paste0("sample", seq_len(5L))
    first <- rbind(
        feature_a = c(1, 2, 3, 5, 4),
        feature_b = c(4, 1, 5, 2, 3),
        feature_c = c(2, 5, 1, 4, 3)
    )
    second <- rbind(
        feature_d = c(5, 4, 3, 2, 1),
        feature_e = c(1, 3, 5, 2, 4)
    )
    colnames(first) <- samples
    colnames(second) <- samples
    list(first = first, second = second)
}

make_cross_association_mae <- function(inputs) {
    second_order <- c("sample3", "sample1", "sample5", "sample2", "sample4")
    second <- inputs$second[, second_order, drop = FALSE]
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            first = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$first)
            ),
            second = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = second)
            )
        )
    )
}

test_that("getCrossAssociation matches stats::cor for list input", {
    inputs <- make_cross_association_inputs()

    result <- getCrossAssociation(list(inputs$first, inputs$second))
    expected <- stats::cor(t(inputs$first), t(inputs$second))

    expect_true(is.matrix(result))
    expect_equal(result, expected)
})

test_that("getCrossAssociation aligns shared MAE samples", {
    inputs <- make_cross_association_inputs()
    mae <- make_cross_association_mae(inputs)

    result <- getCrossAssociation(
        mae,
        experiments = c("first", "second"),
        assay.types = c("counts", "counts")
    )
    expected <- stats::cor(t(inputs$first), t(inputs$second))

    expect_equal(result, expected)
})

test_that("getCrossAssociation supports SingleCellExperiment selection", {
    inputs <- make_cross_association_inputs()
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = inputs$first)
    )
    SingleCellExperiment::altExp(sce, "second") <-
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = inputs$second)
        )

    result <- getCrossAssociation(
        sce,
        experiments = c("main", "second"),
        assay.types = c("counts", "counts")
    )
    expected <- stats::cor(t(inputs$first), t(inputs$second))

    expect_equal(result, expected)
})

test_that("getCrossAssociation validates inputs at the public boundary", {
    inputs <- make_cross_association_inputs()
    mae <- make_cross_association_mae(inputs)

    expect_error(
        getCrossAssociation(inputs$first),
        "must be a list"
    )
    expect_error(
        getCrossAssociation(list(inputs$first)),
        "length must be 2"
    )
    expect_error(
        getCrossAssociation(list(as.data.frame(inputs$first), inputs$second)),
        "include matrices"
    )
    expect_error(
        getCrossAssociation(
            list(inputs$first, inputs$second[, -1L, drop = FALSE])
        ),
        "Sample names"
    )
    expect_error(
        getCrossAssociation(list(inputs$first, inputs$second), method = "bad"),
        "one of the following options"
    )
    expect_error(
        getCrossAssociation(
            mae,
            experiments = c("first", "missing"),
            assay.types = c("counts", "counts")
        ),
        "experiments"
    )
    expect_error(
        getCrossAssociation(
            mae,
            experiments = c("first", "second"),
            assay.types = c("counts", "missing")
        ),
        "Cannot find assay"
    )
})
