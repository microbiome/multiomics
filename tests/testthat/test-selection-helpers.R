make_selection_mae_for_utils <- function() {
    samples <- paste0("sample", seq_len(3L))
    feature_names <- c("feature_a", "feature_b")

    make_sce <- function(offset) {
        counts <- matrix(
            offset + seq_len(6L),
            nrow = 2L,
            dimnames = list(feature_names, samples)
        )
        scaled <- counts + 100
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = counts, scaled = scaled)
        )
    }

    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            first = make_sce(0),
            second = make_sce(10)
        )
    )
}

make_selection_sce_for_utils <- function() {
    samples <- paste0("sample", seq_len(3L))
    feature_names <- c("feature_a", "feature_b")
    main_counts <- matrix(
        seq_len(6L),
        nrow = 2L,
        dimnames = list(feature_names, samples)
    )
    alt_counts <- main_counts + 20

    x <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = main_counts, scaled = main_counts + 100)
    )
    SingleCellExperiment::altExp(x, "outcomes") <-
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = alt_counts, scaled = alt_counts + 100)
        )
    x
}

test_that("MAE experiment selection preserves requested order", {
    mae <- make_selection_mae_for_utils()

    named <- multiomics:::.select_experiments(mae, c("second", "first"))
    indexed <- multiomics:::.select_experiments(mae, c(2L, 1L))

    expect_identical(names(named), c("second", "first"))
    expect_identical(names(indexed), c("second", "first"))
    expect_identical(
        SummarizedExperiment::assay(named[[1L]], "counts"),
        SummarizedExperiment::assay(mae[["second"]], "counts")
    )
})

test_that("MAE experiment selection rejects invalid selections", {
    mae <- make_selection_mae_for_utils()

    expect_error(
        multiomics:::.select_experiments(mae, c("missing", "first")),
        "must specify names or index"
    )
    expect_error(
        multiomics:::.select_experiments(mae, c(0L, 1L)),
        "must specify names or index"
    )
    expect_error(
        multiomics:::.select_experiments(mae, c(1.5, 2)),
        "must specify names or index"
    )
})

test_that("SCE experiment selection includes main and preserves order", {
    sce <- make_selection_sce_for_utils()

    named <- multiomics:::.select_experiments_from_tse(
        sce,
        c("outcomes", "main")
    )
    indexed <- multiomics:::.select_experiments_from_tse(sce, c(2L, 1L))

    expect_identical(names(named), c("outcomes", "main"))
    expect_identical(names(indexed), c("main", "outcomes"))
    expect_identical(
        SummarizedExperiment::assay(named[["main"]], "counts"),
        SummarizedExperiment::assay(sce, "counts")
    )
    expect_identical(
        SummarizedExperiment::assay(named[["outcomes"]], "counts"),
        SummarizedExperiment::assay(
            SingleCellExperiment::altExp(sce, "outcomes"),
            "counts"
        )
    )
})

test_that("SCE experiment selection rejects invalid selections", {
    sce <- make_selection_sce_for_utils()

    expect_error(
        multiomics:::.select_experiments_from_tse(
            sce,
            c("missing", "main")
        ),
        "must specify names or index"
    )
    expect_error(
        multiomics:::.select_experiments_from_tse(sce, c(0L, 1L)),
        "must specify names or index"
    )
    expect_error(
        multiomics:::.select_experiments_from_tse(sce, c(1.5, 2)),
        "must specify names or index"
    )
})

test_that("MAE assay selection retains only requested assays", {
    mae <- make_selection_mae_for_utils()

    result <- multiomics:::.select_assays_from_mae(
        mae,
        c("scaled", "counts")
    )

    expect_identical(SummarizedExperiment::assayNames(result[[1L]]), "scaled")
    expect_identical(SummarizedExperiment::assayNames(result[[2L]]), "counts")
    expect_identical(
        SummarizedExperiment::assay(result[[1L]], "scaled"),
        SummarizedExperiment::assay(mae[[1L]], "scaled")
    )
})

test_that("MAE assay selection rejects missing assays", {
    mae <- make_selection_mae_for_utils()

    expect_error(
        multiomics:::.select_assays_from_mae(
            mae,
            c("missing", "counts")
        ),
        "Cannot find assay"
    )
})

test_that("SCE assay selection retains only requested assays", {
    sce <- make_selection_sce_for_utils()
    selected <- multiomics:::.select_experiments_from_tse(
        sce,
        c("main", "outcomes")
    )

    result <- multiomics:::.select_assays_from_tse_list(
        selected,
        c("scaled", "counts")
    )

    expect_identical(names(result), c("main", "outcomes"))
    expect_identical(SummarizedExperiment::assayNames(result[[1L]]), "scaled")
    expect_identical(SummarizedExperiment::assayNames(result[[2L]]), "counts")
    expect_identical(
        SummarizedExperiment::assay(result[[1L]], "scaled"),
        SummarizedExperiment::assay(sce, "scaled")
    )
})

test_that("SCE assay selection rejects missing assays", {
    sce <- make_selection_sce_for_utils()
    selected <- multiomics:::.select_experiments_from_tse(
        sce,
        c("main", "outcomes")
    )

    expect_error(
        multiomics:::.select_assays_from_tse_list(
            selected,
            c("missing", "counts")
        ),
        "Cannot find assay"
    )
})
