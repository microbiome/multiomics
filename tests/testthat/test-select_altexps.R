test_that("numeric altExp selection accepts valid indices", {
    samples <- paste0("sample", seq_len(4L))

    make_sce <- function(main, alternative) {
        x <- SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = main)
        )
        altExp(x, "alternative") <-
            SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = alternative)
            )
        x
    }

    first_alt <- matrix(
        c(1, 2, 4, 8, 3, 5, 7, 11),
        nrow = 2L,
        dimnames = list(c("first_a", "first_b"), samples)
    )
    second_alt <- matrix(
        c(2, 1, 5, 9, 4, 6, 10, 12),
        nrow = 2L,
        dimnames = list(c("second_a", "second_b"), samples)
    )
    first <- make_sce(
        matrix(seq_len(12L), nrow = 3L, dimnames = list(NULL, samples)),
        first_alt
    )
    second <- make_sce(
        matrix(seq_len(12L) + 12L, nrow = 3L, dimnames = list(NULL, samples)),
        second_alt
    )
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            first = first,
            second = second
        )
    )

    result <- getRDA(
        mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts"),
        altexps = c(1L, 1L)
    )

    expect_s3_class(result, "rda")
    expect_equal(ncol(result$Y), nrow(first_alt))
})
