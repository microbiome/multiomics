make_procrustes_inputs <- function() {
    samples <- paste0("sample", seq_len(8L))
    first <- rbind(
        feature_a = seq_len(8L),
        feature_b = c(2, 1, 4, 3, 6, 5, 8, 7),
        feature_c = c(8, 3, 6, 1, 7, 2, 5, 4)
    )
    second <- rbind(
        feature_d = c(3, 1, 4, 2, 5, 7, 6, 8),
        feature_e = c(8, 6, 4, 2, 1, 3, 5, 7),
        feature_f = c(2, 5, 1, 6, 3, 8, 4, 7)
    )
    colnames(first) <- samples
    colnames(second) <- samples
    list(first = first, second = second)
}

test_that("getProcrustes matches vegan::protest for list input", {
    inputs <- make_procrustes_inputs()

    set.seed(123)
    result <- getProcrustes(
        list(inputs$first, inputs$second),
        dist.methods = c("euclidean", "euclidean"),
        npermutations = 99L
    )

    set.seed(123)
    expected <- vegan::protest(
        stats::cmdscale(vegan::vegdist(t(inputs$first), "euclidean")),
        stats::cmdscale(vegan::vegdist(t(inputs$second), "euclidean")),
        permutations = 99L
    )

    expect_s3_class(result, "protest")
    expect_equal(result$t0, expected$t0)
    expect_equal(result$signif, expected$signif)
})

test_that("getProcrustes preserves results through a MultiAssayExperiment", {
    inputs <- make_procrustes_inputs()
    mae <- MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            first = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$first)
            ),
            second = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$second)
            )
        )
    )

    set.seed(123)
    direct_result <- getProcrustes(
        list(inputs$first, inputs$second),
        dist.methods = c("euclidean", "euclidean"),
        npermutations = 99L
    )

    set.seed(123)
    mae_result <- getProcrustes(
        mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts"),
        dist.methods = c("euclidean", "euclidean"),
        npermutations = 99L
    )

    expect_s3_class(mae_result, "protest")
    expect_equal(mae_result$t0, direct_result$t0)
    expect_equal(mae_result$signif, direct_result$signif)
})

test_that("getProcrustes checks input", {
    inputs <- make_procrustes_inputs()
    expect_error(getProcrustes(list(inputs$first)))
    expect_error(
        getProcrustes(list(inputs$first, as.data.frame(inputs$second)))
    )
    second <- inputs$second
    colnames(second)[[1L]] <- "other"
    expect_error(getProcrustes(list(inputs$first, second)))
})
