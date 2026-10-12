test_that("getPairwiseAssociation uses TreeSE colData for group.by", {
    samples <- paste0("sample", seq_len(8L))
    first <- rbind(
        feature_a = c(1, 2, 4, 3, 5, 7, 6, 8),
        feature_b = c(8, 3, 6, 1, 7, 2, 5, 4)
    )
    second <- rbind(
        feature_c = c(2, 5, 1, 7, 3, 8, 4, 6),
        feature_d = c(8, 4, 6, 2, 1, 7, 3, 5)
    )
    colnames(first) <- samples
    colnames(second) <- samples
    col_data <- data.frame(
        group = rep(c("x", "y"), 4L),
        row.names = samples
    )
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = first),
        colData = S4Vectors::DataFrame(col_data)
    )
    SingleCellExperiment::altExp(sce, "second") <-
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = second)
        )

    result <- getPairwiseAssociation(
        sce,
        experiments = c("main", "second"),
        assay.types = c("counts", "counts"),
        group.by = "group"
    )
    expected <- getPairwiseAssociation(
        list(first, second),
        col.data = col_data,
        group.by = "group"
    )

    expect_equal(result, expected)
})
