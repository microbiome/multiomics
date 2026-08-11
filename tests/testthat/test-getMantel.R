make_mantel_inputs <- function() {
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

make_mantel_mae <- function(inputs) {
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            first = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$first)
            ),
            second = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$second)
            )
        )
    )
}

test_that("getMantel matches vegan::mantel for list input", {
    inputs <- make_mantel_inputs()

    set.seed(123)
    result <- getMantel(
        list(inputs$first, inputs$second),
        dist.methods = c("euclidean", "euclidean"),
        method = "spearman",
        npermutations = 99L
    )

    set.seed(123)
    expected <- vegan::mantel(
        vegan::vegdist(t(inputs$first), method = "euclidean"),
        vegan::vegdist(t(inputs$second), method = "euclidean"),
        method = "spearman",
        permutations = 99L
    )

    expect_s3_class(result, "mantel")
    expect_equal(result$statistic, expected$statistic)
    expect_equal(result$signif, expected$signif)
    expect_equal(result$perm, expected$perm)
})

test_that("getMantel preserves results through a MultiAssayExperiment", {
    inputs <- make_mantel_inputs()
    mae <- make_mantel_mae(inputs)

    set.seed(123)
    direct_result <- getMantel(
        list(inputs$first, inputs$second),
        dist.methods = c("euclidean", "euclidean"),
        method = "spearman",
        npermutations = 99L
    )

    set.seed(123)
    mae_result <- getMantel(
        mae,
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts"),
        dist.methods = c("euclidean", "euclidean"),
        method = "spearman",
        npermutations = 99L
    )

    expect_s3_class(mae_result, "mantel")
    expect_equal(mae_result$statistic, direct_result$statistic)
    expect_equal(mae_result$signif, direct_result$signif)
    expect_equal(mae_result$perm, direct_result$perm)
})

test_that("getMantel validates list inputs and arguments", {
    inputs <- make_mantel_inputs()
    valid_args <- list(
        dist.methods = c("euclidean", "euclidean"),
        npermutations = 9L
    )

    expect_error(
        do.call(getMantel, c(list(inputs$first), valid_args)),
        "must be list"
    )
    expect_error(
        do.call(getMantel, c(list(list(inputs$first)), valid_args)),
        "length must be 2"
    )
    expect_error(
        do.call(
            getMantel,
            c(list(list(as.data.frame(inputs$first), inputs$second)), valid_args)
        ),
        "include matrices"
    )
    expect_error(
        do.call(
            getMantel,
            c(
                list(list(inputs$first, inputs$second)),
                list(dist.methods = "euclidean")
            )
        ),
        "dist.methods"
    )
    expect_error(
        do.call(
            getMantel,
            c(
                list(list(inputs$first, inputs$second)),
                list(dist.methods = c("not-a-distance", "euclidean"))
            )
        )
    )
    expect_error(
        getMantel(
            list(inputs$first, inputs$second),
            dist.methods = c("euclidean", "euclidean"),
            method = 1,
            npermutations = 9L
        ),
        "method"
    )
    expect_error(
        getMantel(
            list(inputs$first, inputs$second),
            dist.methods = c("euclidean", "euclidean"),
            npermutations = 0L
        ),
        "npermutations"
    )
    expect_error(
        getMantel(
            list(inputs$first, inputs$second),
            dist.methods = c("euclidean", "euclidean"),
            npermutations = 3.5
        ),
        "npermutations"
    )
    expect_error(
        getMantel(
            list(inputs$first, inputs$second),
            dist.methods = c("euclidean", "euclidean"),
            npermutations = 9L,
            na.rm = "FALSE"
        ),
        "na.rm"
    )
    expect_error(
        getMantel(
            list(inputs$first, inputs$second),
            dist.methods = c("euclidean", "euclidean"),
            npermutations = 9L,
            group.by = "stratum"
        ),
        "group.by"
    )

    mismatched <- inputs$second
    colnames(mismatched)[[1L]] <- "not-a-sample"
    expect_error(
        getMantel(
            list(inputs$first, mismatched),
            dist.methods = c("euclidean", "euclidean"),
            npermutations = 9L
        ),
        "Sample names"
    )
})
