make_mirkat_inputs <- function() {
    samples <- paste0("sample", seq_len(8L))
    microbiome <- rbind(
        feature_a = c(1, 2, 4, 3, 5, 7, 6, 8),
        feature_b = c(8, 3, 6, 1, 7, 2, 5, 4),
        feature_c = c(2, 5, 1, 6, 3, 8, 4, 7)
    )
    outcomes <- rbind(
        outcome_a = c(2, 5, 1, 7, 3, 8, 4, 6),
        outcome_b = c(8, 4, 6, 2, 1, 7, 3, 5)
    )
    colnames(microbiome) <- samples
    colnames(outcomes) <- samples
    list(microbiome = microbiome, outcomes = outcomes)
}

make_mirkat_mae <- function(inputs) {
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            microbiome = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$microbiome)
            ),
            outcomes = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$outcomes)
            )
        )
    )
}

make_mirkat_sce <- function(inputs) {
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = inputs$microbiome)
    )
    SingleCellExperiment::altExp(sce, "outcomes") <-
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = inputs$outcomes)
        )
    sce
}

mirkat_reference <- function(inputs, nperm = 9L) {
    dist.methods <- c("euclidean", "manhattan", "canberra")
    kernels <- lapply(dist.methods, function(method) {
        distances <- vegan::vegdist(
            t(inputs$microbiome),
            method = method
        )
        MiRKAT::D2K(as.matrix(distances))
    })
    names(kernels) <- dist.methods

    rows <- lapply(seq_len(nrow(inputs$outcomes)), function(i) {
        upstream <- MiRKAT::MiRKAT(
            y = inputs$outcomes[i, ],
            Ks = kernels,
            method = "davies",
            omnibus = "cauchy",
            out_type = "C",
            nperm = nperm,
            returnKRV = TRUE,
            returnR2 = TRUE
        )
        values <- unlist(upstream)
        names(values) <- tolower(gsub("\\.", "_", names(values)))
        values
    })
    names(rows) <- rownames(inputs$outcomes)
    expected <- as.data.frame(do.call(rbind, rows))
    expected[["omnibus_p_adj"]] <- stats::p.adjust(
        expected[["omnibus_p"]],
        method = "fdr"
    )
    expected
}

test_that("getMiRKAT matches MiRKAT for list input", {
    inputs <- make_mirkat_inputs()

    result <- getMiRKAT(list(inputs$microbiome, inputs$outcomes), nperm = 9L)
    expected <- mirkat_reference(inputs)

    expect_s3_class(result, "data.frame")
    expect_identical(rownames(result), rownames(expected))
    expect_identical(colnames(result), colnames(expected))
    expect_equal(result, expected)
})

test_that("getMiRKAT preserves results through MultiAssayExperiment and SingleCellExperiment", {
    inputs <- make_mirkat_inputs()
    direct_result <- getMiRKAT(
        list(inputs$microbiome, inputs$outcomes)
    )

    mae_result <- getMiRKAT(
        make_mirkat_mae(inputs),
        experiments = c(1L, 2L),
        assay.types = c("counts", "counts")
    )
    sce_result <- getMiRKAT(
        make_mirkat_sce(inputs),
        experiments = c("main", "outcomes"),
        assay.types = c("counts", "counts")
    )

    expect_equal(mae_result, direct_result)
    expect_equal(sce_result, direct_result)
})

test_that("getMiRKAT validates list inputs and arguments", {
    inputs <- make_mirkat_inputs()
    valid_args <- list(
        dist.methods = "euclidean",
        nperm = 9L
    )

    expect_error(
        do.call(getMiRKAT, c(list(inputs$microbiome), valid_args)),
        "must be list"
    )
    expect_error(
        do.call(
            getMiRKAT,
            c(list(list(inputs$microbiome)), valid_args)
        ),
        "length must be 2"
    )
    expect_error(
        do.call(
            getMiRKAT,
            c(
                list(list(as.data.frame(inputs$microbiome), inputs$outcomes)),
                valid_args
            )
        ),
        "include matrices"
    )
    expect_error(
        do.call(
            getMiRKAT,
            c(
                list(list(inputs$microbiome, inputs$outcomes)),
                list(dist.methods = 1)
            )
        ),
        "dist.methods"
    )
    expect_error(
        getMiRKAT(
            list(inputs$microbiome, inputs$outcomes),
            dist.methods = "not-a-distance",
            nperm = 9L
        )
    )
    expect_error(
        getMiRKAT(
            list(inputs$microbiome, inputs$outcomes),
            dist.methods = "euclidean",
            method = 1,
            nperm = 9L
        ),
        "method"
    )
    expect_error(
        getMiRKAT(
            list(inputs$microbiome, inputs$outcomes),
            dist.methods = "euclidean",
            npermutations = 3.5
        ),
        "npermutations"
    )
    expect_error(
        getMiRKAT(
            list(inputs$microbiome, inputs$outcomes),
            dist.methods = "euclidean",
            omnibus = 1,
            nperm = 9L
        ),
        "omnibus"
    )
    expect_error(
        getMiRKAT(
            list(inputs$microbiome, inputs$outcomes),
            dist.methods = "euclidean",
            p.adjust.method = "not-a-method",
            nperm = 9L
        ),
        "p.adjust.method"
    )

    mismatched <- inputs$outcomes
    colnames(mismatched)[[1L]] <- "not-a-sample"
    expect_error(
        do.call(
            getMiRKAT,
            c(list(list(inputs$microbiome, mismatched)), valid_args)
        ),
        "Sample names"
    )

    mae <- make_mirkat_mae(inputs)
    expect_error(
        getMiRKAT(
            mae,
            experiments = c(1L, 2L),
            assay.types = "counts",
            nperm = 9L
        ),
        "assay.types"
    )
    expect_error(
        getMiRKAT(
            mae,
            experiments = c(1L, 2L),
            assay.types = c("missing", "counts"),
            nperm = 9L
        ),
        "Cannot find assay"
    )
})
