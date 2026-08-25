make_joint_rpca_inputs <- function() {
    set.seed(1)
    samples <- paste0("sample", seq_len(10L))
    first <- matrix(
        rnorm(40L),
        nrow = 4L,
        dimnames = list(paste0("microbe", seq_len(4L)), samples)
    )
    second <- matrix(
        rnorm(30L),
        nrow = 3L,
        dimnames = list(paste0("metabolite", seq_len(3L)), samples)
    )
    list(first = first, second = second, samples = samples)
}

make_joint_rpca_mae <- function(inputs, reorder_second = FALSE) {
    second <- inputs$second
    if (reorder_second) {
        second <- second[, rev(seq_len(ncol(second))), drop = FALSE]
    }
    MultiAssayExperiment::MultiAssayExperiment(
        experiments = MultiAssayExperiment::ExperimentList(
            microbiome = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = inputs$first)
            ),
            metabolites = SingleCellExperiment::SingleCellExperiment(
                assays = list(counts = second)
            )
        )
    )
}

test_that("getJointRPCA returns the mia result structure for MAE input", {
    inputs <- make_joint_rpca_inputs()
    mae <- make_joint_rpca_mae(inputs)

    result <- getJointRPCA(
        mae,
        experiments = c("microbiome", "metabolites"),
        assay.types = c("counts", "counts"),
        test.set = inputs$samples[1:2],
        ncomponents = 2L,
        max.iterations = 1L
    )

    expect_s3_class(result, "JointRPCA")
    expect_equal(dim(result), c(10L, 2L))
    expect_equal(rownames(result), inputs$samples)
    expect_equal(colnames(result), c("PC1", "PC2"))
    expect_true(all(c(
        "rotation", "percentVar", "lower_dim", "distance",
        "reconstruct_error", "n_features"
    ) %in% names(attributes(result))))
    expect_equal(unname(attr(result, "n_features")), c(4, 3))
})

test_that("getJointRPCA matches mia for the migrated behavior", {
    skip_if_not_installed("mia")
    if (!("getJointRPCA" %in% getNamespaceExports("mia"))) {
        skip("installed mia does not export getJointRPCA")
    }
    inputs <- make_joint_rpca_inputs()
    mae <- make_joint_rpca_mae(inputs)
    args <- list(
        experiments = c("microbiome", "metabolites"),
        assay.types = c("counts", "counts"),
        test.set = inputs$samples[1:2],
        ncomponents = 2L,
        max.iterations = 1L
    )

    result <- do.call(getJointRPCA, c(list(x = mae), args))
    reference <- do.call(mia::getJointRPCA, c(list(x = mae), args))

    expect_equal(result, reference)
})

test_that("getJointRPCA selects a test set by default", {
    inputs <- make_joint_rpca_inputs()
    mae <- make_joint_rpca_mae(inputs)

    result <- getJointRPCA(
        mae,
        experiments = c("microbiome", "metabolites"),
        assay.types = c("counts", "counts"),
        ncomponents = 2L,
        max.iterations = 5L,
        test.ratio = 0.2
    )

    expect_s3_class(result, "JointRPCA")
    expect_equal(dim(result), c(10L, 2L))
    expect_equal(rownames(result), inputs$samples)
})

test_that("getJointRPCA aligns shared MAE samples", {
    inputs <- make_joint_rpca_inputs()
    ordered <- make_joint_rpca_mae(inputs)
    reordered <- make_joint_rpca_mae(inputs, reorder_second = TRUE)

    expected <- getJointRPCA(
        ordered,
        experiments = c("microbiome", "metabolites"),
        assay.types = c("counts", "counts"),
        test.set = inputs$samples[1:2],
        ncomponents = 2L,
        max.iterations = 1L
    )
    result <- getJointRPCA(
        reordered,
        experiments = c("microbiome", "metabolites"),
        assay.types = c("counts", "counts"),
        test.set = inputs$samples[1:2],
        ncomponents = 2L,
        max.iterations = 1L
    )

    expect_equal(result, expected)
})

test_that("getJointRPCA supports SingleCellExperiment altExp selection", {
    inputs <- make_joint_rpca_inputs()
    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = inputs$first)
    )
    SingleCellExperiment::altExp(sce, "metabolites") <-
        SingleCellExperiment::SingleCellExperiment(
            assays = list(counts = inputs$second)
        )

    result <- getJointRPCA(
        sce,
        experiments = c("main", "metabolites"),
        assay.types = c("counts", "counts"),
        test.set = inputs$samples[1:2],
        ncomponents = 2L,
        max.iterations = 1L
    )

    expect_s3_class(result, "JointRPCA")
    expect_equal(dim(result), c(10L, 2L))
    expect_equal(rownames(result), inputs$samples)
})

test_that("getJointRPCA validates selection and test-set inputs", {
    inputs <- make_joint_rpca_inputs()
    mae <- make_joint_rpca_mae(inputs)

    expect_error(
        getJointRPCA(
            mae,
            experiments = "microbiome",
            assay.types = "counts"
        ),
        "length"
    )
    expect_error(
        getJointRPCA(
            mae,
            experiments = c("microbiome", "metabolites"),
            assay.types = c("counts", "counts"),
            test.set = 1L
        ),
        "test.set"
    )
    expect_error(
        getJointRPCA(
            mae,
            experiments = c("microbiome", "missing"),
            assay.types = c("counts", "counts")
        ),
        "experiments"
    )
    expect_error(
        getJointRPCA(
            mae,
            experiments = c("microbiome", "metabolites"),
            assay.types = c("counts", "missing"),
            test.set = inputs$samples[1:2]
        ),
        "Cannot find assay"
    )
    expect_error(
        getJointRPCA(
            mae,
            experiments = c("microbiome", "metabolites"),
            assay.types = c("counts", "counts"),
            test.set = inputs$samples[1:2],
            ncomponents = 10L,
            max.iterations = 1L
        ),
        "ropt"
    )
    expect_error(
        getJointRPCA(
            mae,
            experiments = c("microbiome", "metabolites"),
            assay.types = c("counts", "counts"),
            ncomponents = 2L,
            max.iterations = 5L,
            test.ratio = 1
        ),
        "test.ratio"
    )
})
