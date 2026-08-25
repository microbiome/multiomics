test_that("plotLASSO returns a coefficient heatmap", {
    skip_if_not_installed("ggplot2")
    coefficients <- matrix(
        c(1, NA, -2, 3, 4, NA),
        nrow = 2L,
        dimnames = list(
            c("taxon_a", "taxon_b"),
            c("metabolite_a", "metabolite_b", "metabolite_c")
        )
    )

    result <- plotLASSO(coefficients)

    expect_s3_class(result, "ggplot")
    expect_identical(
        names(result$data),
        c("Taxa", "Metabolite", "Coefficient")
    )
    expect_equal(
        result$data$Taxa,
        rep(c("taxon_a", "taxon_b"), times = 3L)
    )
    expect_equal(
        result$data$Metabolite,
        rep(
            c("metabolite_a", "metabolite_b", "metabolite_c"),
            each = 2L
        )
    )
    expect_equal(result$data$Coefficient, as.vector(coefficients))
})

test_that("plotLASSO validates its input", {
    expect_error(
        plotLASSO(data.frame(coefficient = 1:2)),
        "'x' must be a numeric matrix"
    )
    expect_error(
        plotLASSO(matrix(as.character(1:4), nrow = 2L)),
        "'x' must be a numeric matrix"
    )
})
