test_that("getPairwiseAssociation aligns samples of list input", {
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

    result <- getPairwiseAssociation(list(first, second[, rev(samples)]))
    expected <- getPairwiseAssociation(list(first, second))

    expect_equal(result, expected)
})
