context("test-classifiers.R")

test_that("da_tss returns all valid hits", {
    result <- da_tss(tss, genes)
    ## Count hits.
    expect_equal(length(result), 6)
    ## Check distances.
    distance <- rep(c(0L, 100L, 300L), 2)
    expect_equal(mcols(result)$distance, distance)
    ## Check from-to (# nolint) mapping.
    expected <- Hits(from = idx_da_tss,
                     to = rep(1:2, each = 3),
                     nLnode = length(tss),
                     nRnode = length(genes),
                     distance = distance)
    expect_equal(expected, result)
})

test_that("da_tss validates that genes do not overlap", {
    expect_error(da_tss(tss, genes_overlapping), "overlap")
})

test_that("da_tss validates min parameter", {
    expect_error(da_tss(tss, genes, min = 400, max = 300), "min")
})

test_that("annotate classifies all features", {
    result <- annotate(tss, genes)
    is_da_tss <- mcols(result)$class %in% "daTSS"
    expect_equal(sum(is_da_tss), 6)
    expect_equal(which(is_da_tss), idx_da_tss)
})
