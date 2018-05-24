context("test-checks.R")

test_that("check_overlap requires disjoint GRanges", {
    expect_true(check_disjoint(genes))
    expect_error(check_disjoint(genes_overlapping), "overlap")
})

test_that("check_overlap uses non-standard evaluation", {
    expect_error(check_disjoint(genes_overlapping), "genes_overlapping")
})

test_that("check_lt requires strictly less numeric values", {
    expect_true(check_lt(4000, 4001))
    expect_true(check_lt(4000.0, 4001.0))
    expect_error(check_lt(4000, 4000), "less")
    expect_error(check_lt(4000.0, 4000.0), "less")
    expect_error(check_lt(4000, 3999), "less")
})

test_that("check_lt uses non-standard evaluation", {
    expect_error(check_lt(5001, 4001), "5001")
    expect_error(check_lt(5001, 4001), "4001")
})
