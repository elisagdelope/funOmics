library(testthat)

# Simulation data for testing
g <- 10000
s <- 20
X <- matrix(abs(rnorm(g * s)), nrow = g, dimnames = list(paste0("g", 1:g), paste0("s", 1:s)))
pathways <- as.list(sample(10:100, size = 100, replace = TRUE))
pathways <- lapply(pathways, function(s, g) paste0("g", sample(1:g, size = s, replace = FALSE)), g)
names(pathways) <- paste0("pathway", seq_along(pathways))

test_that("summarize_pathway_level works as expected", {
    # Test for successful aggregation with default parameters
    result <- summarize_pathway_level(X, pathways)
    expect_true(!is.null(result), "Should return a non-null result.")

    # Test for unsuccessful aggregation due to minsize constraint
    result_short <- summarize_pathway_level(X, pathways, minsize = 101)
    expect_true(is.null(result_short), "Should return NULL due to minsize constraint.")

    # Test for successful aggregation with a dimension-reduction score type (e.g., NMF)
    result_nmf <- summarize_pathway_level(X, pathways, type = "nmf")
    expect_true(!is.null(result_nmf), "Should return a non-null result for nmf aggregation.")

    # Test for successful aggregation with a hypothesis test score type (e.g., ttest)
    result_ttest <- summarize_pathway_level(X, pathways, type = "ttest")
    expect_true(!is.null(result_ttest), "Should return a non-null result for ttest aggregation.")
})

test_dir("tests/")