# The success path cannot be pinned down: passing validation is the absence of an
# error, so an empty function body satisfies any assertion that can be written
# here -- expect_no_error and expect_null both pass on one. expect_no_error is
# used anyway because it states what the test is actually for; the return value
# is incidental and no caller reads it.
#
# The failure-path test is therefore what pins the behaviour, and it anchors the
# whole message with ^...$ so that both halves are checked: the rejected
# reduction name and the listing of what the object actually has.

test_that("ValidateReduction passes silently when the reduction exists", {
  expect_no_error(ValidateReduction(pbmc_small, "pca"))
})

test_that("ValidateReduction stops and lists the available reductions", {
  # pbmc_small_v5 ships only the tsne reduction, no umap
  expect_error(ValidateReduction(pbmc_small_v5, "umap"),
               "^Reduction not found: 'umap'\\. Available: tsne$")
})
