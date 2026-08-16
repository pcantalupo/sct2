# The success path cannot be pinned down: passing validation is the absence of an
# error, so an empty function body satisfies any assertion that can be written
# here -- expect_no_error and expect_null both pass on one. expect_no_error is
# used anyway because it states what the test is actually for; the return value
# is incidental and no caller reads it.
#
# The failure-path test is therefore what pins the behaviour, and it anchors the
# whole message with ^...$. An unanchored regex matches a substring, so
# "nope1, nope2" would still pass if setdiff() regressed to listing every column
# it was given -- the message would read "orig.ident, nope1, nope2" and match.

test_that("ValidateMetadataCols passes silently when all columns exist", {
  expect_no_error(ValidateMetadataCols(pbmc_small, c("orig.ident", "groups")))
})

test_that("ValidateMetadataCols stops and lists every missing column", {
  expect_error(ValidateMetadataCols(pbmc_small, c("orig.ident", "nope1", "nope2")),
               "^Metadata column\\(s\\) not found in Seurat object: nope1, nope2$")
})
