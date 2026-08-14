# inst/scripts/seurat_downsample.R resolves the default --outfile (and its
# _ds<tag>) before pacman::p_load(sct2, Seurat) and ReadSeurat(), so the tag can
# be observed cheaply: pre-create the file the tag is expected to produce and
# read the name back out of the "--outfile already exists" guard. The input only
# has to exist -- these runs stop before it is ever read.
#
# The tag scheme exists to keep a fraction and a count that share digits (0.3 and
# 30) from resolving to the same output file: a fraction tags as a percentage and
# a count tags in thousands with a trailing "k".

script_path = function() {
  script = system.file("scripts", "seurat_downsample.R", package = "sct2")
  # A failed lookup is a broken install, not a reason to skip silently.
  if (!nzchar(script)) {
    stop("seurat_downsample.R not found in the installed package")
  }
  return(script)
}

# Runs the script in a scratch directory holding an empty seurat.qs2 plus the
# expected output name, and returns the combined stdout/stderr.
run_downsample = function(downsample, expected_outfile) {
  dir = tempfile()
  dir.create(dir)
  input = file.path(dir, "seurat.qs2")
  file.create(input)
  file.create(file.path(dir, expected_outfile))
  old = setwd(dir)
  on.exit({
    setwd(old)
    unlink(dir, recursive = TRUE)
  })
  # system2() warns on a non-zero exit; the status is asserted by the caller.
  out = suppressWarnings(
    system2("Rscript", c(script_path(), "--seurat", input,
                         "--downsample", downsample),
            stdout = TRUE, stderr = TRUE)
  )
  return(out)
}

# The clobber guard is the probe every tag test below reads the name from, so it
# is anchored first: if the guard stopped naming the file, those tests would pass
# while checking nothing.
expect_tag = function(downsample, expected_outfile) {
  out = run_downsample(downsample, expected_outfile)
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"),
               paste0("--outfile already exists: ", expected_outfile),
               fixed = TRUE)
}

test_that("a fraction tags as a percentage", {
  expect_tag(0.05, "seurat_ds05.qs2")
  expect_tag(0.3, "seurat_ds30.qs2")
})

test_that("a count tags in thousands with a trailing k", {
  expect_tag(30, "seurat_ds0.03k.qs2")
  expect_tag(1500, "seurat_ds1.5k.qs2")
  expect_tag(50000, "seurat_ds50k.qs2")
})

test_that("a fraction and a count with the same digits do not collide", {
  # The whole point of the "k": without it both of these resolve to seurat_ds30.
  expect_tag(0.3, "seurat_ds30.qs2")
  expect_tag(30, "seurat_ds0.03k.qs2")
})

test_that("--downsample 1 tags as the whole object, not one cell", {
  # 1 is a fraction, so DownsampleObject() returns every cell; tagging it ds1
  # would name a full-size object as if a single cell had been kept.
  expect_tag(1, "seurat_ds100.qs2")
})

# The object-size check is the one guard that cannot fire before the load, so
# this runs the script against a real object rather than an empty placeholder.
test_that("seurat_downsample.R rejects a count larger than the object", {
  dir = tempfile()
  dir.create(dir)
  input = file.path(dir, "seurat.qs2")
  WriteSeurat(pbmc_small, input)
  old = setwd(dir)
  on.exit({
    setwd(old)
    unlink(dir, recursive = TRUE)
  })
  out = suppressWarnings(
    system2("Rscript", c(script_path(), "--seurat", input, "--downsample", 500),
            stdout = TRUE, stderr = TRUE)
  )
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"),
               paste0("--downsample 500 exceeds the object size ", ncol(pbmc_small)),
               fixed = TRUE)
  # The point of the guard: no misleadingly-named full-size object is left behind.
  expect_equal(list.files(dir, pattern = "_ds"), character(0))
})

test_that("seurat_downsample.R rejects a --downsample that is not positive", {
  out = run_downsample(0, "unused.qs2")
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"),
               "--downsample must be > 0", fixed = TRUE)
  # The guard runs in the options block; reaching print(opts) would mean it had
  # moved below the object load.
  expect_false(any(grepl("^\\$seurat", out)))
})
