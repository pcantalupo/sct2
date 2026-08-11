# The --outfile guards in inst/scripts/seurat_update_object.R run in the options
# block, before pacman::p_load(sct2, Seurat) and ReadSeurat(). These tests drive
# the script as a subprocess and assert both the message and that it stopped
# before print(opts), the first statement after the guards -- a guard that moved
# below the object load would still error, but would print the options first.

run_update_object = function(outfile) {
  script = system.file("scripts", "seurat_update_object.R", package = "sct2")
  skip_if(script == "", "seurat_update_object.R not found")
  # The guards fire before the object is read, so the input only has to exist.
  input = tempfile(fileext = ".qs2")
  file.create(input)
  on.exit(unlink(input))
  # system2() warns on a non-zero exit; the status is asserted by the caller.
  suppressWarnings(
    system2("Rscript", c(script, "--seurat", input, "--outfile", outfile),
            stdout = TRUE, stderr = TRUE)
  )
}

reached_options_print = function(out) {
  any(grepl("^\\$seurat", out))
}

test_that("seurat_update_object.R rejects an --outfile with an unsupported extension", {
  out = run_update_object(file.path(tempdir(), "obj.qs"))
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"),
               "--outfile must end in .rds or .qs2", fixed = TRUE)
  expect_false(reached_options_print(out))
})

test_that("seurat_update_object.R rejects an --outfile in a directory that does not exist", {
  out = run_update_object(file.path(tempdir(), "sct2-no-such-dir", "obj.qs2"))
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"),
               "--outfile directory does not exist", fixed = TRUE)
  expect_false(reached_options_print(out))
})

test_that("seurat_update_object.R rejects an --outfile in a directory that is not writable", {
  readonly = file.path(tempdir(), "sct2-readonly")
  dir.create(readonly)
  Sys.chmod(readonly, "500")
  on.exit({
    Sys.chmod(readonly, "700")
    unlink(readonly, recursive = TRUE)
  })
  skip_if(file.access(readonly, mode = 2) == 0, "directory is still writable (running as root?)")
  out = run_update_object(file.path(readonly, "obj.qs2"))
  expect_equal(attr(out, "status"), 1L)
  expect_match(paste(out, collapse = "\n"),
               "--outfile directory is not writable", fixed = TRUE)
  expect_false(reached_options_print(out))
})

test_that("seurat_update_object.R accepts a writable directory and an uppercase .RDS", {
  # The script loads the installed sct2 once validation passes; it then fails on
  # the dummy input, which is past the point these tests cover and is not asserted.
  skip_if_not_installed("sct2")
  out = run_update_object(tempfile(fileext = ".RDS"))
  expect_true(reached_options_print(out))
})
