library(testthat)
library(GenomicRanges)
library(IRanges)

test_that("Standard collapsing works (Keep First strategy)", {
  gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = c(10, 20, 50), width = 1),
    strand = "*"
  )
  names(gr) <- c("cg1_A", "cg1_B", "cg2_A")
  gr$val <- c(1, 2, 3) 
  
  res <- .collapse_gr_to_pfx(gr, collapse_sesame_control_probes = FALSE)
  
  expect_equal(length(res), 2)
  expect_equal(names(res), c("cg1", "cg2"))
  
  # Check "Keep First"
  expect_equal(start(res["cg1"]), 10)
  expect_equal(res["cg1"]$val, 1)
})


test_that("Control probe logic works correctly", {
  gr <- GRanges(
    seqnames = "chr1",
    ranges = IRanges(start = 1:4, width = 1)
  )
  names(gr) <- c("cg1_A", "cg1_B", "ctl_X_1", "ctl_Y_2")
  
  # 1. FALSE (Default)
  res_off <- .collapse_gr_to_pfx(gr, collapse_sesame_control_probes = FALSE)
  expect_equal(names(res_off), c("cg1", "ctl")) 
  expect_equal(length(res_off), 2)
  
  # 2. TRUE (Keep controls)
  res_on <- .collapse_gr_to_pfx(gr, collapse_sesame_control_probes = TRUE)
  expect_true("cg1" %in% names(res_on))
  expect_true("ctl_X_1" %in% names(res_on))
  expect_true("ctl_Y_2" %in% names(res_on))
  expect_false("ctl" %in% names(res_on))
  expect_equal(length(res_on), 3)
})


test_that("Error handling for unnamed inputs", {
  gr_unnamed <- GRanges(seqnames = "chr1", ranges = IRanges(1, 1))
  expect_error(.collapse_gr_to_pfx(gr_unnamed), "must have names")
})


test_that("Metadata columns are preserved", {
  gr <- GRanges(seqnames = "chr1", ranges = IRanges(1, 1))
  names(gr) <- "cg1_A"
  gr$gene <- "GeneX"
  
  res <- .collapse_gr_to_pfx(gr)
  expect_equal(res$gene, "GeneX")
})