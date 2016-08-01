library(ffTrack)
library(GenomicRanges)
context("ffTrack ops")

gr <- GRanges(1, IRanges(1,10))

test_that("initialize boolean", {
  ff <- ffTrack(gr, file.name="test.boolean.rds", overwrite = TRUE, vmode='boolean')
})
