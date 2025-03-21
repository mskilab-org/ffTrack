#' @import GenomicRanges
#' @import rtracklayer
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevels
#' @importFrom data.table rbindlist data.table setkey :=
#' @importFrom gUtils gr.sub seg2gr gr.stripstrand si2gr rle.query gr.fix gr.chr gr.tile grl.unlist gr.findoverlaps gr.dice hg_seqlengths
#' @importFrom ff ff is.readonly
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevels
"_PACKAGE"