
library(ffTrack)
library(testthat)
library(BSgenome)
library(gUtils)

context('ffTrack operations')



### ffTrack
### --- initalize    XX
### --- show         XX  
### --- size         XX
### --- vmode        XX
### --- len          XX
### --- levels       XX
### --- set_levels   
### --- ffranges
### --- filename
### --- cp
### --- ffseqlengths
### --- seqinfo
### --- seqlevels
### --- del
### --- writeable_status
### --- [
### --- writeable
### --- mv
### --- [<-
### --- get_seq


### get_seq
### bw2fft
### wig2fft
### seq2fft
### fftab

#system("mkdir tmp")

##system("mkdir tmp2")



test_that('ffTrack', {
    
    expect_error(ffTrack(vmode = 'error')) ## Error: Incorrect argument "vmode". Allowable modes are boolean, byte, character, complex, double, integer, logical, nibble, quad, raw, short, single, ubyte, ushort
    expect_error(ffTrack(gr = GRanges()))  ## Error: Trying to create ffTrack with empty GRanges
    ## check 'vmode'
    ## boolean (1 bit logical)
    ## logical (2 bit logical + NA)
    ## quad (2 bit unsigned integer without NA)
    ## nibble (4 bit unsigned integer without NA)
    ## byte (8 bit signed integer with NA)
    ## ubyte (8 bit unsigned integer without NA)
    ## short (16 bit signed integer with NA)
    ## ushort (16 bit unsigned integer without NA)
    ## integer (32 bit signed integer with NA)
    ## single (32 bit float)
    ## double (64 bit float)
    ## raw (8 bit unsigned char)
    gr = GRanges('1:10000-20000')
    expect_error(ffTrack(gr, file.name = 'test.boolean.rds', overwrite = TRUE, vmode = 'boolean'), NA)  ## trick to check code doesn't throw error
    expect_error(ffTrack(gr, file.name = 'test.logical.rds', overwrite = TRUE, vmode = 'logical'), NA)
    expect_error(ffTrack(gr, file.name = 'test.quad.rds', overwrite = TRUE, vmode = 'quad'), NA)
    expect_error(ffTrack(gr, file.name = 'test.nibble.rds', overwrite = TRUE, vmode = 'nibble'), NA)
    expect_error(ffTrack(gr, file.name = 'test.byte.rds', overwrite = TRUE, vmode = 'byte'), NA)
    expect_error(ffTrack(gr, file.name = 'test.ubyte.rds', overwrite = TRUE, vmode = 'ubyte'), NA)
    expect_error(ffTrack(gr, file.name = 'test.short.rds', overwrite = TRUE, vmode = 'short'), NA)
    expect_error(ffTrack(gr, file.name = 'test.ushort.rds', overwrite = TRUE, vmode = 'ushort'), NA)
    expect_error(ffTrack(gr, file.name = 'test.integer.rds', overwrite = TRUE, vmode = 'integer'), NA)
    expect_error(ffTrack(gr, file.name = 'test.single.rds', overwrite = TRUE, vmode = 'single'), NA)
    expect_error(ffTrack(gr, file.name = 'test.double.rds', overwrite = TRUE, vmode = 'double'), NA)
    expect_error(ffTrack(gr, file.name = 'test.raw.rds', overwrite = TRUE, vmode = 'raw'), NA)
    expect_error(ffTrack(gr, file.name = 'test.complex.rds', overwrite = TRUE, vmode = 'complex'))  ## Message: vmode 'complex' not implemented
    expect_error(ffTrack(gr, file.name = 'test.character.rds', overwrite = TRUE, vmode = 'character')) ## Message: vmode 'character' not implemented
    ## test ffTrack methods
    testff = ffTrack(gr, file.name = 'test.boolean.rds', overwrite = TRUE, vmode = 'boolean')
    ## if (!(vmode[1] %in% MODES)){
    expect_error(ffTrack(gr, file.name = 'test.failure.rds', overwrite = TRUE, vmode = 'failure'))
    ## if (length(gr) == 0 | !any(width(gr) > 0))
    expect_error(ffTrack(GRanges(), file.name = 'test.empty.rds', overwrite = TRUE, vmode = 'boolean'))
    ## if (!grepl('\\.rds$',  file.name) & !grepl('\\.RDS$',  file.name)){
    expect_error(ffTrack(gr, file.name = 'test.empty', overwrite = TRUE, vmode = 'boolean'), NA)
    ## if ((file.exists(file.name) | file.exists(ff.filename)) & !overwrite){
    expect_error(ffTrack(gr, file.name = 'test.boolean.rds', overwrite = FALSE, vmode = 'boolean'))
    ## ISSUE
    ## 'complex', 'character' not implemented
    ##
    ## > (ffTrack(gr, file.name = 'test.character.rds', overwrite = TRUE, vmode = 'character'))
    ## Error in ff(default.val, length = pmin(len, .Object@.blocksize), vmode = .Object@.vmode,  : 
    ##   vmode 'character' not implemented
    ## 
    ## if (verbose){
    foobar = ffTrack(gr, file.name = 'test.boolean.rds', overwrite = TRUE, vmode = 'boolean', verbose=TRUE)
    ## show
    expect_equal(basename(filename(foobar))[1], 'test.boolean.ffdata')
    expect_equal(basename(filename(foobar))[2], 'test.boolean.rds')
    ## 
    ## size
    expect_equal(ffTrack::size(testff), 0.001252)
    ## vmode 
    expect_match(ffTrack::vmode(testff), 'boolean')
    ## length
    expect_equal(length(foobar), 10001)
    ## levels
    ##
    ## expect_equal(as.logical(levels(foobar)), NA)
    ## setMethod('levels<-', 'ffTrack', function(x, value){
    ## 
    expect_error(levels(foobar), NA) ## checks runs
    levels(foobar) = 'chr2'
    ## expect_match(levels(foobar), 'chr2')
    ## ffranges
    expect_true(is(ffTrack::ranges(testff), 'GRanges'))
    expect_equal(width(ffTrack::ranges(testff)), 10001)
    expect_equal(ffTrack::ranges(testff)$ix.s, 1)
    ## filename
    expect_equal(basename(ffTrack::filename(testff)[1]), 'test.boolean.ffdata')
    expect_equal(basename(ffTrack::filename(testff)[2]), 'test.boolean.rds')    
    
})




## print('/home/travis/build/mskilab/ffTrack/tmp2/:   ')
## print(list.files('/home/travis/build/mskilab/ffTrack/tmp2/'))
## print('/home/travis/build/mskilab/ffTrack/:   ')
## print(list.files("/home/travis/build/mskilab/ffTrack/"))

### [<-




test_that('checking [<- ', {

    gr = GRanges('2:10000-20000')
    testff = ffTrack(gr, file.name = 'test.ffdatapopulator.rds', overwrite = TRUE, vmode = 'boolean')
    expect_error((testff[gr.start(gr)] = 15000), NA)
    ## > testff[gr.start(gr)] = rep(15, width(gr))
    ## Error in .local(x, i, ..., value) : 
    ##   Error: value must be list or vector of same length as GRanges input "i", or if full = TRUE a vector of same length as sum(width(granges))
    ## > testff[gr.start(gr)] = rep(15, width(gr)+1)
    ## Error in .local(x, i, ..., value) : 
    ##   Error: value must be list or vector of same length as GRanges input "i", or if full = TRUE a vector of same length as sum(width(granges))
    ### as the width(gr) > 1
    ## > testff[gr] = runif(sum(width(gr)))
    ## Error in .local(x, i, ..., value) : 
    ##    Error: value must be list or vector of same length as GRanges input "i", or if full = TRUE a vector of same length as sum(width(granges))

})




## 
## but simple example is somethiing like
## myff = ffTrack(granges, path)
## myff[gr.start(granges)] = 10
## test
## myff[granges] = 10
## myff[granges] = rep(10, width(granges))
## but this should fail
## myff[granges] = rep(10, width(granges)+1)
## 
## myff[granges] = runif(width(granges))
## well that's for a length 1 granges
## if it's width > 1
## then something like myff[granges] = runif(sum(width(granges)))
## ie the right hand side needs to be a vector whose length is the same as the summed width of the granges argument
## the granges in the ffTrack instantiator is just the universe of values that you're allowed to populate
## e.g. just exons
## so I think it should error out if you try to populate anything outside of that territory
## myff[granges+10] = 10 should error or at least warn that it can't populate values outsi
## 










##  hg A BSgenome or and ffTrack object with levels = c('A','T','G','C','N')
##  gr GRanges object to define the ranges
##  unlist logical whether to unlist the final output into a single DNAStringSet. Default TRUE
##  mc.cores Optional multicore call. Default 1
##  mc.chunks Optional define how to chunk the multicore call. Default mc.cores
## as.data.table boolean
## verbose Increase verbosity
test_that('get_seq', {

    library(BSgenome)
    hg19 = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
    gr = GRanges('1:10000-20000')
    seqlevelsStyle(gr) = 'UCSC'  ### convert to chr1
    ## default args
    expect_equal(substr(as.character(get_seq(hg19, gr)[[1]]), 1, 10), 'NTAACCCTAA')
    ## unlist
    expect_equal(substr(as.character(get_seq(hg19, gr, unlist=FALSE)[[1]]), 1, 10), 'NTAACCCTAA')    
    ## mc.cores 
    expect_equal(substr(as.character(get_seq(hg19, gr, mc.cores=2)[[1]]), 1, 10), 'NTAACCCTAA')
    ## mc.chunks
    expect_equal(substr(as.character(get_seq(hg19, gr, mc.chunks=2)[[1]]), 1, 10), 'NTAACCCTAA')   
    ## verbose
    expect_equal(substr(as.character(get_seq(hg19, gr, verbose=TRUE)[[1]]), 1, 10), 'NTAACCCTAA')
    ##  if (inherits(gr, 'GRangesList')){
    gr1 = GRanges('chr1:30000-30005')
    gr2 = GRanges('chr2:30000-30005')
    gr3 = GRanges('chr3:30000-30005')
    grlfoo = GRangesList(gr1, gr2, gr3)
    expect_equal(as.character(get_seq(hg19, grlfoo, verbose=TRUE)[[1]]), 'TGGGGA')
    ##  if (is(hg, 'ffTrack')){

})








### bw2fft
test_that('bw2fft', {

    smallbw = 'output.bw'
    foobar = bw2fft(smallbw, mc.cores=2, overwrite = TRUE, verbose=TRUE)
    expect_equal(basename(filename(foobar))[1], 'output.ffdata')
    expect_equal(basename(filename(foobar))[2], 'output.rds')

})







### wig2fft
# test_that('wig2fft', {
#
# })








## seq2fft
## 
## 
##  Creates ffTrack object from BSGenome or FASTA (coming soon) file
## 
##  ## will either convert (1) raw sequence (2) k-nucleotide context centered around base or (3) motifs defined by some dictionary (anchored at first base) into leveled ffTrack (i.e. integer track with populated levels field)
## 
## seq BSGenome object, ffTrack object representing genomic sequence, or (not yet supported) FASTA file
## nnuc how many nucleotides to left and right to enumerate
## dict this should be a character vector or DNAStringSet, overrides nnuc arg if not null
## chrsub whether to sub in / sub out 'chr' when accessing seq file
## neg whether to analyze sequence data on negative strand (i.e. motifs will be analyzed in rev complement)
## region GRanges specifying regions to limit ffTrack computation to (default is whole genome, ie seqnames of BigWig file)
## mc.cores currently mc.cores can only be one (weird mclapply bug when running)
## verbose logical flag
## buffer integer size of how big of a buffer to use when transferring data from BigWig to ffTrack object; number of bases to access at a time
## skip.sweep logical flag (default FALSE) if TRUE will skip the sweep of "region" for the portions that have non-NA values; if TRUE will not sweep for covered region, just make a whole genome file or a file across provided regions
## vmode  character specifyhing vmode to use for encoding (by default double)
## min.gapwidth  minimum gap-width with which to merge reference adjacent intervals, this will mildly increase the file size but reduce the range complexity of the GRanges object; flank (to reduce the range complexity of the ffdata skeleton, but increase file size)
## ffTrack object corresponding to the data in the BigWig file
## 
##
## seq2fft = function(seq, nnuc = 0, dict = NULL, chrsub = TRUE, neg = FALSE, region = NULL, mc.cores = 1, verbose = FALSE,
##     buffer = 1e5, skip.sweep = FALSE, vmode = 'ubyte', min.gapwidth = 1e3)



test_that('seq2fft works', {

    ## default
    hg19 = getBSgenome('BSgenome.Hsapiens.UCSC.hg19')
    gr = GRanges('1:10000-20000')
    ## seq2fft(hg19)  fails
    ## testff = ffTrack(hg19, file.name = 'test.hg19.rds', overwrite = TRUE, vmode = 'boolean')
    
    ###
    ## if (!inherits(seq, 'BSgenome') & !is(seq, 'ffTrack'))
    expect_error(seq2fft(GRanges()))
    ### Error in is(fftpath, "ffTrack") : object 'fftpath' not found
    expect_error(seq2fft(hg19))

})




test_that('fftab works', {

    ## default
    gr = GRanges('1:10000-20000')
    gr2 = GRanges(1, IRanges(c(1,9), c(6,14)), strand=c('+','-'), seqinfo=Seqinfo("1", 25), field=c(1,2))
    testff = ffTrack(gr, file.name = 'test.boolean.rds', overwrite = TRUE, vmode = 'boolean') 
    integer = GRanges(1, IRanges(c(3,7,13), c(5,9,16)), strand=c('+','-','-'), seqinfo=Seqinfo("1", 25), name=c("A","B","C"))
    foo1 = fftab(testff, integer)
    foo2 = fftab(testff, gr2)
    expect_true(is(foo1, 'GRanges'))
    expect_equal(length(foo1), 3)
    expect_equal(foo1$name[3], 'C')
    ## signature
    foosig = fftab(testff, integer, signatures=list(3, 9, 27))
    expect_equal(length(foosig), 3)
    expect_equal(unique(foosig$sig1), 0)
    expect_equal(unique(foosig$sig2), 0)
    expect_equal(unique(foosig$sig3), 0)
    ## FUN
    expect_error(fftab(testff, integer, FUN=abs))   ## Error in FUN(dat, na.rm = na.rm) : 2 arguments passed to 'abs' which requires 1
    func1 = fftab(testff, integer, FUN=mean) 
    func2 = fftab(testff, integer, FUN=rep) 
    expect_equal(length(func1), 3)
    expect_equal(unique(as.logical(func1$score)), NA)  ## nothing calculated; document this a bit more
    expect_equal(length(func2), 3)
    expect_equal(unique(as.logical(func2$score)), NA)
    ## grep 
    ## mc.cores 
    expect_equal(length(fftab(testff, integer, mc.cores=2)), 3)
    expect_equal(unique(fftab(testff, integer, mc.cores=2)$score), 0)
    ## chunksize 
    expect_equal(length(fftab(testff, integer, chunksize=2)), 3)
    expect_equal(unique(fftab(testff, integer, chunksize=2)$score), 0)
    ## verbose 
    fftab(testff, integer, verbose=FALSE)
    expect_equal(length(fftab(testff, integer, verbose=FALSE)), 3)
    expect_equal(unique(fftab(testff, integer, verbose=FALSE)$score), 0)
    ## na.rm 
    footrue = fftab(testff, integer, na.rm=TRUE)
    foofalse = fftab(testff, integer, na.rm=FALSE)
    expect_equal(unique(footrue$score), 0)
    expect_equal(unique(as.logical(foofalse$score)), NA)
    ## errors
    ## if (!is(ff, 'ffTrack') & !is(ff, 'RleList')){
    expect_error(fftab(GRanges(), integer))  ## Error in fftab(GRanges(), intergr) : Error: Input ff should be ffTrack or RleList
    ##  if (length(intervals)==0){
    expect_error(fftab(testff, GRanges())) 
    ##  if (!is.list(signatures)){
    expect_error(fftab(testff, integer, signatures=matrix())) ##  Error: Signatures must be a named list of arbitrary length character or length 1 or 2 numeric vectors

})












