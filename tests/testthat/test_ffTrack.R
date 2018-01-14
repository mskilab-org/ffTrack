library(ffTrack)

context('ffTrack operations')


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
    ## show
    ## 
    ## size
    expect_equal(size(testff), 0.001252)
    ## vmode 
    expect_match(vmode(testff), 'boolean')
    ## len
    expect_equal(len(testff), 10001)
    ## levels
    ## set_levels
    ## ranges
    expect_true(is(ranges(testff), 'GRanges'))
    expect_equal(width(ranges(testff)), 10001)
    ## filename
    ## cp ## then check it exists as done
    ## seqlengths
    ## seqinfo
    ## seqlevels
    ## del 
    ## writeable_status
    ## writeable
    ## mv ## then check it exists as done


})














