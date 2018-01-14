#############################################################################r
## Marcin Imielinski
## The Broad Institute of MIT and Harvard / Cancer program.
## marcin@broadinstitute.org
##
## Weill-Cornell Medical College
## mai9037@med.cornell.edu
##
## New York Genome Center
## mimielinski@nygenome.org
##
## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

setOldClass('ff_vector')

#' @name ffTrack-class
#' @title ffTrack-class
#' @rdname ffTrack-class
#' @description
#'
#' class::ffTrack
#'
#' Class \code{ffTrack} a pointer for rapid GRanges-based access to genomic data on disk.
#'
#' @import GenomicRanges
#' @importFrom ff ff is.readonly close.ff open.ff
#' @importFrom methods setClass setGeneric setMethod setRefClass
#' @importFrom S4Vectors Rle
#' @importFrom Biostrings DNAStringSet
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevels
#' @exportClass ffTrack
#' @author Marcin Imielinski
setClass('ffTrack', representation(.ff = 'ff_vector', ## primary ff object
                                   .ffaux = 'list', ## list of auxiliary ff objects
                                   .length = 'numeric', ## total length of object
                                   .levels = 'vector', ## vector of unique values (in which case stored values interpreted as 0-based factor)
                                   .blocksize = 'numeric', ## maximum addressable blocksize (.Machine$integer.max))
                                   .gr = 'GRanges', .vmode = 'character', .ff.filename = 'character', .rds.filename = 'character'))
##
## default.val = NA, overwrite = FALSE, levels = NULL, verbose = FALSE
setMethod('initialize', 'ffTrack', function(.Object,
                                            gr, ## GRanges of input ranges
                                            file.name, ## file.name should have .rds suffix, if not, one will be appended
                                            overwrite = FALSE, ## whether to overwrite
                                            levels = NULL, ## vector of unique values (only for vmode of integer or integer-like)
                                            default.val = NA,
                                            verbose = FALSE,
                                            vmode = 'double') ## data mode (see above), scalar character
{

    if (missing(gr) | missing(file.name)){
        stop("Error: arguments 'gr' and 'file.name' are both required for 'ffTrack'. Please see documentation for details.")
    }

    MODES = c('boolean', 'byte', 'character', 'complex', 'double', 'integer',
        'logical', 'nibble', 'quad', 'raw', 'short', 'single', 'ubyte', 'ushort')

    if (!(vmode[1] %in% MODES)){
        stop(sprintf('Error: Incorrect argument "vmode". Allowable modes are %s', paste(MODES, collapse = ', ')))
    }

    if (length(gr) == 0 | !any(width(gr) > 0)){
        stop('Error: Trying to create ffTrack with empty GRanges.')
    }

    ## start and end indices of range in ff vector
    gr$ix.s = c(1, 1 + cumsum(as.numeric(width(gr)[-length(gr)])))

    .Object@.blocksize = .Machine$integer.max
    .Object@.gr = gr.stripstrand(gr[, 'ix.s'])
    .Object@.vmode = vmode[1]

    if (!grepl('\\.rds$',  file.name) & !grepl('\\.RDS$',  file.name)){
        file.name = paste(file.name, '.rds', sep = '')
    }

    ff.filename = gsub('rds$', 'ffdata', file.name)

    ## touch the files

    if ((file.exists(file.name) | file.exists(ff.filename)) & !overwrite){
        stop('Error: Target files already exist, to overwrite use overwrite = TRUE')
    }

    writeLines('', file.name)
    writeLines('', ff.filename)

    .Object@.rds.filename = normalizePath(file.name)
    .Object@.ff.filename = normalizePath(ff.filename)

    len = sum(as.numeric(width(gr)))

    ## allocate
    .Object@.ff = ff(default.val, length = pmin(len, .Object@.blocksize), vmode = .Object@.vmode, filename = .Object@.ff.filename, overwrite = overwrite)
    .Object@.ffaux = list()
    .Object@.length = len

    if (len > .Object@.blocksize){
        newlen = len - .Object@.blocksize;
        i = 1
        while (newlen>0){
            .Object@.ffaux[[i]] = ff(default.val, length = pmin(newlen, .Object@.blocksize), vmode = .Object@.vmode, filename = paste(.Object@.ff.filename, '.', i, sep = ''), overwrite = overwrite)
            newlen = newlen - .Object@.blocksize;
            i = i + 1
        }
    }

    ## don't allow factors for character vmode (remember vmode refers to the actual stored value, characters can be represented
    ## as factors by providing a 'character' level
    if (!(vmode %in% c('character')) & !is.null(levels)){
        .Object@.levels = levels
    }
    else{
        .Object@.levels = NA
    }

    saveRDS(.Object, .Object@.rds.filename)

    validObject(.Object)

    if (verbose){

        file.size = file.info(.Object@.ff.filename)$size / 1e6

        if (length(.Object@.ffaux) > 0){
            file.size = file.size + file.info(sapply(.Object@.ffaux, file.name))$size/1e6
        }

        cat('Created ffTrack object with .rds file %s and .ffdata base file %s spanning %s block(s) occupying %sM of disk\n',
            .Object@.rds.filename, .Object@.ff.filename, file.size)
    }

    print(.Object)

    return(.Object)

})


setValidity('ffTrack', function(object){

    problems = c();

    if (!file.exists(object@.ff.filename)){
        problems = c('ffdata file is missing')
    }
    else if (normalizePath(object@.ff.filename) != normalizePath(ff::filename(object@.ff))){
        problems = c('ffTrack filename does not match .ffdata filename')
    }

    if (length(object@.ffaux)>0){

        for (i in 1:length(object@.ffaux)){

            if (!file.exists(ff::filename(object@.ffaux[[i]]))){
                problems = c(problems, 'ffaux file is missing')
            }
            else if (normalizePath(paste(object@.ff.filename, '.', i, sep = '')) != normalizePath(ff::filename(object@.ffaux[[i]]))){
                problems = c(problems, 'ffaux object filename not compatible with .ff.filename')
            }
        }
    }

    if (!file.exists(object@.rds.filename)){
        problems = c('rds file is missing')
    }

    if (!grepl('\\.ffdata$', object@.ff.filename)){
        problems = c('ffTrack ff filename does not end in .ffdata')
    }

    if (!grepl('\\.rds$', object@.rds.filename)){
        problems = c('ffTrack filename does not end in .rds')
    }

    if (is.null(object@.gr$ix.s)){
        problems = c('internal gRanges missing ix.s meta data field')
    }

    if (!is.null(object@.gr$ix.s)){
        if (any(is.na(object@.gr$ix.s))){
            problems = c('internal gRanges ix.s data field is corrupt')
        }
    }

    if (!(length(problems)==0)){
        paste(problems, collapse = '\n')
    }
})


setMethod('show', 'ffTrack', function(object){
    validObject(object)
    fn = object@.ff.filename

    if (length(object@.ffaux) > 0){
        fn = paste(c(fn, sapply(object@.ffaux, filename)), collapse = ', ')
    }
    
    cat(sprintf('ffTrack object of vmode %s of ffdata filename(s) %s comprising %sM of disk space and %s GRanges: \n', vmode(object), fn, round(size(object), 2), length(object@.gr)))
})




#' @name size
#' @title size
#' @description
#'
#' Determine size in MB for this object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('size', function(object) standardGeneric('size'))
setMethod('size', 'ffTrack', function(object)
{
    fn = object@.ff.filename

    sz = as.numeric(file.info(fn)$size / 1e6)

    if (length(object@.ffaux) > 0){
        fn = sapply(object@.ffaux, filename)
        sz = sz + sum(as.numeric(file.info(fn)$size)) / 1e6
    }

    return(sz)
})




#' @name vmode
#' @title vmode
#' @description
#'
#' vmode of ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('vmode', function(object) standardGeneric('vmode'))
setMethod('vmode', 'ffTrack', function(object){
    object@.vmode
})




#' @name len
#' @title len
#' @description
#'
#' len of ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('len', function(object) standardGeneric('len'))
setMethod('len', 'ffTrack', function(object){
    sum(as.numeric(width(object@.gr)))
})




#' @name levels
#' @title levels
#' @description
#'
#' get levels of ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('levels', function(object) standardGeneric('levels'))
setMethod('levels', 'ffTrack', function(object){
    object@.levels
})




#' @name set_levels
#' @title set_levels
#' @description
#'
#' set levels of ffTrack object
#'
#' @param value param info
#' @export
#' @author Marcin Imielinski
setGeneric('set_levels', function(object, value) standardGeneric('set_levels'))
setMethod('set_levels', 'ffTrack', function(object, value)
{
    if (!all(is.na(object@.levels))){
        stop('Error: Levels not defined in original ffTrack instantiation, please re-instantiate to add levels')
    }

    if (is.vector(value) & length(value) == length(object@.levels)){
        object@.levels = value
    }
    else{
        stop('Error: Replacement levels must be a vector of same length as current set of levels')
    }

    validObject(object)
    return(object)

})




#' @name ranges
#' @title ranges
#' @description
#'
#' ranges underlying ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('ranges', function(object) standardGeneric('ranges'))
setMethod('ranges', 'ffTrack', function(object){
    object@.gr
})




#' @name filename
#' @title filename
#' @description
#'
#' filename associated with ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('filename', function(object) standardGeneric('filename'))
setMethod('filename', 'ffTrack', function(object){
    c(ff = object@.ff.filename, rds = object@.rds.filename)
})





#' @name cp
#' @title copy ffTrack object to new path on the file system (and all data files)
#' @description
#'
#' copy ffTrack object to new path on the file system (and all data files)
#'
#' @export
#' @author Marcin Imielinski
setGeneric('cp', function(.Object, path, overwrite = FALSE) standardGeneric('cp'))
setMethod('cp', 'ffTrack', function(.Object, path, overwrite = FALSE){
    return(mv(.Object, path, overwrite = overwrite, keep.original = TRUE))
})




#' @name seqlengths
#' @title seqlengths
#' @description
#'
#' seqlengths of ffTrack object
#'
#' @export
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevels
#' @author Marcin Imielinski
setGeneric('seqlengths', function(object) standardGeneric('seqlengths'))
setMethod('seqlengths', 'ffTrack', function(object){
    seqlengths(object@.gr)
})




#' @name seqinfo
#' @title seqinfo
#' @description
#'
#' seqinfo of ffTrack object
#'
#' @export
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevels
#' @author Marcin Imielinski
setGeneric('seqinfo', function(object) standardGeneric('seqinfo'))
setMethod('seqinfo', 'ffTrack', function(object){
    seqinfo(object@.gr)
})




#' @name seqlevels
#' @title seqlevels of ffTrack object
#' @description
#'
#' seqlevels of ffTrack object
#'
#' @importFrom GenomeInfoDb seqlengths seqinfo seqlevels
#' @export
#' @author Marcin Imielinski
setGeneric('seqlevels', function(object) standardGeneric('seqlevels'))
setMethod('seqlevels', 'ffTrack', function(object){
    seqlevels(object@.gr)
})



#' @name del
#' @title delete ffTrack object (and all data files)
#' @description
#'
#' delete ffTrack object (and all data files)
#'
#' @export
#' @author Marcin Imielinski
setGeneric('del', function(.Object, path, overwrite = FALSE) standardGeneric('del'))
setMethod('del', 'ffTrack', function(.Object, path, overwrite = FALSE){

    fdel = c(.Object@.ff.filename, .Object@.rds.filename)

    if (length(.Object@.ffaux) > 0){
        fdel = c(fdel, sapply(.Object@.ffaux, filename))
    }

    if (any(file.exists(fdel))){
        i = sapply(fdel, function(x) system(sprintf('rm -f %s', x)))
        cat(sprintf('Removed files: %s\n', paste(fdel, collapse = ', ')))
    }
    else{
        cat('Object files already deleted')
    }
})



#' @name writeable_status
#' @title toggle writeable status of ffTrack object
#' @description 
#'
#' Toggle writeable status of ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('writeable_status', function(.Object, value) standardGeneric('writeable_status'))
setMethod('writeable_status', 'ffTrack', function(.Object, value)
{
    if (!is.logical(value)){
        stop('Error: value must be logical')
    }

    if (value){
        system(paste('chmod +w ', .Object@.ff.filename))
        close.ff(.Object@.ff)
        open.ff(.Object@.ff, readOnly = FALSE)
    }
    else{
        close.ff(.Object@.ff)
        open.ff(.Object@.ff, readOnly = TRUE)
    }

    return(.Object)

})




#' @name [
#' @title [
#' @description
#'
#' Takes as input either GRanges or GRangesList "i", and returns a vector or vector list (respectively) of
#' data from the corresponding ranges.
#'
#' Strand is taken into account here - i.e. a negative range will yield reversed data (note: not reverse complement)
#'
#' @param x param info
#' @param i GRanges or GRangesList
#' @param gr logical flag whether to return GRanges with field $Score populated with values (=FALSE); if T will return GRanges with field $score populated w values
#' @param raw logical flag whether to convert raw data to levels (if levels exist) (=FALSE); if T will not convert raw data to levels (if levels exist)
#' @return either vector (if i is a GRanges) or vector list (if i is a GRangesList)
#' @export
#' @author Marcin Imielinski
## setGeneric('[', function(x, i, gr, raw) standardGeneric('['))
setMethod('[', 'ffTrack', function(x, i, gr = FALSE, raw = FALSE)
{
    if (inherits(i, 'GRangesList')){
        query = grl.unlist(i)
    }
    else if (inherits(i, 'GRanges')){
        query = i
    }
    else{
        stop('Error: ffTrack accessor index must be a GRanges or GRangesList')
    }

    out = rep(NA, sum(as.numeric(width(query))))

    ov = gr.findoverlaps(query, x@.gr)

    if (length(ov) > 0){

        q.ix.s = c(1, 1 + cumsum(width(query))[-length(query)])
        q.ix1 = start(ov) - start(query)[ov$query.id] + q.ix.s[ov$query.id]
        q.ix2 = end(ov) - start(query)[ov$query.id] + q.ix.s[ov$query.id]

        s.ix1 = start(ov) - start(x@.gr)[ov$subject.id] + x@.gr$ix.s[ov$subject.id]
        s.ix2 = end(ov) - start(x@.gr)[ov$subject.id] + x@.gr$ix.s[ov$subject.id]

        q.ix = do.call('c', lapply(1:length(q.ix1), function(j) q.ix1[j]:q.ix2[j]))
        s.ix = do.call('c', lapply(1:length(s.ix1), function(j) s.ix1[j]:s.ix2[j]))

        aux.ix = s.ix > x@.blocksize

        out[q.ix[!aux.ix]] = x@.ff[s.ix[!aux.ix]]

        if (any(aux.ix)){
            aux.chunk = floor(s.ix[aux.ix] / x@.blocksize)

            for (j in unique(aux.chunk)){
                tmp.ix = which(aux.chunk == j)
                out[q.ix[aux.ix][tmp.ix]] = x@.ffaux[[j]][s.ix[aux.ix][tmp.ix]-j*x@.blocksize]
            }
        }

        if (!raw & !all(is.na(x@.levels))){
            out = as.integer(out)
            out[out==0] = NA ## 0 = NA for types where NA is cast to 0 (e.g. ubyte)
            out = x@.levels[out]
        }

        ## reverse data for negative strand queries
        if (any(ix = as.logical(strand(query)=='-'))){
            w = width(query)
            q.id = unlist(lapply(1:length(query), function(x) rep(x, w[x])))
            q.l = split(1:length(out), q.id)

            for (j in q.l[ix]){
                out[j] = rev(out[j])
            }
        }

        if (gr){
            tmp.out = gr.dice(i[, c()])
            tmp.out$score = out
            tmp.out = tmp.out[!is.na(tmp.out$score)]
            out = tmp.out
        }
    }

    if (inherits(i, 'GRangesList')){
        out = split(out, rep(query$grl.ix, width(query)))
    }

    return(out)

})




#' @name writeable
#' @title writeable
#' @description
#'
#' Access writeable status of ffTrack object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('writeable', function(.Object) standardGeneric('writeable'))
setMethod('writeable', 'ffTrack', function(.Object) file.access(.Object@.ff.filename, 2)==0 & !is.readonly(.Object@.ff))




#' @name mv
#' @title  mv
#' @description
#'
#' moves location of .Object to new filepath, returns new updated object
#'
#' @export
#' @author Marcin Imielinski
setGeneric('mv', function(.Object, path, overwrite = FALSE, keep.original = FALSE) standardGeneric('mv'))
setMethod('mv', 'ffTrack', function(.Object, path, overwrite = FALSE, keep.original = FALSE)
{
    new.root.path = path
    path = gsub('(^|(.*\\/))?([^\\/]*)$', '\\2', path)
    .file.name <- function(paths) return(gsub('(^|(.*\\/))?([^\\/]*)', '\\3', paths)) ## put here so don't rely on skitools
    if (!file.exists(.Object@.ff.filename)){
        if (all(file.exists(paste(new.root.path, .file.name(sapply(c(list(.Object@.ff), .Object@.ffaux), filename)), sep = '/')))){
             warning('Warning: not finding source .ff filenames, but finding .ffdata files that the ffTrack in this object used to point to.  
                      Will try rebuilding, by linking the GRanges and vmode info in this object with these .ff files')
        }
    }

    if (!grepl('\\/', path)){
        path = paste('./', path, sep = '')
    }

    if (file.exists(path)){
        if (file.info(path)$isdir){
            path = paste(path, .file.name(.Object@.rds.filename), sep = '/')
        }
    }
    else if (grepl("\\/$", path)){
        stop('Error: directory does not exist, please create before moving / copying')
    }

    if (!grepl('\\.rds$', path)){
        path = paste(path, '.rds', sep = '')
    }

    ff.path = gsub('\\.rds', '\\.ffdata', path)

    if ((file.exists(ff.path) | file.exists(path)) & !overwrite){
        stop('Error: One or more of the target paths exist, rerun with overwrite = FALSE to overwrite')
    }

    ff.aux.path = c()
    if (length(.Object@.ffaux)){
        for (i in 1:length(.Object@.ffaux)){
            ff.aux.path[i] = paste(ff.path, '.', i, sep = '')
            if (file.exists(ff.aux.path[i]) & !overwrite){
                stop('Error: One or more of the target paths exist, rerun with overwrite = FALSE to overwrite')
            }
        }
    }

    fstring = 'cp %s %s'

    newobj = .Object
    newobj@.rds.filename = path
    newobj@.ff.filename = ff.path
    ff::filename(newobj@.ff) = newobj@.ff.filename

    if (keep.original){
        system(sprintf(fstring, ff.path, .Object@.ff.filename))  ## reverse copy since ff already "moves" for us
    }

    if (length(.Object@.ffaux)>0){
        for (i in 1:length(.Object@.ffaux)){
            ff::filename(newobj@.ffaux[[i]]) = ff.aux.path[i]
            if (keep.original){
                system(sprintf(fstring, ff.aux.path[i], ff::filename(.Object@.ffaux[[i]]))) ## reverse copy since ff already "moves" for us
            }
        }
    }

    saveRDS(newobj, path)
    validObject(newobj)

    return(newobj)

})




#' @name [<-
#' @title [<-
#' @description
#'
#' ffTrack data populator
#'
#' Takes as input only a GRanges object "i" and vector or list "value" of length(i) which is interpreted as follows
#' (1) if vector, then each value[j] corresponds is assigned to (entire) range i[j]
#' (if if list of vectors, then length(value[[j]]) must be equal to length(i[j])
#'
#' vmode of values must be also compatible or coercible to vmode(x)
#'
#' Will throw a warning if ranges in "i" are out of the "universe" of the ffTrack object.
#'
#' If raw = T and fft has .levels, then will populate entries directly (without factorizing)
#'
#' if vmode is numeric op can equal "+", "-", "*", "/" .. and results in an update of the current entries in the
#' ffTrack with op and value eg
#' ff[i, op="+"] = value is same as ff[i, full = TRUE] = ff[i] + rep(value, width(i))
#' ff[i, op="-"] = value is same as ff[i, full = TRUE] = ff[i] + rep(value, width(i))
#' and so on.
#'
#' @param i GRanges specifying intervals of ffTrack to populate
#' @param value if i is GRanges then value is length(i) vector of data values, if i is GRangesList then value[[i]] is length (width(i)), but if full = TRUE, then value is a length(sum(width(i))) vector of data values, or if i is a GRangesList, then values[[j]] is a sum(width(i[[j]])) vector.
#' @param op operation can be either "+", "-", "*", "/", and only for numeric ffTrack
#' @param raw output raw integer data for factor as opposed to character
#' @param full logical flag whether the replacement value is a single vector whose sum is the summed width of the ranges "i"
#' @export
#' @author Marcin Imielinski
## setGeneric('[<-', function(x, i, value, op, raw, full) standardGeneric('[<-'))
setMethod('[<-', 'ffTrack', function(x, i, value, op = NULL, raw = TRUE, full = FALSE)
{
    query = i;

    if (!is.null(op)){
      
        ALLOWABLE.OPS = c('+', '-', '*', '/')
                    
        if (!(op %in% ALLOWABLE.OPS)){
            stop(sprintf('Error: argument "op" must be one of the following: %s', paste(ALLOWABLE.OPS, collapse = ',')))
        }

        if (vmode(x) %in% c('character') | !is.na(levels(x))){
            stop('Error: argument "op" can only be specified for numeric ffTrack, this track is either character or factor track')
        }
                      
    }

    if (!writeable(x)){
        stop('Error: oject is read-only, please make writeablb by setting argument "writeable" to TRUE')
    }

    if (!is(query, 'GRanges')){
        stop('Error: ffTrack index must be a GRanges')
    }

    if (full){
        if (length(value) != sum(width(i))){
            stop('Error: if full = TRUE then value must be of length = sum(width(ranges))')
        }
    }
    else if (is.list(value)){
        if (any(width(i) != sapply(value, length))){
           stop('Error: Mismatch between widths of input GRanges and value list')
        }
    }
    else{
        if (length(value) != length(i) & length(value) != 1){
            stop('Error: value must be list or vector of same length as GRanges input "i", or if full = TRUE a vector of same length as sum(width(granges))')
        }
    }

    ov = gr.findoverlaps(query, x@.gr)

    if (any(ix = (start(ov) != start(query)[ov$query.id] | end(ov) != end(query)[ov$query.id]))){
        warning(sprintf('Warning: Parts of %s ranges ignored', sum(ix)))
    }

    if (is.list(value) | full){
        values = unlist(value)
    }
    else{
        values = rep(value, width(i))
    }

    if (length(ov)>0){
        q.ix.s = c(1, 1 + cumsum(width(query))[-length(query)])
        q.ix1 = start(ov) - start(query)[ov$query.id] + q.ix.s[ov$query.id]
        q.ix2 = end(ov) - start(query)[ov$query.id] + q.ix.s[ov$query.id]

        s.ix1 = start(ov) - start(x@.gr)[ov$subject.id] + x@.gr$ix.s[ov$subject.id]
        s.ix2 = end(ov) - start(x@.gr)[ov$subject.id] + x@.gr$ix.s[ov$subject.id]

        q.ix = do.call('c', lapply(1:length(q.ix1), function(j) q.ix1[j]:q.ix2[j]))
        s.ix = do.call('c', lapply(1:length(s.ix1), function(j) s.ix1[j]:s.ix2[j]))

        aux.ix = s.ix > x@.blocksize

        ## reverse values for negative strand queries
        if (any(ix = as.logical(strand(query)=='-'))){
            w = width(query)
            q.id = unlist(lapply(1:length(query), function(x) rep(x, w[x])))
            q.l = split(1:length(values), q.id)
            for (j in q.l[ix]){
                values[j] = rev(values[j])
            }
        }
                 
        ## populate as factor if levels exist and raw = FALSE
        if (!(all(is.na(x@.levels))) & !raw){
            x@.ff[s.ix[!aux.ix]] = factor(values[q.ix[!aux.ix]], x@.levels)
        } 
        else{
            if (is.null(op)){
              x@.ff[s.ix[!aux.ix]] = values[q.ix[!aux.ix]]
            }
            else{
                if (op == '+'){
                    x@.ff[s.ix[!aux.ix]] = x@.ff[s.ix[!aux.ix]] + values[q.ix[!aux.ix]]
                }
                else if (op == '-'){
                    x@.ff[s.ix[!aux.ix]] = x@.ff[s.ix[!aux.ix]] - values[q.ix[!aux.ix]]
                }
                else if (op == '*'){
                    x@.ff[s.ix[!aux.ix]] = x@.ff[s.ix[!aux.ix]] * values[q.ix[!aux.ix]]
                }
                else if (op == '/'){
                    x@.ff[s.ix[!aux.ix]] = x@.ff[s.ix[!aux.ix]] / values[q.ix[!aux.ix]]
                }
            }
        }

        if (any(aux.ix)){
            aux.chunk = floor(s.ix[aux.ix] / x@.blocksize)

            for (j in unique(aux.chunk)){
                tmp.ix = which(aux.chunk == j)
                
                ## populate as factor if levels exist
                if (!all(is.na(x@.levels)) & !raw){
                    x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] = factor(as.vector(values[q.ix[aux.ix][tmp.ix]]), x@.levels)
                }
                else{
                    if (is.null(op)){
                        x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] = as.vector(values[q.ix[aux.ix][tmp.ix]])
                    }
                    else{
                        if (op == '+'){
                            x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] = x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] + as.vector(values[q.ix[aux.ix][tmp.ix]])
                        }
                        else if (op == '-'){
                            x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] = x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] - as.vector(values[q.ix[aux.ix][tmp.ix]])
                        }
                        else if (op == '*'){
                            x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] = x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] * as.vector(values[q.ix[aux.ix][tmp.ix]])
                        }
                        else if (op == '/'){
                            x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] = x@.ffaux[[j]][s.ix[aux.ix][tmp.ix] - j*x@.blocksize] / as.vector(values[q.ix[aux.ix][tmp.ix]])
                        }
                    }
                }
            }
        }
    }

    return(x)
          
})





## Constructor Functions

#' @name ffTrack
#' @title ffTrack
#' @description
#'
#' Constructs new \code{ffTrack} for offline storage of genomic data.  Allocates memory to store one of several data "modes"
#' (e.g. \code{numeric}, \code{byte}, \code{character}) data types across a fixed interval set (\code{GRanges}).
#' Useful for numeric (e.g. conservation track) or character (e.g. human genome sequence) data. 
#' Physical instantiation will result in the creation of one or more "heavy" *.ffData files and a lightweight *.rds pointer 
#' which is the ffTrack object that is returned by this function.  That object can be read or written to using GRanges indices.
#'
#' Initialization requires (1) a filename (2) a set of GRanges corresponding to the "space"
#' (3) a vmode (one of the following:
#' \itemize{
#' \item boolean (1 bit logical)
#' \item logical (2 bit logical + NA)
#' \item quad (2 bit unsigned integer without NA)   (NOTE: quad' allows efficient storage of genomic data as an 'A','T','G','C' factor. See 'ff' documentation.)
#' \item nibble (4 bit unsigned integer without NA)
#' \item byte (8 bit signed integer with NA)
#' \item ubyte (8 bit unsigned integer without NA)
#' \item short (16 bit signed integer with NA)
#' \item ushort (16 bit unsigned integer without NA)
#' \item integer (32 bit signed integer with NA)
#' \item single (32 bit float)
#' \item double (64 bit float)
#' \item raw (8 bit unsigned char)
#' }
#' Initialization will create two files (1) an .rds object meta data (2) *.ffdata binary ff object
#' These files should have static paths (i.e. should not be moved outside of R) - otherwise will break.
#' However the object will still be functional if the .rds file is moved to another location and loaded
#' from there.
#'
#' optional argument "levels" will (by default) convert stored values to levels prior to accessing, and populated
#' values to integers prior to storing.  In the file, levels will be indexed with 0-based indices (i.e. 0 will refer to the
#' the first level item)
#'
#' @param gr \code{GRanges} of input ranges
#' @param file.name Filename to store the genomic data
#' @param default.val default val
#' @param overwrite Whether to overwrite existing in the filename
#' @param levels Optional argument to specify unique levels for storage of factors,
#' @param verbose Set verbosity
#' @param vmode vmode
#' @return ffTrack object
#' @export
#' @author Marcin Imielinski
ffTrack = function(gr, file.name, default.val = NA, overwrite = FALSE, levels = NULL, verbose = FALSE, vmode = 'double', ...){
    new('ffTrack', gr = gr, file.name = file.name, default.val = default.val, overwrite = overwrite, levels = levels, verbose = verbose, vmode = vmode, ...)
}




#' @name bw2fft
#' @title bw2fft
#' @description
#'
#' Creates ffTrack object from input bigwig file.
#'
#' @import rtracklayer
#' @param bwpath path to BigWig
#' @param fftpath path to ffTrack .rds that will be created by this (by default .bw is replaced by .rds)
#' @param region GRanges specifying regions to limit ffTrack computation to (default is whole genome, ie seqnames of BigWig file); whether to limit to certain region (instead of whole genome)
#' @param chrsub whether to sub out 'chr' in seqnames / seqlevels of BigWig object; whether to sub in / sub out 'chr' when accessing bigwig file
#' @param verbose logical flag
#' @param buffer integer size of how big of a buffer to use when transferring data from BigWig to ffTrack object; number of bases to access at a time
#' @param skip.sweep logical flag (default FALSE) if TRUE will skip the sweep of "region" for the portions that have non-NA values; if TRUE will not sweep for covered region, just make a whole genome file or a file across provided regions
#' @param vmode character specifying vmode to use for encoding (default == 'double')
#' @param resume logical flag specifying whether to resume the populatino of an already existing ffTrack object (default FALSE)
#' @param min.gapwidth  minimum gap-width with which to merge reference adjacent intervals, this will mildly increase the file size but reduce the range complexity of the GRanges object; flank (to reduce the range complexity of the ffdata skeleton, but increase file size)
#' @return ffTrack object corresponding to the data in the BigWig file
#' @export
#' @author Marcin Imielinski
bw2fft = function(bwpath, fftpath = gsub('(\\.bw.*)|(\\.bigwig.*)', '.rds', bwpath), region = NULL, chrsub = TRUE, 
    verbose = FALSE, buffer = 1e5, skip.sweep = FALSE, vmode = 'double', resume = FALSE,  min.gapwidth = 1e3)
{

    mc.cores = 1;

    if (!is.null(region)){
        if (chrsub){
            region = gr.fix(gr.chr(region), seqinfo(BigWigFile(normalizePath(bwpath))), drop = TRUE)
        }
    }

    if (!skip.sweep){
        if (!is.null(region)){
            tiles = gr.tile(region, buffer)
        }
        else{
            tiles = gr.tile(si2gr(seqinfo(BigWigFile(normalizePath(bwpath)))), buffer)
        }

        if (verbose){
            cat(sprintf('\nInput path %s, Output path %s, Buffer %s, min.gapwidth %s', bwpath, fftpath, buffer, min.gapwidth))
        }

        if (verbose){
            cat(sprintf('\nSweeping BigWig file for covered positions across %s tiles covering %s bases with buffer size %s \n', length(tiles), sum(as.numeric(width(tiles))), buffer))
        }

        ## first sweep file to find all "covered" ranges (in order to make skeleton ffTrack object)
        covered = reduce(do.call('c', mclapply(1:length(tiles), function(x){
            if (verbose){
              cat(x, ' ')
            }
            reduce(import.ucsc(bwpath, selection = tiles[x], chrsub = FALSE), min.gapwidth = min.gapwidth)
        }, mc.cores = mc.cores)), min.gapwidth = min.gapwidth)
    }
    else if (!is.null(region)){
        covered = region
    }

    ## assume entire seqinfo is covered
    else{
        covered = si2gr(seqinfo(BigWigFile(normalizePath(bwpath))))
    } 
    
    if (chrsub){
        covered = gr.sub(covered, 'chr', '')
    }

    if (resume){
        fft = readRDS(fftpath)
    }
    else{
        fft = ffTrack(covered, fftpath)
    }

    if (verbose){
        cat(sprintf('\t.ffdata file %s has size %sM\n', filename(fft)['ff'], round(file.info(filename(fft)['ff'])$size/1e6, 2)))
    }

    covered.tile = gr.tile(covered, buffer)

    if (verbose){
        cat(sprintf('\nPopulating ffTrack at %s tiles covering %s bases with buffer size %s \n', length(covered.tile), sum(as.numeric(width(covered.tile))), buffer))
    }

    mclapply(1:length(covered.tile), function(x){
        if (verbose){
            cat(x, ' ')
        }
        tmp = import.ucsc(bwpath, selection = covered.tile[x], chrsub = chrsub)
        fft[tmp] = tmp$score
        gc()
    }, mc.cores = mc.cores)

    if (verbose){
        cat('\n')
    }

    return(fft)

}




#' @name wig2fft 
#' @title wig2fft 
#' @description
#'
#' Creates ffTrack object from input .wig file
#'
#' @import rtracklayer
#' @param wigpath path to Wig
#' @param fftpath path to ffTrack .rds that will be created by this (by default .bw is replaced by .rds)
#' @param chrsub whether to sub out 'chr' in seqnames / seqlevels of Wig object
#' @param verbose logical flag
#' @param buffer integer size of how big of a buffer to use when transferring data from Wig to ffTrack object; number of bases to access at a time
#' @param skip.sweep logical flag (default FALSE) if TRUE will skip the sweep of "region" for the portions that have non-NA values; if TRUE will not sweep for covered region, just make a whole genome file or a file across provided regions
#' @param seqlengths info
#' @param vmode  character specifyhing vmode to use for encoding (by default double)
#' @param gz stuff
#' @param bz2 stuff
#' @param min.gapwidth  minimum gap-width with which to merge reference adjacent intervals, this will mildly increase the file size but reduce the range complexity of the GRanges object; flank (to reduce the range complexity of the ffdata skeleton, but increase file size)
#' @return ffTrack object corresponding to the data in the Wig file
#' @export
wig2fft = function(wigpath, fftpath = gsub('(\\.wig.*)', '.rds', wigpath), chrsub = TRUE, verbose = FALSE,
    buffer = 1e5, skip.sweep = FALSE, seqlengths = hg_seqlengths(), vmode = 'double', gz = grepl('\\.gz$', wigpath),
    bz2 = grepl('\\.bz2$', wigpath), min.gapwidth = 1e3 ){

    if (gz){
        grepstr = sprintf('gunzip -c %s | grep -nP "\\S+Step" ', wigpath)
        wcstr = sprintf('gunzip -c %s | wc -l', wigpath)
    }
    else if (bz2){
        grepstr = sprintf('bunzip2 -c %s | grep -nP "\\S+Step" ', wigpath)
        wcstr = sprintf('gunzip -c %s | wc -l', wigpath)
    }
    else{
        grepstr = sprintf('grep -nP "\\S+Step" %s', wigpath)
        wcstr = sprintf('wc -l %s', wigpath)
    }

    p = pipe(grepstr); steps = readLines(p); close(p)
    p = pipe(wcstr); nlines = as.numeric(strsplit(readLines(p), ' ')[[1]][1]); close(p)

    type = grepl('variableStep', steps)

    if (all(type)){
        type = 'var' # TODO: implement variable step
        stop('Error: only fixedstep WIG currently supported')
    }
    else if (all(!type)){
        type = 'fixed'
    }
    else{
        stop('Error: Input format not supported or WIG file corrupt: input WIGS must be either all variable or all fixed step')
    }

    ## obtain ranges and step sizes for data
    if (verbose){
        cat(sprintf('Parsing %s step wig file %s with %s lines and %s ranges\n', type, wigpath, nlines, length(steps)))
    }

    if (type == 'fixed'){
        tmp = strsplit(steps[1], '(\\:)|( )')[[1]]
        ncol = length(tmp)
        col.names = gsub('(\\w+)\\=.*', '\\1', tmp)
        col.names[1:2] = c('line', 'type')

        if (verbose){
            cat(sprintf('Converting wig declarations to matrix\n'))
        }

        if (ncol == 4){
            tmp = tryCatch(matrix(unlist(strsplit(steps, '(\\:)|( )')), ncol = 4, byrow = TRUE, dimnames = list(c(), col.names)), error = function(x) NULL)

            if (is.null(steps)){
                stop('Error: FixedStep WIG file corrupt: make sure that every declaration line in file has same format with 4, 5, or 6 columns (+/- step, span) according to UCSC website')
            }

            tab = data.frame(chr = gsub('chrom\\=', '', tmp[, 'chrom']),
                start = as.numeric(gsub('start\\=', '', tmp[, 'start'])),
                line = as.numeric(gsub('line\\=', '', tmp[, 'line'])),
                step = 1,
                span = 1, 
                stringsAsFactors = FALSE)
        }
        else if (ncol == 5){

            tmp = tryCatch(matrix(unlist(strsplit(steps, '(\\:)|( )')), ncol = 5, byrow = TRUE, dimnames = list(c(), col.names)), error = function(x) NULL)

            if (is.null(steps)){
                stop('Error: FixedStep WIG file corrupt: make sure that every declaration line has same format with 5 or 6 columns (+/- span) according to UCSC website')
            }

            tab = data.frame(chr = gsub('chrom\\=', '', tmp[, 'chrom']),
                start = as.numeric(gsub('start\\=', '', tmp[, 'start'])),
                line = as.numeric(gsub('line\\=', '', tmp[, 'line'])),
                step = as.numeric(gsub('step\\=', '', tmp[, 'step'])),
                span = 1, 
                stringsAsFactors = FALSE)
            }
        else if (ncol == 6){

            tmp = tryCatch(matrix(unlist(strsplit(steps, '(\\:)|( )')), ncol = 5, byrow = TRUE, dimnames = list(c(), col.names)), error = function(x) NULL)

            if (is.null(tmp)){
                stop('Error: FixedStep WIG file corrupt: make sure that every declaration line has same format with 5 or 6 columns (+/- span) according to UCSC website')
            }

            tab = data.frame(chr = gsub('chrom\\=', '', tmp[, 'chrom']),
                start = as.numeric(gsub('start\\=', '', tmp[, 'start'])),
                line = as.numeric(gsub('line\\=', '', tmp[, 'line'])),
                step = as.numeric(gsub('step\\=', '', tmp[, 'step'])),
                span = as.numeric(gsub('span\\=', '', tmp[, 'span'])),
                stringsAsFactors = FALSE)
        }
        else{
            stop('Error: WIG file corrupt')
        }

        tab$width = (diff(c(tab$line, nlines+1))-1) * tab$step ## step tells us what is the (maximum) footprint of each sub-interval corresponding to a line
        tab$end = tab$start + tab$width -1
        tab$line.start = tab$line+1  ## beginning and ends of line
        tab$line.end = tab$line.start + tab$width -1

        if (chrsub){
            tab$chrom = gsub('chr', '', tab$chr)
        }

        if (verbose){
            cat(sprintf('Creating ffData with vmode %s for %s ranges spanning %s bases of sequence\n', vmode, nrow(tab), sum(tab$width)))
        }

        fft = ffTrack(reduce(gr.fix(seg2gr(tab, seqlengths = NULL), seqlengths), min.gapwidth = min.gapwidth), fftpath, vmode = vmode)

        if (verbose){
            cat(sprintf('\t.ffdata file %s has size %sM\n', filename(fft)['ff'], round(size(fft))))
        }

        ## now populate fft
        con = file(wigpath, 'r')
        tmp = readLines(con, tab$line[1])
        scores = NULL;
        curbuf = 0
        last.dump = 0
        w = tab$width;
        st = tab$step
        sp = tab$span;

        if (verbose){
            numpoints = 100
            last.point = 0
            cat('\nPopulating .. \n\nProgress bar:\n')
            cat(paste(rep('*', numpoints), collapse = ''), '\n')
        }

        for (i in 1:nrow(tab)){
            tmp = readLines(con, w[i])
            if (vmode != 'character'){
                tmp = as.numeric(tmp)
            }

            if (st[i] > 1){
                tmp2 = rep(NA, length(tmp)*st[i])
                tmp.st = ((1:length(tmp))-1)*st[i] + 1
                tmp.ix = unlist(mapply(function(s, e) s:e, tmp.st, tmp.st + sp[i] - 1))
                tmp2[tmp.ix] = tmp
                tmp = tmp2
            }
            else if (sp[i]>1){
                tmp = rep(tmp, each = sp[i])
            }

            scores = c(scores, list(tmp))
            curbuf = curbuf + length(tmp)

            if (curbuf > buffer){

                if (verbose){
                    if (((i/nrow(tab))*numpoints - last.point) > 1){
                        cat('*')
                        last.point = last.point + 1
                    }
                }

                fft[seg2gr(tab[(last.dump+1):i, ], seqlengths = NULL)] = scores
                last.dump = i;

                curbuf = 0
                scores = NULL
            }

            tmp = readLines(con, 1)

        }

        if (curbuf > 0){
            fft[seg2gr(tab[(last.dump+1):i, ], seqlengths = NULL)] = scores
        }

        if (verbose){
            if (((i/nrow(tab))*numpoints - last.point)>1){
                cat('*')
            }
            cat('\n')
        }

        close(con)

        return(fft)

    }
}




#' @name get_seq
#' @title get_seq 
#' @description
#'
#' Retrieve genomic sequenes
#'
#' Wrapper around getSeq which does the "chr" and seqnames conversion if necessary
#' also handles GRangesList queries
#'
#' @param hg A BSgenome or and ffTrack object with levels = c('A','TRUE','G','C','N')
#' @param gr GRanges object to define the ranges
#' @param unlist logical whether to unlist the final output into a single DNAStringSet. Default TRUE
#' @param mc.cores Optional multicore call. Default 1
#' @param mc.chunks Optional define how to chunk the multicore call. Default mc.cores
#' @param as.data.table boolean
#' @param verbose Increase verbosity
#' @return DNAStringSet of sequences
#' @export
get_seq = function(hg, gr, unlist = TRUE, mc.cores = 1, mc.chunks = mc.cores,
                   as.data.table = FALSE, verbose = FALSE)
{
    if (inherits(gr, 'GRangesList')){
        grl = gr;
        old.names = names(grl);
        gr = unlist(grl);
        names(gr) = unlist(lapply(1:length(grl), function(x) rep(x, length(grl[[x]]))))
        seq = get_seq(hg, gr, mc.cores = mc.cores, mc.chunks = mc.chunks, verbose = verbose)
        cl = class(seq)
        out = split(seq, names(gr))
        out = out[order(as.numeric(names(out)))]
        if (unlist){
            out = do.call('c', lapply(out, function(x) do.call(cl, list(unlist(x)))))
        }
        names(out) = names(grl)

        return(out)
    }
    else{

        if (is(hg, 'ffTrack')){
            if (!all(sort(levels(hg)) == sort(c('A', 'T', 'G', 'C', 'N')))){
                cat("ffTrack not in correct format for get_seq, levels must contain only: 'A', 'T', 'G', 'C', 'N'\n")
            }
        }
        ## only sub in 'chr' if hg is a BSGenome
        else{

            if (!all(grepl('chr', as.character(seqnames(gr))))){
                gr = gr.chr(gr)
            }

            gr = gr.fix(gr, hg)
  
            if (mc.cores > 1){

                ix = suppressWarnings(split(1:length(gr), 1:mc.chunks))

                if (is(hg, 'ffTrack')){

                    mcout = mclapply(ix, function(x){
                        tmp = hg[gr[x]]

                        if (any(is.na(tmp))){
                            stop('Error: ffTrack corrupt: has NA values, cannot convert to DNAString')
                        }

                        if (!as.data.table) {
                            bst = Biostrings::DNAStringSet(sapply(split(tmp, as.vector(Rle(1:length(x), width(gr)[x]))), function(y) paste(y, collapse = '')))
                            names(bst) = names(gr)[x]
                        } 
                        else {
                            bst = data.table(seq=sapply(split(tmp, as.vector(Rle(1:length(x), width(gr)[x]))), function(y) paste(y, collapse='')))
                            bst[, names:=names(gr)[x]]
                        }   

                        if (any(strand(gr)[x]=='-')){

                            ix.tmp = as.logical(strand(gr)[x]=='-')
                            if (!as.data.table){
                                bst[ix.tmp] = Biostrings::complement(bst[ix.tmp])
                            }
                            else{
                                bst$seq[ix.tmp] = as.character(Biostrings::complement(DNAStringSet(bst$seq[ix.tmp])))
                            }
                        }

                        if (verbose){
                            cat('.')
                        }

                        return(bst)

                    }, mc.cores = mc.cores)

                    if (!as.data.table){

                        if (length(mcout)>1){
                            tmp = c(mcout[[1]], mcout[[2]])
                        }

                    out = do.call('c', mcout)[order(unlist(ix))]

                    }
                    else{
                        out = rbindlist(mcout)
                    }
                }
                else{

                    out = do.call('c', mclapply(ix, function(x){

                        if (verbose){
                            cat('.')
                        }

                        return(getSeq(hg, gr[x]))

                    }, mc.cores = mc.cores))[order(unlist(ix))]

                    if (verbose){
                        cat('\n')
                    }
                }
            }
            else{

                if (is(hg, 'ffTrack')){

                    tmp = hg[gr]

                    tmp[is.na(tmp)] = 'N'

                    if (any(is.na(tmp))){
                        stop('Error: ffTrack corrupt: has NA values, cannot convert to DNAString')
                    }

                    if (as.data.table) {
                        bst = data.table(seq=sapply(split(tmp, as.vector(Rle(1:length(gr), width(gr)))), function(x) paste(x, collapse='')))
                        bst[, names:=names(gr)]
                    } 
                    else {
                        bst = DNAStringSet(sapply(split(tmp, as.numeric(Rle(1:length(gr), width(gr)))), function(x) paste(x, collapse = '')))
                        names(bst) = names(gr)
                    }

                    if (any(as.character(strand(gr))=='-')){
                        ix = as.logical(strand(gr)=='-')
                    }

                    if (!as.data.table) {
                        bstc = as.character(bst)
                        bstc[ix] = as.character(Biostrings::complement(bst[ix]))
                        bst = Biostrings::DNAStringSet(bstc)  ## BIZARRE bug with line below
                        #bst[ix] = Biostrings::complement(bst[ix])
                    } 
                    else {
                        bst$seq[ix] = as.character(Biostrings::complement(DNAStringSet(bst$seq[ix])))
                    }

                    return(bst)
                }
                else{
                    out = getSeq(hg, gr)
                }
            }

            return(out)

        }
    }
}




#' @name seq2fft
#' @title seq2fft
#' @description
#'
#' Creates ffTrack object from BSGenome or FASTA (coming soon) file
#'
#' ## will either convert (1) raw sequence (2) k-nucleotide context centered around base or (3) motifs defined by some dictionary (anchored at first base) into leveled ffTrack (i.e. integer track with populated levels field)
#'
#' @import rtracklayer
#' @param seq BSGenome object, ffTrack object representing genomic sequence, or (not yet supported) FASTA file
#' @param nnuc how many nucleotides to left and right to enumerate
#' @param dict this should be a character vector or DNAStringSet, overrides nnuc arg if not null
#' @param chrsub whether to sub in / sub out 'chr' when accessing seq file
#' @param neg whether to analyze sequence data on negative strand (i.e. motifs will be analyzed in rev complement)
#' @param region GRanges specifying regions to limit ffTrack computation to (default is whole genome, ie seqnames of BigWig file)
#' @param mc.cores currently mc.cores can only be one (weird mclapply bug when running)
#' @param verbose logical flag
#' @param buffer integer size of how big of a buffer to use when transferring data from BigWig to ffTrack object; number of bases to access at a time
#' @param skip.sweep logical flag (default FALSE) if TRUE will skip the sweep of "region" for the portions that have non-NA values; if TRUE will not sweep for covered region, just make a whole genome file or a file across provided regions
#' @param vmode  character specifyhing vmode to use for encoding (by default double)
#' @param min.gapwidth  minimum gap-width with which to merge reference adjacent intervals, this will mildly increase the file size but reduce the range complexity of the GRanges object; flank (to reduce the range complexity of the ffdata skeleton, but increase file size)
#' @return ffTrack object corresponding to the data in the BigWig file
#' @export
#'
seq2fft = function(seq, nnuc = 0, dict = NULL, chrsub = TRUE, neg = FALSE, region = NULL, mc.cores = 1, verbose = FALSE,
    buffer = 1e5, skip.sweep = FALSE, vmode = 'ubyte', min.gapwidth = 1e3)
{
    if (!inherits(seq, 'BSgenome') & !is(seq, 'ffTrack')){
        stop('Error: Only BSGenome and ffTrack input for seq currently supported')
    }

    if (is.null(region)){

        region = si2gr(seq)

        if (chrsub){
            region = gr.sub(region, 'chr', '')
        }

        if (neg){
            strand(region) = '-'
        }
    }

    tiles = gr.tile(region, buffer)

    context = FALSE
    
    if (is.null(dict)){

        DNA_BASES = c('A', 'G', 'C', 'T', 'N')
        dict = Biostrings::DNAStringSet(Biostrings::mkAllStrings(DNA_BASES, nnuc*2 + 1))
        context = TRUE
    }
    else{
        if (!is(dict, 'character'))
          dict = Biostrings::DNAStringSet(dict)
    }


    if (is(fftpath, 'ffTrack')){
      
        if (verbose){
            cat(sprintf('Populating ffTrack with filename %s with %s MB of sequence\n', filename(fftpath)['rds'], round(sum(as.numeric(width(region)))/1e6, 2)))
        }

        fft = fftpath ## append to existing fftpath
    }
    else{

        if (verbose){
            cat(sprintf('Making ffTrack for genome %s spanning %s MB of sequence\n', attributes(seq)$seqs_pkgname, round(sum(as.numeric(width(region)))/1e6, 2)))
        }

        fft = ffTrack(region, fftpath, levels = as.character(dict), vmode = vmode)

    }

    if (verbose){
        cat(sprintf('\t.ffdata file %s has size %sM\n', filename(fft)['ff'], round(size(fft))))
    }

    print(tiles)

    mclapply(1:length(tiles), function(x){

        if (verbose){
            cat(x, '\n')
        }

        tmp.tile = tiles[x]

        ## get flanking sequence if dict was null (i.e. we are looking for k-nucleotide context)
        if (context){
            tmp.tile = suppressWarnings(tmp.tile + nnuc)
        }

        tmp = get_seq(seq, tmp.tile)[[1]];

        ix = match.bs(tmp, dict)

        ## get flanking sequence if dict was null (i.e. we are looking for k-nucleotide context)
        if (context) {

            flank.left = start(tiles)[x] - start(tmp.tile)
            flank.right = end(tmp.tile) - end(tiles)[x]

            ix = ix[1:(length(ix)-2*nnuc)]

            ## pad left
            if ((pad = nnuc - flank.left) > 0){
                ix = c(rep(NA, pad), ix)
            }

            ## pad left
            if ((pad = nnuc - flank.right) > 0){
                ix = c(ix, rep(NA, pad))
            }
        }

        fft[tiles[x], raw = TRUE] = list(ix)


    }, mc.cores = mc.cores)

    return(fft)

}




#' @name fftab
#' @title Tabulate data in an \code{ffTrack}
#' @description
#'
#' Tabulates data in ffTrack file across a set of intervals (GRanges) by counting the number of positions 
#' matching a given "signature" or applying FUN to aggregate data.  
#' Returns the input GRanges populated with one or more meta data columns of counts or averages.
#'
#' Similar to gr.val in gUtils
#'
#' ff can be an ffTrack but also an RleList from same genome as intervals.
#'
#' returns a GRanges with additional columns for metadata counts
#'
#' @param ff  ffTrack or RleList to pull data from
#' @param intervals intervals
#' @param signatures Signatures is a named list that specify what is to be tallied.  Each signature (ie list element)
#' consist of an arbitrary length character vector specifying strings to %in% (grep = FALSE)
#' or length 1 character vector to grepl (if grep = TRUE)
#' or a length 1 or 2 numeric vector specifying exact value or interval to match (for numeric data)
#'
#' Every list element of signature will become a metadata column in the output GRanges
#' specifying how many positions in the given interval match the given query
#' @param FUN function to aggregate with (default is sum)
#' @param grep logical flag (default FALSE), if TRUE will treat the strings in signature as inputs to grep (instead of exact matches if FALSE)
#' @param mc.cores how many cores (default 1)
#' @param chunksize chunk of FF to bring into memory (i.e. the width of interval), decrease if memory becomes an issue
#' @param verbose logical flag
#' @param na.rm logical flag whether to remove na during aggregation.
#' @importFrom data.table rbindlist data.table setkey :=
#' @importFrom gUtils gr.sub seg2gr gr.stripstrand si2gr rle.query gr.fix gr.chr gr.tile grl.unlist gr.findoverlaps gr.dice hg_seqlengths
#' @export
fftab = function(ff, intervals, signatures = NULL, FUN = sum, grep = FALSE, mc.cores = 1, chunksize = 1e6, verbose = TRUE, na.rm = TRUE)
{
    id = ix = NULL ## NOTE fix

    if (!is(ff, 'ffTrack') & !is(ff, 'RleList')){
        stop('Error: Input ff should be ffTrack or RleList\n')
    }

    if (length(intervals)==0){
        stop('Error: Must provide non empty interavl input as GRanges')
    }

    if (!is.null(signatures)){
        if (!is.list(signatures)){
            stop('Error: Signatures must be a named list of arbitrary length character or length 1 or 2 numeric vectors')
        }

        if (is.null(names(signatures))){
            names(signatures) = paste('sig', 1:length(signatures), sep = '')
        }

        check = sapply(signatures, function(x){
            if (is.numeric(x)){
                return(length(x)>=1 & length(x)<=2)
            }
            if (is.character(x)){
                if (grep){
                    return(length(x)==1)
                }
            }
            return(TRUE)
        })

        if (!all(check)){
            stop('Error: signatures input is malformed, should be either length 1 or 2 numeric, length 1 character (if grep = TRUE), or atbitrary length character otherwise)')
        }
    }
    else{
        signatures = list(score = numeric()) ## we are just scoring bases
    }

    ## generate command that will be executed at each access
    cmd = paste('list(', paste(sapply(names(signatures), function(x){

        sig = signatures[[x]]
        if (is.numeric(sig)){

            if (length(sig)==0){
                cmd = sprintf('%s = FUN(dat, na.rm = na.rm)', x)
            }
            else if (length(sig)==1){
                cmd = sprintf('%s = FUN(dat == %s, na.rm = na.rm)', x, sig[1])
            }
            else{
                cmd = sprintf('%s = FUN(dat > %s & dat< %s, na.rm = na.rm)', x, sig[1], sig[2])
            }

        }
        else{

            if (grep){
                cmd = sprintf('%s = FUN(grepl("%s", dat), na.rm = na.rm) ', x, sig[1])
            }
            else{
                cmd = paste(x, '= FUN(dat %in%', paste('c(', paste("\"", sig, "\"", sep = '', collapse = ','), '), na.rm = na.rm)', sep = ''))
            }
        }
    }), collapse = ', ', sep = ''), ')', sep = '')

    val = values(intervals)
    intervals$ix = 1:length(intervals)

    if (verbose){
        cat('Made command\n')
    }

    ## sorting will hopefully make data access more efficient
    gr = sort(intervals[, 'ix'])

    if (verbose){
        cat('Sorted intervals\n')
    }

    ## tailor the chunking to the size of the individual segments
    gr$chunk.id = ceiling(cumsum(as.numeric(width(gr)))/chunksize)
    gr$num = 1:length(gr)

    ## get down to business
    chunks = split(1:length(gr), gr$chunk.id)
    if (verbose){
        cat('Split intervals\n')
    }

    out = rbindlist(mclapply(chunks, function(ix){

        chunk = gr[ix]
        if (verbose){
            cat(sprintf('Intervals %s to %s of %s, total width %s, starting\n', chunk$num[1], chunk$num[length(chunk)], length(gr), sum(width(chunk)))) 
        }

        if (is(ff, 'ffTrack')){
            tmp = data.table(dat = ff[chunk], id = rep(1:length(chunk), width(chunk)) )
        }
        ## also can handle rle data
        else{
            tmp = data.table(dat = as.numeric(rle.query(ff, chunk)), id = rep(1:length(chunk), width(chunk)) )
        }

        setkey(tmp, id)

        if (verbose){
            cat(sprintf('Intervals %s to %s of %s: read in ff data\n', chunk$num[1], chunk$num[length(chunk)], length(gr)))
        }

        tab = tmp[, eval(parse(text=cmd)), keyby = id]

        if (verbose){
            cat(sprintf('Intervals %s to %s of %s: tabulated\n', chunk$num[1], chunk$num[length(chunk)], length(gr)))
        }

        ix = chunk$ix
        out = tab[1:length(chunk), ]
        out$id = NULL
        out$ix = ix
        if (verbose){
            cat(sprintf('Intervals %s to %s of %s: FINISHED\n', chunk$num[1], chunk$num[length(chunk)], length(gr)))
        }

        return(out)

    }, mc.cores = mc.cores))

    setkey(out, ix)
    out = as.data.frame(out[list(1:length(intervals)), ])
    out$ix = NULL
    values(intervals) = cbind(val, out)

    return(intervals)

}




#' @name match.bs
#' @title 
#' @description
#'
#' Identify matches between query and dictionary
#'
#' Wrapper around matchPdict to identify matches between a query
#' string query and dictionary dict (both BString objects or subclasses)
#'
#' @param query Query
#' @param dict Dictionary
#' @param midpoint boolean Flag for output the coordinates of the match as the location, where the midpoint of the dict string matches the given query. Default FALSE
#' @return a vector of indices of length width(query) that contains indices of the (starting) dictionary in the query string
#' @export
match.bs = function(query, dict, midpoint = FALSE)
{
    names(dict) = as.character(1:length(dict))

    tmp = sort(unlist(Biostrings::matchPDict(dict, query)))
    out = rep(NA, length(query))

    if (!midpoint){
      out[start(tmp)] = as.numeric(names(tmp))
    }
    else{
      out[floor((start(tmp)+end(tmp))/2)] = as.numeric(names(tmp))
    }

    return(out)
}





