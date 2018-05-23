#' @importFrom BiocGenerics invertStrand
#' @importFrom GenomicRanges distanceToNearest isDisjoint mcols mcols<-
#'     promoters psetdiff resize
#' @importFrom IRanges CharacterList findOverlaps
#' @importFrom S4Vectors from
NULL

#' Classify all gene positions.
#'
#' \code{\link{annotate}} is a convenience wrapper for running all other
#' classification functions in this package using default values.  See the
#' package vignette for the decision tree used to determine the logical steps of
#' classification.
#'
#' @param tss GRanges of experimentally found active gene region start sites to
#'     be classified.
#' @param genes GRanges of known genes around which to classify \code{tss}
#'     regions.
#' @param min_start Minimum integer distance away from gene start to find
#'     feature.
#' @param max_start Maximum integer distance away from gene start to find
#'     feature.
#' @return annotate returns GRanges of same size as \code{tss} with an
#'     additional \code{class} column of best position classification, and
#'     associated \code{gene} column of the corresponding gene.
annotate <- function(tss, genes) {
    hits <- da_tss(tss, genes)
    mcols(tss)$class <- NA
    mcols(tss)$class[from(hits)] <- "daTSS"
    tss
}

#' \code{da_tss} tests for downstream-antisense TSS.
#'
#' @rdname annotate
#' @return \code{\link{da_tss}} returns a Hits object mapping gene with
#'     downstream-antisense TSS with a column for the tss \code{query} index
#'     (queryHits), genes \code{subject} index (subjectHits) and the
#'     \code{distance} between the pair.
da_tss <- function(tss, genes, min_start = 0, max_start = 300) {
    ## Input checks.
    if (! isDisjoint(genes)) {
        stop(sQuote("genes"), " has overlapping ranges.",
             " This will likely create invalid results.",
             " Please ensure ", sQuote("genes"), " is disjoint.")
    }
    ## More useful error message than downstream complaints of invalid ranges.
    if (min_start >= max_start) {
        stop(sQuote("min_start"), " must be less than ", sQuote("max_start"))
    }

    look_in <- promoters(genes, downstream = max_start)
    look_in <- psetdiff(
        look_in,
        promoters(genes, downstream = min_start))
    ## The distance calculation methods don't have an option to "look on the
    ## opposite strand", therefore invert the strand as a workaround.
    tss_inv <- invertStrand(tss)
    ol_proximal <- findOverlaps(resize(tss_inv, 0),
                                look_in,
                                ## Overlap includes touching each other.
                                maxgap = 0)
    idx_proximal <- from(ol_proximal)
    dist <- distanceToNearest(tss_inv[idx_proximal],
                              resize(look_in, 0))
    ## Match indices of tss input.
    ##
    ## Using S4Vectors::remapHits() (# nolint) doesn't seem to address this use
    ## case.  Will discuss with upstream.
    dist@from <- idx_proximal[from(dist)]
    dist@nLnode <- length(tss)
    dist
}
