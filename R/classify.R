#' @importFrom BiocGenerics invertStrand
#' @importFrom GenomicRanges distanceToNearest isDisjoint promoters psetdiff
#'     resize
#' @importFrom IRanges CharacterList findOverlaps
#' @importFrom S4Vectors from to Hits mcols mcols<-
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
#' @return \code{\link{annotate}} returns GRanges of same size as \code{tss}
#'     with an additional \code{class} column of best position classification,
#'     and associated \code{gene} column of the corresponding gene.
#' @export
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
#' @export
da_tss <- function(tss, genes, min_start = 0, max_start = 300) {
    ## Input checks.
    check_disjoint(genes)
    ## More useful error message than downstream complaints of invalid ranges.
    check_lt(min_start, max_start)

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
    Hits(from = idx_proximal[from(dist)],
         to = to(dist),
         nLnode = length(tss),
         nRnode = length(genes),
         distance = mcols(dist)$distance)
}

#' \code{ua_rna} tests for upstream-antisense RNA.  The distance at which
#' neighboring promoters appear is ~300 bp according to Fig 2d in Chen et al
#' 2016 \url{http://dx.doi.org/10.1038/ng.3616}
#'
#' @rdname annotate
#' @return \code{\link{ua_rna}} returns a Hits object mapping gene with
#'     downstream-antisense RNA with a column for the tss \code{query} index
#'     (queryHits), genes \code{subject} index (subjectHits) and the
#'     \code{distance} between the pair.
#' @export
ua_rna <- function(tss, genes, min_start = 0, max_start = 300) {
    ## Input checks.
    check_disjoint(genes)
    ## More useful error message than downstream complaints of invalid ranges.
    check_lt(min_start, max_start)

    look_in <- promoters(genes, upstream = max_start, downstream = -min_start)
    ## The distance calculation methods don't have an option to "look on the
    ## opposite strand", therefore invert the strand as a workaround.
    tss_inv <- invertStrand(tss)
    ol_proximal <- findOverlaps(resize(tss_inv, 0, fix = "end"),
                                look_in,
                                ## Overlap includes touching each other.
                                maxgap = 0)
    idx_proximal <- from(ol_proximal)
    dist <- distanceToNearest(tss_inv[idx_proximal],
                              resize(look_in, 0, fix = "end"))
    ## Match indices of tss input.
    Hits(from = idx_proximal[from(dist)],
         to = to(dist),
         nLnode = length(tss),
         nRnode = length(genes),
         distance = mcols(dist)$distance)
}

#' Object validation for nascentrna package.
#'
#' \code{\link{check_disjoint}} throws an error if the input has overlapping
#' regions.
#'
#' @rdname validation
#' @param gr GRanges object to validate.
#' @return \code{\link{check_disjoint}} invisibly returns TRUE if the input is
#'     disjoint.
check_disjoint <- function(gr) {
    if (! isDisjoint(gr)) {
        obj <- sQuote(substitute(gr))
        stop(obj, " has overlapping ranges.",
             " This will likely create invalid results.",
             " Please ensure ", obj, " is disjoint.")
    }
    invisible(TRUE)
}

#' \code{\link{check_lt}} throws an error if the \code{lesser} argument is equal
#' to or exceeds the \code{greater} argument.
#'
#' @rdname validation
#' @param lesser Numeric value that should be less than \code{greater}.
#' @param greater Numeric value used to validate \code{lesser}.
#' @return \code{\link{check_lt}} invisibly returns TRUE if \code{lesser} is
#'     strictly less than \code{greater}
check_lt <- function(lesser, greater) {
    if (lesser >= greater) {
        obj_lesser <- sQuote(substitute(lesser))
        obj_greater <- sQuote(substitute(greater))
        stop(obj_lesser, " must be less than ", obj_greater)
    }
    invisible(TRUE)
}
