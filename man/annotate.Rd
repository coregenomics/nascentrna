% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classify.R
\name{annotate}
\alias{annotate}
\alias{da_tss}
\alias{ua_rna}
\alias{enhancer}
\title{Classify all gene positions.}
\usage{
annotate(tss, genes)

da_tss(tss, genes, min_start = 0, max_start = 300)

ua_rna(tss, genes, min_start = 0, max_start = 300)

enhancer(tre, genes)
}
\arguments{
\item{tss}{GRanges of experimentally found active gene region start sites to
be classified.}

\item{genes}{GRanges of known genes around which to classify \code{tss}
regions.}

\item{min_start}{Minimum integer distance away from gene start to find
feature.}

\item{max_start}{Maximum integer distance away from gene start to find
feature.}

\item{tre}{GRanges of experimentally found active putative transcription
regulatory regions to be classified.}
}
\value{
\code{\link{annotate}} returns GRanges of same size as \code{tss}
    with an additional \code{class} column of best position classification,
    and associated \code{gene} column of the corresponding gene.

\code{\link{da_tss}} returns a Hits object mapping gene with
    downstream-antisense TSS with a column for the tss \code{query} index
    (queryHits), genes \code{subject} index (subjectHits) and the
    \code{distance} between the pair.

\code{\link{ua_rna}} returns a Hits object mapping gene with
    downstream-antisense RNA with a column for the tss \code{query} index
    (queryHits), genes \code{subject} index (subjectHits) and the
    \code{distance} between the pair.

\code{\link{enhancer}} subsets \code{tre} to values outside genes.
}
\description{
\code{\link{annotate}} is a convenience wrapper for running all other
classification functions in this package using default values.  See the
package vignette for the decision tree used to determine the logical steps of
classification.

\code{da_tss} tests for downstream-antisense TSS.

\code{ua_rna} tests for upstream-antisense RNA.  The distance at which
neighboring promoters appear is ~300 bp according to Fig 2d in Chen et al
2016 \url{http://dx.doi.org/10.1038/ng.3616}

\code{enhancer} tests for transcription regulatory elements outside of
annotated genes.
}
