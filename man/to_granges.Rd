% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotatr.r
\name{to_granges}
\alias{to_granges}
\title{Create a GenomicRanges object from a tibble}
\usage{
to_granges(x)
}
\arguments{
\item{x}{a tibble with columns (chrom | chr) and (start,end | pos). 
Optionally, a strand column. All other columns will be included
as metadata.}
}
\description{
Create a GenomicRanges object from a tibble
}
