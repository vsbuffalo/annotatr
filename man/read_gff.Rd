% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotatr.r
\name{read_gff}
\alias{read_gff}
\title{Read a GFF file}
\usage{
read_gff(file)
}
\arguments{
\item{file}{the GFF file (can be gzipped) to read in.

Read in a GFF file. This will not parse the attributes 
column, as lazyily extracting what we need from this is 
faster.}
}
\description{
Read a GFF file
}
