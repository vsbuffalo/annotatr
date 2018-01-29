## annotatr.r -- tidy annotation helpers

extract_attr <- function(x, key)  {
  has_key <- regexpr(key, x)
  value <- gsub(sprintf('.*%s=([^;]+).*', key), '\\1', x)
  value[has_key == -1] <- NA
  out <- strsplit(value, ',')
  multi_value <- any(sapply(out, length) > 1)
  if (!multi_value) return(unlist(out))
  out
}

#' Read a GFF file
#' 
#' @param file the GFF file (can be gzipped) to read in.
#'
#' Read in a GFF file. This will not parse the attributes 
#' column, as lazyily extracting what we need from this is 
#' faster.
#' 
#' @export
read_gff <- function(file) {
  gff_cols <- readr::cols(chrom=readr::col_character(),
                   source=readr::col_character(),
                   feature=readr::col_character(),
                   start=readr::col_integer(),
                   end=readr::col_integer(),
                   score=readr::col_double(),
                   strand=readr::col_character(),
                   frame=readr::col_character(),
                   attr=readr::col_character())

  out <- readr::read_tsv(file, comment='#', 
                  col_names=names(gff_cols$cols),
                  col_types=gff_cols,
                  na='.')
  # Convert repetitive columns to factors. Specifying 
  # factor columns with readr can cause issues when reading in large files.
  dplyr::mutate(out, chrom=factor(chrom), source=factor(source), 
         feature=factor(feature))
}

#' Tidy a GFF tibble
#'
#' @param x a tibble of a GFF file.
#' @param chroms a character vector of chromosomes to include.
#' @param features a character vector of features
#'    of the file to include (others filtered out).
#'
#' @export
tidy_gff <- function(x, chroms=NULL,
                     features=c('five_prime_UTR', 'three_prime_UTR', 
                                'exon', 'intron', 'gene', 'CDS')) {
  out <- dplyr::select(dplyr::mutate(dplyr::filter(x, feature %in% features), 
                id=extract_attr(attr, 'ID'),
                parent=extract_attr(attr, 'Parent'),
                parent_type=extract_attr(attr, 'parent_type')),
        chrom, start, end, strand, feature, id, parent, parent_type)
  if (!is.null(chroms)) {
    out[out$chrom %in% chroms, ]
  }
  out
}

#' Create a GenomicRanges object from a tibble
#'
#' @param x a tibble with columns (chrom | chr) and (start,end | pos). 
#'   Optionally, a strand column. All other columns will be included
#'   as metadata.
#'
#' @export
to_granges <- function(x) {
 if (class(x) == 'GRanges') return(x)
 if (!('chrom' %in% colnames(x))) {
    if ('chr' %in% colnames(x)) {
      # rename column name
      colnames(x)[which(colnames(x) == 'chr')] <- 'chrom'
    } else {
      stop("'x' must contain either a 'chrom' or 'chr' column")
    }
  }
  if (!('start' %in% colnames(x) && 'end' %in% colnames(x))) {
    if ('pos' %in% colnames(x)) {
      x$start <- x$pos
      x$end <- x$pos
      x <- x[,-which(colnames(x) == 'pos')]
    } else {
      msg <- "'x' must have either 'start' and 'end', or 'pos' in column names"
      stop(msg)
    }
  }
  out <- GenomicRanges::GRanges(seqnames=as.character(x$chrom),
                                        IRanges::IRanges(x$start, x$end))
  if ('strand' %in% colnames(out))
    strand(out) <- x$strand
 
  range_cols <- c('chrom', 'start', 'end', 'strand')
  metacols <- x[, !(colnames(x) %in% range_cols)]
  GenomicRanges::mcols(out) <- S4Vectors::DataFrame(metacols)
  out
}

annot_table <- function(x, features) {
  tbl <- table(factor(x, levels=features))
  tibble::as_tibble(t(as.matrix(tbl)))
}

na_to_0 <- function(x) ifelse(is.na(x), 0, x)

vmessage <- function(msg, verbose) {
  if (verbose)
    message(msg)
}

#' Annotate the overlaps of dataset using a tibble or GRanges object of features
#'
#' @param x a tibble that can be converted to a GRanges object.
#' @param annot annotation tibble or GRanges object, with a \code{features} 
#'   of the feature labels.
annotate_overlaps <- function(x, annot, verbose=TRUE) {
  if (is.list(annot$parent)) {
    annot <- tidyr::unnest(annot, parent)
  }
  vmessage('converting to GRanges...', verbose)
  # convert to GRanges
  annot_gr <- to_granges(annot)
  x_gr <- to_granges(x)
  vmessage('finding overlaps...', verbose)
  hits <- findOverlaps(x_gr, annot_gr)

  # query the annotation for overlaps
  annot_df <- tibble::tibble(id=queryHits(hits), 
                             feature=annot$feature[subjectHits(hits)])

  # count up hits to features, and spread the columns
  annot_counts <- dplyr::group_by(dplyr::mutate(annot_df, count=1), feature, id)
  # if (reduce)
    # reduce_fun <- function(x) as.integer(x > 0)
  vmessage('summarizing and spreading feature columns...', verbose)
  annot_df <- tidyr::spread(dplyr::summarize(annot_counts, 
                                             count=sum(count)), 
                            feature, count, fill=0)

  # join the new columns to original matrix
  vmessage('joining to other original data...', verbose)
  jdf <- dplyr::left_join(tibble::tibble(id=seq_len(nrow(x))), annot_df)

  # drop id column and convert NA to 0
  jdf <- dplyr::mutate_all(jdf[, -1], funs(na_to_0)) 
  bind_cols(x, jdf)
}
