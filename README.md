## anotator -- simple annotation for population genomic variant data

**Extreme beta, proceed with caution.** This is mostly code I've extracted from
a project because (1) I like working with tibble more than `GRanges` and (2) I
feel like what I (and maybe others) want to do 90% of the time is simple and
deserves a simple way to do it.

The purpose of this package is annotate variants by overlapping genomic
features using a subset of a GFF file. The end result is the exact same
dataframe/tibble you put in, plus some additional columns indicator how many of
each feature are overlapped by a particular variant row. It's much like
`bedtools annotate`, but meant to be party of tidyverse pipelines. It feels
like 90% of annotation problems come down to finding overlapping features and
appending them as columns, yet (AFAIK) there's no set of functions to get us
there.

For extremely thorough variant annotation (e.g.  finding LoF variants, finding
frameshifts, etc.) see GEMINI, SNPeff, or Annovar.

It uses `GenomicRanges` for the overlap finding, but is designed to work with
the tidyverse. This is sort of the next iteration of my experimental
[gplyr](https://github.com/vsbuffalo/gplyr), with a narrower focus on simple
but correct annotation of variants.

```r
# Read in the GFF file. This is lazy in parsing the attributes column
# and extracts only a few key-values later (e.g. ID, parent, parent_type).
> gff <- read_gff('inst/data/dmel-all-r5.30.gff.gz')
> tgff <- gff %>% tidy_gff()  ## only keeps specified features
> tgff
# A tibble: 222,629 x 8
    chrom start   end strand feature                                 id
   <fctr> <int> <int>  <chr>  <fctr>                              <chr>
 1     2L  7529  9484      +    gene                        FBgn0031208
 2     2L  7529  8116      +    exon                      FBgn0031208:1
 3     2L  7680  8116      +     CDS              CDS_FBgn0031208:1_738
 4     2L  8117  8192      +  intron intron_FBgn0031208:1_FBgn0031208:2
 5     2L  8117  8192      +  intron intron_FBgn0031208:1_FBgn0031208:3
 6     2L  8193  9484      +    exon                      FBgn0031208:3
 7     2L  8193  8610      +     CDS              CDS_FBgn0031208:3_738
 8     2L  8193  8589      +     CDS              CDS_FBgn0031208:2_738
 9     2L  8193  8589      +    exon                      FBgn0031208:2
10     2L  8590  8667      +  intron intron_FBgn0031208:2_FBgn0031208:4
# ... with 222,619 more rows, and 2 more variables: parent <list>,
#   parent_type <chr>
```

Note that `parent` is a list-column; you can use `tidyr`'s `unnest(parent)` to
expand out each parent gene.

Now, imagine we have a tibble of data we want to annotate. This should have
columns `chrom` and `pos`, or `chrom`, `start`, and `end`:

```r
# our dataset
> d
# A tibble: 1,503,403 x 10
   chrom   pos   ref sample  time   rep treatment     A     T     C     G     N
   <chr> <int> <chr>  <chr> <int> <int>     <chr> <int> <int> <int> <int> <int>
 1    2L  15012     G  cont_rep1     0     1   control     0     9     0    21     0
 2    2L  15413     C  cont_rep1     0     1   control     0     0     9    23     0
 [...] 

## annotate this dataframe:
> da <- annotate_overlaps(d, tgff)
> a %>% select(chrom, pos, CDS:intron)  # look at the new columns
# A tibble: 14,053,853 x 6
   chrom   pos   CDS  exon  gene intron
   <chr> <int> <dbl> <dbl> <dbl>  <dbl>
 1    2L  15012     2     3     2      0
 2    2L  15413     2     3     2      0
 3    2L  15418     1     1     2      2
 4    2L  15447     1     1     2      2
 [...]
```

Note because we have overlapping exon annotations, we can have a variant
classified as both an intron and exon. By default, the package will **not**
collapse exons (e.g. by using `GenomicRanges` `reduce()`). However, you can use
`dplyr` to do this yourself.

```r
> da_collapsed_exons <- da %>% select(chrom, pos, CDS:intron) %>% 
       mutate(intron=ifelse(gene > 1 & intron > 0 & exon == 0, 1, 0))
```

Note that you can chain annotation. E.g. annotate the dataset with
exons/introns/genes, etc. Then, collapse all genes, expand by 10kb on either
side, and add another column. This package has a function `to_granges()` to
make this easier:

```r
> genes <- tgff %>% filter(feature=='gene') %>% to_granges() %>% reduce()
> genes <- reduce(genes) ## collapse overlapping genes
> start(genes) <- start(genes) - 10e3
> end(genes) <- end(genes) + 10e3
# if you want, you could use trim() if you have chromosome lengths.
> intragenic <- gaps(genes)
> intragenic$feature <- 'intragenic'
> da <- annotate_overlaps(da, intragenic)
```

Now we have annotated regions that are at least 10kb away from a gene region
annotated as intragenic.


## Future Goals

Currently the annotation is just number of overlaps to a particular type of
feature. Using list-columns, it's trivial to extend this such that there are
additional columns containing the ragged list of identifiers.
