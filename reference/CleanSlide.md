# Identify spots isolated from the largest groups of immediately adjacent spots

Identify spots in a 10X Visium Seurat object that are potentially low
quality. Mostly useful for slides where there are large areas without
tissue and standard spot-filtering based on nCounts is insufficient.

1.  Identify neighboring spots on x,y coordinates with SpotGraph()

2.  Performing modularity maximization to identify clusters

3.  Calculate total number of transcripts detected in each cluster
    (nCount)

4.  Set a threshold to identify which clusters should be filtered out

## Usage

``` r
CleanSlide(coord, nCount)
```

## Arguments

- coord:

  A two-column data.frame or matrix, where each column contains x or y
  coordinates.

- nCount:

  A named vector of transcript counts per spot, where names are spot
  barcodes/ids, and each element should match each row of `coord`.

## Value

A data.frame with the following columns:

1.  ig_cluster: igraph modularity maximization results

2.  cluster_nCount: summed transcript counts within each 'ig_cluster'
    group

3.  threshold: whether the respective 'ig_cluster' group passed the
    automatically determiend threshold

4.  barcode: the spot id or barcode sequence

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a coordinate data frame with two isolated groups of points
df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
colnames(df) = c('x', 'y')
nspots = dim(df)[1]
rownames(df) = paste0('spot', 1:nspots)

# Create some vector of transcript counts
spotcounts = rbinom(25, 10000, 0.3)
spotcounts = c(spotcounts, rbinom(9, 100, 0.3))
names(spotcounts) = rownames(df)

# The resulting data.frame should indicate that the smaller
# cluster of 3x3 points did not pass the threshold, while the
# larger group of 5x5 points should have passed the threshold.
res = CleanSlide(coord = df, nCount = spotcounts)
} # }
```
