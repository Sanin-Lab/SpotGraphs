# Identify center of some neighborhood of spots, e.g., the center of a tumor

Score spots using igraph::centr_eigen to identify spots in the center of
some region/neighborhood.

## Usage

``` r
NeighborhoodCenters(igraph_object = NULL, is_neighborhood)
```

## Arguments

- igraph_object:

  An igraph object with named vertices as an attribute called `name`

- is_neighborhood:

  A boolean vector (`TRUE`/`FALSE`) of equal length as the number of
  nodes in `igraph_object`. Centers will be identified for each
  neighborhood of spots labeled `TRUE`.

## Value

a list with the following elements:

1.  `eigen.scores`: a numeric vector with the output from
    igraph::centr_eigen, with a score for each spot, where
    `max(eigen.scores)` within each neighborhood will identify its
    center.

2.  `centers`: a character vector with barcode ids of spots that were
    identified as the center of each neighborhood of spots.

3.  `boundaries`: a character vector with barcode ids of spots that were
    identified be on the boundary of each neighborhood of spots.

## Examples

``` r
if (FALSE) { # \dontrun{
df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
colnames(df) = c('x', 'y')
nspots = dim(df)[1]
rownames(df) = paste0('spot', 1:nspots)
ig = SpotGraph(df)
igraph::V(ig)$name = rownames(df)

igraph::V(ig)$cluster = igraph::cluster_leiden(ig, resolution = 0)$membership
is_cluster = igraph::V(ig)$cluster=='1'
names(is_cluster) = names(igraph::V(ig))

res = NeighborhoodCenters(igraph_object = ig, is_neighborhood = is_cluster)
} # }
```
