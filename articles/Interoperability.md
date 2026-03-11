# Seurat and SpatialExperiment interaction

``` r
library(SpotGraphs)
library(Seurat)
library(TENxVisiumData)
library(SpatialExperiment)
library(dplyr)
library(igraph)
```

## Overview

In this vignette, we demonstrate how SpotGraphs can interact with Seurat
and SpatialExperiment objects. A typical workflow will involve the
following steps:

1.  Extract x,y coordinates from Seurat or SpatialExperiment object.

2.  Create igraph object with
    [`SpotGraph()`](https://potential-adventure-or7z9q9.pages.github.io/reference/SpotGraph.md).

3.  Perform some analysis.

4.  Store results back into the original Seurat or SpatialExperiment
    object.

## Seurat object interaction

We will use the example Seurat object provided in this package to
demonstrate a typical workflow.

``` r
class(scc_s1)
#> [1] "Seurat"
#> attr(,"package")
#> [1] "SeuratObject"
```

First, extract the x,y coordinates

``` r
coord = Seurat::GetTissueCoordinates(scc_s1)
```

Create an igraph object with
[`SpotGraph()`](https://potential-adventure-or7z9q9.pages.github.io/reference/SpotGraph.md)

``` r
ig = SpotGraph(coord[,c('x','y')])
```

Perform some analysis. In this case, we will calculate the number of
edges per spot.

``` r
n_edges = igraph::degree(ig)
```

Finally, store these results back into the metadata of the original
Seurat object

``` r
scc_s1 = AddMetaData(scc_s1, metadata = n_edges, col.name = 'degree')
```

## SpatialExperiment interaction

To demonstrate interaction with a SpatialExperiment object, we will load
an object from ExperimentHub.

``` r
eh = ExperimentHub::ExperimentHub()
q = AnnotationHub::query(eh, 'TENxVisium')
id = q$ah_id[q$title=='HumanGlioblastoma_v3.13']
spe = eh[[id]]
class(spe)
#> [1] "SpatialExperiment"
#> attr(,"package")
#> [1] "SpatialExperiment"
```

Once we have our SpatialExperiment object, we can follow the same
workflow

``` r
coord = as.data.frame(spatialCoords(spe))
ig = SpotGraph(coord)
n_edges = igraph::degree(ig)
colData(spe)$degree = n_edges
```

We could also calculate multiple statistics from several analyses on the
igraph object before storing them in the original
Seurat/SpatialExperiment object. If we store all the vertex/spot-level
results back into the igraph object with
[`V()`](https://r.igraph.org/reference/V.html), we can extract all of
the results at once with
[`igraph::vertex_attr()`](https://r.igraph.org/reference/vertex_attr.html),
and coerce this into a `data.frame`.

``` r
V(ig)$cluster = factor(cluster_louvain(ig)$membership)
V(ig)$ecc = eccentricity(ig)
V(ig)$centrality = centr_eigen(ig)$vector

df = as.data.frame(vertex_attr(ig))
head(df)
#>                 name coord_x coord_y is_boundary cluster ecc   centrality
#> 1 AAACAAGTATCTCCCA-1    7438    9389       FALSE       1  69 3.244117e-01
#> 2 AAACAATCTACTAGCA-1    1802    5327        TRUE       2   6 1.885746e-16
#> 3 AAACACCAATAACTGC-1    8514    3670        TRUE       3  75 3.573227e-02
#> 4 AAACAGAGCGACTCCT-1    3122    8840       FALSE       4  64 8.092829e-02
#> 5 AAACAGTGTTCCTGGG-1   10193    5323       FALSE       5  74 1.289729e-01
#> 6 AAACATTTCCCGGATT-1    8756    9044       FALSE       6  72 3.142856e-01
```

This resulting `data.frame` and then be stored back into the
SpatialExperiment object.

``` r
colData(spe) = cbind(colData(spe), df)
head(colData(spe))
#> DataFrame with 6 rows and 12 columns
#>                                 sample_id in_tissue array_row array_col
#>                               <character> <logical> <integer> <integer>
#> AAACAAGTATCTCCCA-1 HumanGlioblastoma_v3..      TRUE        50       102
#> AAACAATCTACTAGCA-1 HumanGlioblastoma_v3..      TRUE         3        43
#> AAACACCAATAACTGC-1 HumanGlioblastoma_v3..      TRUE        59        19
#> AAACAGAGCGACTCCT-1 HumanGlioblastoma_v3..      TRUE        14        94
#> AAACAGTGTTCCTGGG-1 HumanGlioblastoma_v3..      TRUE        73        43
#> AAACATTTCCCGGATT-1 HumanGlioblastoma_v3..      TRUE        61        97
#>                       degree               name   coord_x   coord_y is_boundary
#>                    <numeric>        <character> <integer> <integer>   <logical>
#> AAACAAGTATCTCCCA-1         6 AAACAAGTATCTCCCA-1      7438      9389       FALSE
#> AAACAATCTACTAGCA-1         2 AAACAATCTACTAGCA-1      1802      5327        TRUE
#> AAACACCAATAACTGC-1         4 AAACACCAATAACTGC-1      8514      3670        TRUE
#> AAACAGAGCGACTCCT-1         6 AAACAGAGCGACTCCT-1      3122      8840       FALSE
#> AAACAGTGTTCCTGGG-1         6 AAACAGTGTTCCTGGG-1     10193      5323       FALSE
#> AAACATTTCCCGGATT-1         6 AAACATTTCCCGGATT-1      8756      9044       FALSE
#>                     cluster       ecc  centrality
#>                    <factor> <numeric>   <numeric>
#> AAACAAGTATCTCCCA-1        1        69 3.24412e-01
#> AAACAATCTACTAGCA-1        2         6 1.88575e-16
#> AAACACCAATAACTGC-1        3        75 3.57323e-02
#> AAACAGAGCGACTCCT-1        4        64 8.09283e-02
#> AAACAGTGTTCCTGGG-1        5        74 1.28973e-01
#> AAACATTTCCCGGATT-1        6        72 3.14286e-01
```
