# Identify centers or boundaries of some region of interest

``` r
library(SpotGraphs)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
```

## Overview

A common task in spatial analysis is to perform some analysis with a
region of interest, such as a tumor. In this vignette, we demonstrate
how to use the
[`NeighborhoodCenters()`](https://sanin-lab.github.io/SpotGraphs/reference/NeighborhoodCenters.md)
function to help us identify characteristics such as the center of a
region of interest or its boundary.

## 1. Load data

We will again use the 10X Visium dataset, downloaded from GSE208253, to
demonstrate
[`NeighborhoodCenters()`](https://sanin-lab.github.io/SpotGraphs/reference/NeighborhoodCenters.md).

``` r
class(scc_s1)
#> [1] "Seurat"
#> attr(,"package")
#> [1] "SeuratObject"
dim(scc_s1)
#> Loading required namespace: SeuratObject
#> [1] 36601  1185

coord = Seurat::GetTissueCoordinates(scc_s1)
coord = coord[,c('x', 'y')]
```

## 2. Select a region of interest

In this vignette, we will use clustering to select some region of
interest in this tissue sample, but in practice, this should instead be
a region of some biological interest. The region of interest, i.e.,
neighborhood, is provided to
[`NeighborhoodCenters()`](https://sanin-lab.github.io/SpotGraphs/reference/NeighborhoodCenters.md)
as a named vector of `TRUE`/`FALSE` values, where the names correspond
to spot barcodes in the given dataset.

``` r
# create an igraph object with SpotGraph() and perform clustering
ig = SpotGraph(coord)
cl = igraph::cluster_fast_greedy(ig)
cl = factor(cl$membership)

# select cluster 1 and name the boolean vector with vertex names
roi = cl==1
names(roi) = names(igraph::V(ig))
```

## 3. Run `NeighborhoodCenters()`

``` r
res = NeighborhoodCenters(coord = coord, is_neighborhood = roi)
```

The output from this function returns a list, where the first element is
a data.frame with the scores from
[`igraph::centr_eigen()`](https://r.igraph.org/reference/centr_eigen.html)
where `is_neighborhood == T`, and 0 where `is_neighborhood == F`.

``` r
print(head(res$eigen.scores))
#>              barcode center_eigen
#> 1 AAACACCAATAACTGC-1 0.0000000000
#> 2 AAACAGGGTCTATATT-1 0.0006409795
#> 3 AAACCGTTCGTCCAGG-1 0.0000000000
#> 4 AAACGAGACGGTTGAT-1 0.0000000000
#> 5 AAACTGCTGGCTCCAA-1 0.0000000000
#> 6 AAAGACTGGGCGCTTT-1 0.0000000000

# Set centrality scores as a vertex attribute in our igraph object
igraph::V(ig)$centr_eigen = res$eigen.scores$center_eigen
igraph::V(ig)$roi_center = names(igraph::V(ig)) %in% res$centers
igraph::V(ig)$roi_boundary = names(igraph::V(ig)) %in% res$boundary

plt1 = SpatialPlotGraph(ig, group.by = 'centr_eigen',
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_viridis() +
  theme_void()
#> Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
#> ℹ Please use tidy evaluation idioms with `aes()`.
#> ℹ See also `vignette("ggplot2-in-packages")` for more information.
#> ℹ The deprecated feature was likely used in the SpotGraphs package.
#>   Please report the issue to the authors.
#> This warning is displayed once every 8 hours.
#> Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
#> generated.
plt2 = SpatialPlotGraph(ig, group.by = 'roi_center',
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_manual(values = c('grey', 'red')) +
  theme_void() +
  theme(legend.position = 'none')
plt3 = SpatialPlotGraph(ig, group.by = 'roi_boundary', 
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_manual(values = c('grey', 'red')) +
  theme_void()

wrap_plots(plt1, plt2, plt3, design = 'A#\nBC')
```

![](NeighborhoodCenters_files/figure-html/Plot%20eigenvector%20centrality%20scores-1.png)

## 4. Select multiple regions of interest

The
[`NeighborhoodCenters()`](https://sanin-lab.github.io/SpotGraphs/reference/NeighborhoodCenters.md)
function can also identify characteristics of multiple regions by
providing a named boolean vector indicating whether a spot is within any
of these regions. This can be useful when there are several tumor
regions within the same tissue section. Here, we select three different
clusters as regions of interest, instead of just one.

``` r
roi = cl %in% c('1', '2', '5')
names(roi) = names(igraph::V(ig))
res = NeighborhoodCenters(coord = coord, is_neighborhood = roi)
```

``` r
# Set centrality scores as a vertex attribute in our igraph object
igraph::V(ig)$centr_eigen = res$eigen.scores$center_eigen
igraph::V(ig)$roi_center = names(igraph::V(ig)) %in% res$centers
igraph::V(ig)$roi_boundary = names(igraph::V(ig)) %in% res$boundary

plt1 = SpatialPlotGraph(ig, group.by = 'centr_eigen', 
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_viridis() +
  theme_void()
plt2 = SpatialPlotGraph(ig, group.by = 'roi_center', 
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_manual(values = c('grey', 'red')) +
  theme_void() +
  theme(legend.position = 'none')
plt3 = SpatialPlotGraph(ig, group.by = 'roi_boundary', 
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_manual(values = c('grey', 'red')) +
  theme_void()

wrap_plots(plt1, plt2, plt3, design = 'A#\nBC')
```

![](NeighborhoodCenters_files/figure-html/multiple%20region%20plotting-1.png)

## 5. Calculate distance between each spot and ROI center

``` r
roi_dist = igraph::distances(ig, weights = NA)[,res$centers] %>% 
  apply(1, min)
boundary_dist = igraph::distances(ig, weights = NA)[,res$boundary] %>% 
  apply(1, min)
# boundary_dist = ifelse(!roi, boundary_dist, 0)

igraph::V(ig)$roi_dist = roi_dist
igraph::V(ig)$boundary_dist = boundary_dist

plt.dist = SpatialPlotGraph(ig, group.by = 'roi_dist',
                 pt.size = 0.5, linewidth = 0.1) +
  scale_color_distiller(palette = 'RdBu', direction = 1) +
  ggtitle('distance to ROI centers') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        legend.position = 'none')
plt.boundarydist = SpatialPlotGraph(ig, group.by = 'boundary_dist',
                 pt.size = 0.5, linewidth = 0.1) +
  scale_color_distiller(palette = 'RdBu', direction = 1) +
  ggtitle('distance to ROI boundaries') +
  guides(color = guide_colorbar(title = 'distance'),
         alpha = guide_none()) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

wrap_plots(plt.dist, plt.boundarydist, nrow = 1)
```

![](NeighborhoodCenters_files/figure-html/ROI%20distances-1.png)
