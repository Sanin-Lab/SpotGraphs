# Identify centers or boundaries

``` r
library(SpotGraphs)
library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)
library(patchwork)
```

## Overview

A common task in spatial analysis is to perform some analysis with a
region of interest, such as a tumor. In this vignette, we demonstrate
how to use the
[`NeighborhoodCenters()`](https://potential-adventure-or7z9q9.pages.github.io/reference/NeighborhoodCenters.md)
function to help us identify characteristics such as the center of a
region of interest or its boundary.

## 1. Load data

We will again use the 10X Visium dataset, downloaded from GSE208253, to
demonstrate
[`NeighborhoodCenters()`](https://potential-adventure-or7z9q9.pages.github.io/reference/NeighborhoodCenters.md).

``` r
data("scc_s1", package = "SpotGraphs")
scc_s1 = UpdateSeuratObject(scc_s1)
#> Validating object structure
#> Updating object slots
#> Ensuring keys are in the proper structure
#> Ensuring keys are in the proper structure
#> Ensuring feature names don't have underscores or pipes
#> Updating slots in Spatial
#> Updating slots in slice1
#> Warning: Not validating Centroids objects
#> Updated Centroids object 'centroids' in FOV 'slice1'
#> Updated boundaries in FOV 'slice1'
#> Validating object structure for Assay5 'Spatial'
#> Validating object structure for VisiumV2 'slice1'
#> Object representation is consistent with the most current Seurat version
class(scc_s1)
#> [1] "Seurat"
#> attr(,"package")
#> [1] "SeuratObject"
dim(scc_s1)
#> [1] 36601  1185

coord = Seurat::GetTissueCoordinates(scc_s1)
coord = data.frame(x = coord$y, y = -coord$x)
head(coord)
#>       x      y
#> 1 16571  -4809
#> 2 13546  -3944
#> 3 14812  -8142
#> 4 10536 -13505
#> 5 13053 -11764
#> 6  9011  -4240
```

## 2. Select a region of interest

In this vignette, we will use clustering to select some region of
interest in this tissue sample, but in practice, this should instead be
a region of some biological interest. The region of interest, i.e.,
neighborhood, is provided to
[`NeighborhoodCenters()`](https://potential-adventure-or7z9q9.pages.github.io/reference/NeighborhoodCenters.md)
as a named vector of `TRUE`/`FALSE` values, where the names correspond
to spot barcodes in the given dataset.

``` r
# create an igraph object with SpotGraph() and perform clustering
ig = SpotGraph(coord)
cl = igraph::cluster_fast_greedy(ig)
cl = factor(cl$membership)

# select cluster 1 and name the boolean vector with vertex names
roi = cl==5
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
#>   barcode center_eigen
#> 1       1            0
#> 2       2            0
#> 3       3            0
#> 4       4            0
#> 5       5            0
#> 6       6            0

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

![](NeighborhoodCenters_files/figure-html/Plot%20eigenvector%20centrality%20scores-1.png)

## 4. Select multiple regions of interest

The
[`NeighborhoodCenters()`](https://potential-adventure-or7z9q9.pages.github.io/reference/NeighborhoodCenters.md)
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
