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
#> Warning: Not validating Centroids objects
class(scc_s1)
#> [1] "Seurat"
#> attr(,"package")
#> [1] "SeuratObject"
dim(scc_s1)
#> [1] 36601  1185

coord = Seurat::GetTissueCoordinates(scc_s1)
coord = data.frame(x = coord$x, y = -coord$y)
```

Perform some spot-level filtering (in the same fashion as in the [Get
Started](https://potential-adventure-or7z9q9.pages.github.io/articles/SpotGraphs.Rmd)
vignette.

``` r
ig = SpotGraph(coord)
cl = igraph::cluster_louvain(ig, resolution = 0)
scc_s1 = AddMetaData(scc_s1, cl$membership %in% c('1','2'), 'cl_threshold')
scc_s1 = subset(scc_s1, subset = cl_threshold)
```

Update the igraph object after filtering the Seurat object

``` r
coord = GetTissueCoordinates(scc_s1)
coord = data.frame(x = coord$x, y = -coord$y)
ig = SpotGraph(coord)
```

## 2. Select a region of interest

In this vignette, we will use clustering to select some region of
interest in this tissue sample, but in practice, this should instead be
a region of some biological interest. The region of interest, i.e.,
neighborhood, is provided to
[`NeighborhoodCenters()`](https://potential-adventure-or7z9q9.pages.github.io/reference/NeighborhoodCenters.md)
as a named vector of `TRUE`/`FALSE` values.

``` r
# perform clustering on the igraph object
cl = igraph::cluster_louvain(ig, resolution = 0.4)

# store clustering results back into igraph object
# - note: V() will access all vertices of an igraph object. 
#        `V() <-` also allows us to store vertex-level attributes 
igraph::V(ig)$cluster = factor(cl$membership)

# select cluster 5 and name the boolean vector with vertex names
roi = cl$membership==5
igraph::V(ig)$is_roi = roi

# visualize clusters and region selection
plt1 = SpatialPlotGraph(ig, group.by = 'cluster', label = T,
                        pt.size = 0.5, linewidth = 0.1) +
  theme_void() +
  theme(legend.position = 'none')
plt2 = SpatialPlotGraph(ig, group.by = 'is_roi', label = F,
                        pt.size = 0.5, linewidth = 0.1) +
  scale_color_manual(values = c('grey','red')) +
  theme_void()

wrap_plots(plt1, plt2)
```

![](NeighborhoodCenters_files/figure-html/Choose%20a%20region-1.png)

## 3. Run `NeighborhoodCenters()`

``` r
res = NeighborhoodCenters(coord = coord, is_neighborhood = roi)
```

The output from this function returns a list, where the first element is
a data.frame with the scores from
[`igraph::centr_eigen()`](https://r.igraph.org/reference/centr_eigen.html)
where `is_neighborhood == T`, and 0 where `is_neighborhood == F`.

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
roi = cl$membership %in% c('2', '5', '7')
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
roi_dist = igraph::distances(ig)[,res$centers] %>% apply(1, min)
boundary_dist = igraph::distances(ig)[,res$boundary] %>% apply(1, min)

# store calculated distances back into igraph objects
igraph::V(ig)$roi_dist = roi_dist
igraph::V(ig)$boundary_dist = boundary_dist

# visualize distances between centers or boundaries
plt.dist = SpatialPlotGraph(ig, group.by = 'roi_dist',
                 pt.size = 0.5, linewidth = 0.1) +
  scale_color_distiller(palette = 'RdBu', direction = -1) +
  ggtitle('distance to ROI centers') +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
        legend.position = 'none')
plt.boundarydist = SpatialPlotGraph(ig, group.by = 'boundary_dist',
                 pt.size = 0.5, linewidth = 0.1) +
  scale_color_distiller(palette = 'RdBu', direction = -1) +
  ggtitle('distance to ROI boundaries') +
  guides(color = guide_colorbar(title = 'distance'),
         alpha = guide_none()) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold'))

wrap_plots(plt.dist, plt.boundarydist, nrow = 1)
```

![](NeighborhoodCenters_files/figure-html/ROI%20distances-1.png)
