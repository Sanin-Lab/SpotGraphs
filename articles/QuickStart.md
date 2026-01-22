# QuickStart

``` r
library(SpotGraphs)
library(TENxVisiumData)
library(SpatialExperiment)
library(dplyr)
```

## Load example dataset

``` r
# retrieve data from TENxVisium
eh = ExperimentHub::ExperimentHub()
q = AnnotationHub::query(eh, 'TENxVisium')
id = q$ah_id[q$title=='HumanGlioblastoma_v3.13']
spe = eh[[id]]
```

## Filter stray spots

``` r
# 1. extract spot x,y-coordinates and create igraph object
coord.xy = as.data.frame(spatialCoords(spe))
ig = SpotGraph(coord.xy)

# 2. perform clustering on graph to identify communities of connected spots
# - store results in igraph object with V() <-
cl = igraph::cluster_louvain(ig, resolution = 0)
igraph::V(ig)$cluster = factor(cl$membership)

# 3. plot clustering results to identify which spots to keep or remove
plt1 = SpatialPlotGraph(ig, group.by = 'cluster', pt.size = 1) +
  theme_void() +
  theme(legend.position = 'none')
plt2 = SpatialPlotGraph(ig, group.by = 'cluster', pt.size = 1, label = T) +
  theme_void() +
  theme(legend.position = 'none')

wrap_plots(plt1, plt2)

# 3. store results in SpatialExperiment object
colData(spe)$igraph_cluster = V(ig)$cluster

# 4. keep spots in igraph clusters 1 and 2
spe = spe[,spe$igraph_cluster %in% c('1', '2')]
```
