# SpotGraphs <!-- badges: start -->  [![R-CMD-check](https://github.com/Sanin-Lab/SpotGraphs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Sanin-Lab/SpotGraphs/actions/workflows/R-CMD-check.yaml)  <!-- badges: end --> <a href="https://sanin-lab.github.io/SpotGraphs/"><img src="man/figures/logo.png" align="right" height="139" alt="SpotGraphs website" /></a>


## Functions to aid in Visium ST data analysis with igraph

### Current functions in this repo:
1. `SpotGraph(coord, dist.buffer = 1.05, max.dist = NULL)`
   - this will create an igraph object given a 2-column data frame or matrix of x,y coordinates
   - each node = a spot, where an edge exists between two spots if they are adjacent to each other
   - returns an igraph object
2. `CleanSlide(coord, nCount)`
   - takes a 2-column data frame or matrix of x,y coordinates and calls SpotGraph to create an igraph object
   - performs modularity maximization to group connected communities of spots
   - adds total transcripts together per-cluster (nCount)
   - calculates a threshold to identify low quality/small communities of spots
   - returns a data frame with per-spot cluster results, total counts per cluster, and whether the cluster passed the automatically detected threshold
3. `SpatialPlotGraph(igraph_object, coord = NULL, group.by, label = T)`
   - use the original x,y coordinates to plot some vertex_attr in the corresponding igraph object
5. `CutEdges(igraph_object, cluster_pairs = NULL, cluster.col = "cluster")`
   - remove all edges between specified groups of spots in an igraph object
