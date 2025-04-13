# Functions to aid in Visium ST data analysis with igraph

Current functions in this repo:
1. SpotGraph(coord, max.dist)
   - this will create an igraph object given a 2-column data frame or matrix
   - uses max.dist as the maximum distance allowed between two adjacent spots
   - each node = a spot, where an edge exists between two spots if they are adjacent to each other
2. CleanSlide(obj)
   - takes a Seurat object with Visium data and calls SpotGraph to create an igraph object
   - performs modularity maximization to group connected communities of spots
   - adds total transcripts together per-cluster (nCount_Spatial)
   - calculates a threshold to identify low quality/small communities of spots
