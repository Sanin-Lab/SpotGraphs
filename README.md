# SpotGraphs [![R-CMD-check](https://github.com/Sanin-Lab/SpotGraphs/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Sanin-Lab/SpotGraphs/actions/workflows/R-CMD-check.yaml) <a href="https://sanin-lab.github.io/SpotGraphs/"><img src="man/figures/logo.png" align="right" height="139" alt="SpotGraphs website" /></a>


## Functions to aid in spatial transcriptomics data analysis

## Motivation
Current spatial transcriptomic analysis pipelines in R focus on pre-processing and visualization, while providing limited tools to interact with the "spatial" aspect of this data. To address this limitation, we provide a set of tools that allow more flexibility when working with the spatial coordinates to enable further filtering of low quality spots on tissue debris, edit spot-level adjacencies, and identify centers or boundaries of user-defined neighborhoods of interest.

This package allows users to simply provide the x,y-coordinates of their spatial data to create an igraph object with the `SpotGraph()` function, from which various graph-based statistics can be calculated and stored as meta data in the user's original Seurat or SpatialExperiment object. 

## Method overview
To construct a graph given x,y-coordinates of spatial data, the SpotGraph() function implements two approaches to identifying neighboring spots to build an adjacency matrix, either based on Euclidean distance or with Delaunay triangulation. See our manuscript for more details. 

## Installation
To install the latest version of our package, run:
`devtools::install_github("Sanin-Lab/SpotGraphs")`

## Usage
To create an igraph object from spatial data, the only required input into the `SpotGraph()` function is a data frame or matrix with two columns corresponding to x,y-coordinates. 
`SpotGraph(coord.df)`

## Tutorials and Applications
See documentation website for in-depth walkthroughs.
