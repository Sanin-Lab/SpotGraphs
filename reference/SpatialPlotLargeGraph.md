# Visualize a large graph on tissue x,y coordinates

Visualize a large graph on tissue x,y coordinates

## Usage

``` r
SpatialPlotLargeGraph(
  igraph_object = NULL,
  only_edges = T,
  aspec.ratio = 1.6,
  group.by = NULL
)
```

## Arguments

- igraph_object:

  An igraph object

- only_edges:

  Boolean to determine if only the edges of the graph should be
  returned. Default is `TRUE`.

- group.by:

  Some vertex attribute to plot onto tissue coordinates. Must be present
  in `vertex_attr(ig)`.

- aspect.ratio:

  Numeric input to constrain aspect ratio of the generated plot.

## Value

a ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a coordinate data frame and create and igraph object
df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
colnames(df) = c('x', 'y')
ig = SpotGraph(df)

# Plot cluster results
SpatialPlotLargeGraph(igraph_object = ig)
} # }
```
