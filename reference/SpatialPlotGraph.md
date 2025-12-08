# Visualize igraph vertex attribute on tissue x,y coordinates

Visualize igraph vertex attribute on tissue x,y coordinates

## Usage

``` r
SpatialPlotGraph(
  igraph_object = NULL,
  coord = NULL,
  group.by = "is_boundary",
  label = FALSE,
  flip.axes = TRUE,
  linewidth = 0.5,
  pt.size = 1.6
)
```

## Arguments

- igraph_object:

  An igraph object

- coord:

  (optional) A two-column data.frame or matrix, where each column
  contains x or y coordinates. if not provided, will look for `coord_x`
  and `coord_y` attributes in the provided igraph object.

- group.by:

  Some vertex attribute to plot onto tissue coordinates. Must be present
  in `vertex_attr(ig)`.

- label:

  Logical, whether to label groups on the plot based on `group.by`.

- flip.axes:

  Default is TRUE. Invert the x and y axes, typically needed with
  coordinates from Visium Seurat objects to align the spots with the
  tissue image.

- linewidth:

  Numeric input to determine width of segments between nodes. Used in
  `geom_segment().`

- pt.size:

  Numeric input to determine size of each point on the plot. Used as the
  size parameter in `geom_point()`.

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
SpatialPlotGraph(igraph_object = ig, group.by = 'is_boundary')
} # }
```
