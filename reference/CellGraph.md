# Creates an igraph object given x,y coordinates

Create an igraph object with x,y coordinates, drawing edges where cells
are adjacent to each other based on a calculated threshold.

## Usage

``` r
CellGraph(coord, dist.buffer = 1.05)
```

## Arguments

- coord:

  A two-column data.frame or matrix, where each column contains x or y
  coordinates.

- dist.buffer:

  Influences the maximum distance a node is allowed to be from another
  node to be considered a neighbor.

## Value

an igraph object, where each vertex (i.e., node) corresponds to each row
in the `coord` input, with un-weighted edges between vertices that are
immediately adjacent to each other.

## Examples

``` r
if (FALSE) { # \dontrun{
# Create a coordinate data frame with two isolated groups of points
df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
colnames(df) = c('x', 'y')

# Identify points that are immediately adjacent to each other and create an igraph object
ig = CellGraph(df)
} # }
```
