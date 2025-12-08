# Creates an igraph object given x,y coordinates

Create an igraph object with x,y coordinates, assuming all adjacent
spots are equidistant from each other, and drawing edges where spots are
immediately adjacent to each other.

## Usage

``` r
SpotGraph(
  coord,
  delaunay = F,
  delaunay.trim = T,
  dist.buffer = 1.05,
  max.dist = NULL
)
```

## Arguments

- coord:

  A two-column data.frame or matrix, where each column contains x or y
  coordinates.

- delaunay:

  Whether to use Delaunay triangulation to identify edges for network
  construction. Default is FALSE.

- delaunay.trim:

  If TRUE (default) and delaunay = TRUE, then edge lengths longer than
  max.dist will be removed after Delaunay triangulation.

- dist.buffer:

  Influences the maximum distance a node is allowed to be from another
  node to be considered a neighbor. Only used if `max.dist = NULL`.

- max.dist:

  The furthest a node should be from another node to be considered
  neighbors. If NULL (default), the shortest distance between any two
  nodes is the minimum distance (`min.dist`) required for two nodes to
  be considered neighbors, and the `max.dist` is calculated as the
  hypotenuse of a triangle with two sides of length `min.dist`, and then
  is multiplied by `dist.buffer`.

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

# Preview the points on a grid with their corresponding x,y coordinates
ggplot(df, aes(x = x,y=y)) +
  geom_point()

# Identify points that are immediately adjacent to each other and create an igraph object
ig = SpotGraph(df)

# Optionally view the network with ggnetwork
# library(ggnetwork)
ggplot(ig, aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_edges() +
  geom_nodes()
} # }
```
