#' Creates an igraph object given x,y coordinates
#' @description
#' Create an igraph object with x,y coordinates, assuming all adjacent
#' spots are equidistant from each other, and drawing edges where spots
#' are immediately adjacent to each other.
#'
#' @param coord A two-column data.frame or matrix, where each column contains x or y coordinates.
#' @param delaunay Whether to use Delaunay triangulation to identify edges for network construction.
#' Default is FALSE.
#' @param delaunay.trim If TRUE (default) and delaunay = TRUE, then edge lengths longer than max.dist
#'  will be removed after Delaunay triangulation.
#' @param dist.buffer Influences the maximum distance a node is allowed to be from another
#' node to be considered a neighbor. Only used if `max.dist = NULL`.
#' @param max.dist The furthest a node should be from another node to be considered neighbors.
#' If NULL (default), the shortest distance between any two nodes is the minimum distance (`min.dist`)
#' required for two nodes to be considered neighbors, and the `max.dist` is calculated as the hypotenuse
#' of a triangle with two sides of length `min.dist`, and then is multiplied by `dist.buffer`.
#
#' @return an igraph object, where each vertex (i.e., node) corresponds to each
#' row in the `coord` input, with un-weighted edges between vertices that are
#' immediately adjacent to each other.
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a coordinate data frame with two isolated groups of points
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#'
#' # Preview the points on a grid with their corresponding x,y coordinates
#' ggplot(df, aes(x = x,y=y)) +
#'   geom_point()
#'
#' # Identify points that are immediately adjacent to each other and create an igraph object
#' ig = SpotGraph(df)
#'
#' # Optionally view the network with ggnetwork
#' # library(ggnetwork)
#' ggplot(ig, aes(x=x, y=y, xend=xend, yend=yend)) +
#'   geom_edges() +
#'   geom_nodes()
#' }
SpotGraph = function(coord,
                     delaunay = F,
                     delaunay.trim = T,
                     dist.buffer = 1.05,
                     max.dist = NULL) {
  # Get coordinates and calculate euclidean distance
  d = dist(coord, method = 'euclidean')
  m = as.matrix(d)

  # Calculate max tolerated distance to define immediately adjacent spots
  min.dist = min(d)
  if (is.null(max.dist)) max.dist = sqrt(2*min.dist^2)*dist.buffer

  # Use Delaunay triangulation if desired
  if (delaunay) {
    tr.obj = interp::tri.mesh(x = coord)

    # extract edge list and create igraph object
    edges = tr.obj$arcs
    ig = igraph::graph_from_edgelist(edges, directed = F)

    # Add x,y coordinates to igraph object
    igraph::V(ig)$name = rownames(coord)

    # remove long edges
    if (delaunay.trim) {
      edge.df = as_ids(E(ig)) %>%
        stringr::str_split('\\|') %>%
        do.call(what = rbind) %>%
        as.data.frame()
      edge.lengths = apply(edge.df, 1, function(xx) m[xx[1],xx[2]])

      edge.del = edge.df[edge.lengths>max.dist,] %>%
        apply(1, paste, collapse = '|')

      ig = igraph::delete_edges(ig, edge.del)
    }
  } else {
    # Create igraph object from edge data frame
    ig = igraph::graph_from_adjacency_matrix(m <= max.dist, mode = "undirected", diag = F)
  }

  # Add x,y coordinates to igraph object
  igraph::V(ig)$coord_x = coord[,1]
  igraph::V(ig)$coord_y = coord[,2]

  # identify boundary nodes
  is_boundary = igraph::degree(ig) < max(igraph::degree(ig))
  igraph::V(ig)$is_boundary = is_boundary

  return(ig)
}


