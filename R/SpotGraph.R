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

  coord_mat <- as.matrix(coord)
  n_points  <- nrow(coord_mat)

  # Use Delaunay triangulation if desired
  if (delaunay) {
    tr.obj = interp::tri.mesh(x= jitter(coord_mat[, 1],0.05), y = jitter(coord_mat[, 2],0.05))

    # extract edge list and create igraph object
    edges = tr.obj$arcs
    ig = igraph::graph_from_edgelist(edges, directed = F)

    # Add x,y coordinates to igraph object
    igraph::V(ig)$name = rownames(coord)

    # remove long edges
    if (delaunay.trim) {
      #extract edges directly as a 2 column matrix
      e_mat = igraph::as_edgelist(ig)

      #vectorized distance calculation for only connected edges
      p1 = coord_mat[e_mat[, 1], , drop = F]
      p2 = coord_mat[e_mat[, 2], , drop = F]
      edge.lengths = sqrt(rowSums((p1 - p2)^2))

      #calculate max.dist
      if (is.null(max.dist)) {
        min.dist = min(edge.lengths)
        max.dist = sqrt(2 * min.dist^2) * dist.buffer
      }

      ig = igraph::delete_edges(ig, igraph::E(ig)[edge.lengths > max.dist])
    }
  } else {
    # Get coordinates and calculate euclidean distance
    #cap k dynamically
    knn_result <- dbscan::kNN(as.matrix(coord), k = min(50, n_points - 1))

    # Calculate max tolerated distance to define immediately adjacent spots
    min.dist = min(knn_result$dist)
    if (is.null(max.dist)) max.dist = sqrt(2*min.dist^2)*dist.buffer

    # Create igraph object from edge data frame
    n <- nrow(knn_result$id)
    k <- ncol(knn_result$id)

    # Vectorized operations
    from <- rep(1:n, times = k)
    to <- as.vector(knn_result$id)
    dist <- as.vector(knn_result$dist)

    # Single mask operation
    mask <- dist <= max.dist & from != to

    # Direct edge list creation
    ig <- igraph::make_undirected_graph(
      edges = rbind(from[mask], to[mask]),
      n = n
    )

    # Remove duplicates
    ig <- igraph::simplify(ig, remove.multiple = TRUE)
    }

  # Add x,y coordinates to igraph object
  igraph::V(ig)$coord_x = coord_mat[, 1]
  igraph::V(ig)$coord_y = coord_mat[, 2]

  # identify boundary nodes
  is_boundary = igraph::degree(ig) < max(igraph::degree(ig))
  igraph::V(ig)$is_boundary = is_boundary

  return(ig)
}


