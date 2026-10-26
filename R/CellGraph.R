#' Creates an igraph object given x,y coordinates
#' @description
#' Create an igraph object with x,y coordinates, drawing edges where cells
#' are adjacent to each other based on a calculated threshold.
#'
#' @param coord A two-column data.frame or matrix, where each column contains x or y coordinates.
#' @param dist.buffer Influences the maximum distance a node is allowed to be from another
#' node to be considered a neighbor.
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
#' # Identify points that are immediately adjacent to each other and create an igraph object
#' ig = CellGraph(df)
#' }
CellGraph <- function(
    coord, dist.buffer = 1.05){

  coord_mat <- as.matrix(coord) #convert to matrix
  n_points <- nrow(coord_mat)

  knn_result <- dbscan::kNN(coord_mat, k = min(50, n_points - 1))  #Memory efficient ED calculation

  dens <- density(knn_result$dist) #Smooth distances
  peaks <- which(diff(sign(diff(dens$y))) == -2) + 1 #find local maxima
  mean.dist <- dens$x[peaks][1] #find first peak
  max.dist = sqrt(2 * mean.dist^2) * dist.buffer #calculate distance threshold

  #Generate graph
  n <- nrow(knn_result$id)
  k <- ncol(knn_result$id)
  from <- rep(1:n, times = k)
  to <- as.vector(knn_result$id)
  dist <- as.vector(knn_result$dist)
  mask <- dist <= max.dist & from != to
  ig <- igraph::make_undirected_graph(edges = rbind(from[mask],to[mask]), n = n)
  ig <- igraph::simplify(ig, remove.multiple = TRUE)

  igraph::V(ig)$coord_x = coord_mat[, 1]
  igraph::V(ig)$coord_y = coord_mat[, 2]

  return(ig)
}

