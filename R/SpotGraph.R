#' Creates an igraph object given x,y coordinates
#' @description
#' Create an igraph object with x,y coordinates, assuming all adjacent
#' spots are equidistant from each other, and drawing edges where spots
#' are immediately adjacent to each other.
#'
#' @param coord A two-column data.frame or matrix, where each column contains x or y coordinates.
#' @param cluster Whether clustering should be performed on the igraph object.
#' @param resolution If cluster = T, a higher number will return more clusters and
#' a lower number will return fewer clusters. Clustering results are stored in the
#' returned igraph object under the name 'iglouvain_cluster' and can be accessed
#' by vertex_attr(igraph, 'iglouvain_cluster').
#
#' @return an igraph object, where each vertex (i.e., node) corresponds to each
#' row in the coord input, with un-weighted edges between vertices that are
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


SpotGraph = function(coord, dist.buffer = 1.05, max.dist = NULL,
                     cluster = F, resolution = 0.5) {
  # Get coordinates and calculate euclidean distance
  d = dist(coord, method = 'euclidean')
  m = as.matrix(d)

  # Find the distance between two adjacent spots
  # - buffer the min.dist to define adjacency by dist.buffer (default 5%)
  # - buffer the max.dist by hypotenuse if max.dist is not provided by user
  # - a spot that is immediately adjacent needs to be within min.dist and max.dist
  min.dist = min(d)
  if (is.null(max.dist)) max.dist = sqrt(2*min.dist^2)*dist.buffer

  # Create igraph object from edge data frame
  ig = igraph::graph_from_adjacency_matrix(m <= max.dist, mode = "undirected", diag = F)

  # Add x,y coordinates to igraph object
  igraph::V(ig)$coord_x = coord[,1]
  igraph::V(ig)$coord_y = coord[,2]

  # identify boundary nodes
  is_boundary = degree(ig) < max(degree(ig))
  igraph::V(ig)$is_boundary = is_boundary

  # calculate edge weights based on connections to boundary nodes
  boundary_nodes = names(is_boundary)[is_boundary]
  weights = as_ids(E(ig)) %>% str_split('\\|') %>%
    sapply(function(nodes) {
      xx = ifelse(nodes %in% boundary_nodes, 0.5, 1)
      return(prod(xx))
    })

  # Add weights to graph
  igraph::E(ig)$weight <- weights

  # if (cluster) {
  #   # perform clustering
  #   cluster_res = igraph::cluster_louvain(ig, weights = weights, resolution = resolution)
  #   cluster_res = factor(cluster_res$membership)
  #   ig = igraph::set_vertex_attr(ig, name = 'iglouvain_cluster', value = cluster_res)
  # }

  return(ig)
}


