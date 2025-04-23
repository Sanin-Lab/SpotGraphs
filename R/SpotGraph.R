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
#' # Create a coordinate data frame with two isolated groups of points
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#'
#' # Preview the points on a grid with their corresponding x,y coordinates
#' ggplot(df, aes(x = x,y=y)) + geom_point()
#'
#' # Identify points that are immediately adjacent to each other and create an igraph object
#' ig = SpotGraph(df)
#'
#' # Optionally view the network with ggnetwork
#' library(ggnetwork)
#' ggplot(ig, aes(x=x, y=y, xend=xend, yend=yend)) +
#'   geom_edges() +
#'   geom_nodes()

SpotGraph = function(coord, cluster = F, resolution = 0.5) {
  # Get coordinates and calculate euclidean distance
  d = dist(coord, method = 'euclidean')
  m = as.matrix(d)

  # Find the distance between two adjacent spots
  # - assumes all adjacent spots are equidistant
  max.dist = min(d) + 1

  # Identify whether a spot is immediately neighboring another spot
  # - this algorithm seems robust to a range of distance thresholds
  #   to identify neighbors, as long as this threshold is larger than
  #   the minimum distance between spots.
  neighbors = apply(m < max.dist, 1, function(spot) colnames(m)[spot])

  # Create edge data frame
  # - two columns, each column indicating end nodes
  # - each spot is a node in ST data
  edges = lapply(names(neighbors), function(xx) {
    data.frame(adj.spot = neighbors[xx] %>%
                 unlist(use.names = F),
               spot = names(neighbors[xx]))
  }) %>% do.call(what = rbind) %>%
    filter(spot != adj.spot)

  # Remove duplicated edges
  edges$edgeid = apply(edges, 1, function(xx) {
    str_flatten(sort(c(xx[1], xx[2])))
  })
  edges = edges %>%
    filter(!duplicated(edgeid)) %>%
    select(!matches('edgeid')) %>%
    as.matrix

  # Create igraph object from edge data frame
  ig = graph_from_edgelist(edges, directed = F)

  if (cluster) {
    cluster_res = factor(cluster_louvain(ig, resolution = resolution)$membership)
    ig = set_vertex_attr(ig, name = 'iglouvain_cluster', value = cluster_res)
  }

  return(ig)
}


