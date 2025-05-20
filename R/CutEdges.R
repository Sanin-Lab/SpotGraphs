#' Remove all edges that lie between pairs of clusters
#'
#' @param igraph_object An igraph object with cluster results stored.
#' @param cluster_pairs A list containing two-element vectors indicating
#' pairs of cluster names where we want to remove edges between.
#' @param cluster.col A character value with the vertex attribute name in
#' the igraph object with clustering results.
#'
#' @return an updated igraph object with some edges removed
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a coordinate data frame
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#'
#' # Create an igraph object with these coordinates
#' ig = SpotGraph(df)
#'
#' # Perform clustering and add those cluster assignments to the igraph object
#' cl = igraph::cluster_louvain(ig)$membership
#' ig = igraph::set_vertex_attr(ig, 'cluster', value = factor(cl))
#'
#' # Examine graph before removing edges
#' SpatialPlotGraph(ig, group.by = 'cluster')
#'
#' # Remove edges between clusters 1 and 2, and between clusters 3 and 4
#' ig = CutEdges(igraph_object = ig,
#'               cluster_pairs = list(c(1,2), c(1,3)),
#'               cluster.col = 'cluster')
#'
#' # Examine graph after removing edges
#' SpatialPlotGraph(ig, group.by = 'cluster')
#' }
CutEdges = function(igraph_object, cluster_pairs = NULL, cluster.col = 'cluster') {
  ig = igraph_object

  # get data frame matching vertices with cluster ids
  vert.df = data.frame(vertices = igraph::as_ids(igraph::V(ig)),
                       cluster = igraph::vertex_attr(ig, cluster.col))

  # === update code to use igraph::crossing(), then identify edges to remove
  # from cluster_pairs argument 4/29/2025 ===

  # create edge data frame
  edge.df = igraph::as_ids(igraph::E(ig)) %>%
    stringr::str_split('\\|') %>%
    do.call(what = rbind) %>%
    as.data.frame()
  colnames(edge.df) = c('node1', 'node2')

  # match each node to their respective cluster id
  edge.df = edge.df %>%
    dplyr::mutate(node1_cl = vert.df$cluster[match(node1, vert.df$vertices)],
                  node2_cl = vert.df$cluster[match(node2, vert.df$vertices)])
  edge.df$edge_cl = apply(edge.df, 1, function(i) {
    # sort the cluster ids as characters; this is more reliable
    # than coercing to numeric since node ids could be characters.
    node_clusterid = as.character(c(i['node1_cl'], i['node2_cl']))
    stringr::str_flatten(sort(node_clusterid), collapse = '_')
  })

  # remove edges between cluster pairs
  cluster_pairs = lapply(cluster_pairs, function(x) {
    # coerce to character to match edge.df
    pair_of_clusters = as.character(x)
    sort(pair_of_clusters) %>% stringr::str_flatten(collapse = '_')
  }) %>% unlist
  edge.df = edge.df %>% dplyr::filter(edge_cl %in% cluster_pairs)
  ig = igraph::delete_edges(ig, paste(edge.df$node1, edge.df$node2, sep = '|'))
  return(ig)
}
