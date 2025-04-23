# Remove all edges that lie between pairs of clusters
# - cluster_pairs should be a list containing two-element vectors
#   indicating pairs of cluster names where we want to remove edges
# - output is an updated igraph object with those edges removed
CutEdges = function(igraph_object, cluster_pairs = NULL, cluster.col = 'cluster') {
  ig = igraph_object

  # get data frame matching vertices with cluster ids
  vert.df = data.frame(vertices = as_ids(V(ig)),
                       cluster = vertex_attr(ig, cluster.col))

  # create edge data frame
  edge.df = data.frame(edges = as_ids(E(ig)))
  edge.df = as_ids(E(ig)) %>% strsplit('\\|') %>%
    do.call(what = rbind) %>% as.data.frame()
  colnames(edge.df) = c('node1', 'node2')

  # match each node to their respective cluster id
  edge.df = edge.df %>%
    mutate(node1_cl = vert.df$cluster[match(node1, vert.df$vertices)],
           node2_cl = vert.df$cluster[match(node2, vert.df$vertices)])
  edge.df$edge_cl = apply(edge.df, 1, function(i) {
    str_flatten(sort(c(i['node1_cl'], i['node2_cl'])), collapse = '_')
  })

  # remove edges between cluster pairs
  cluster_pairs = lapply(cluster_pairs, function(x) {
    sort(x) %>% str_flatten(collapse = '_')
  }) %>% unlist
  edge.df = edge.df %>% filter(edge_cl %in% cluster_pairs)
  ig = igraph::delete_edges(ig, paste(edge.df$node1, edge.df$node2, sep = '|'))
  return(ig)
}
