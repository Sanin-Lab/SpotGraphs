# Define a function to create an igraph object, given tissue coordinates
# - input should be a two-column dataframe or matrix
# - each column corresponds to x and y coordinates
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


