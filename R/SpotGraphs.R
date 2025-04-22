# Load libraries required for functions to work
library(igraph)

# Remove spots that are not part of the major tissue mass on each
# slide by identifying which spots are connected by:
# 1. Identify neighboring spots based on euclidean distance (on x/y coordinates)
# 2. Performing modularity maximization to identify clusters
# 3. Calculating total number of transcripts detected in each cluster
# 4. Setting a threshold to identify which clusters should be filtered out
CleanSlide = function(obj) {
  # Get coordinates, calculate euclidean distance, create igraph object
  coord = GetTissueCoordinates(obj)
  ig = SpotGraph(coord[,1:2], max.dist = 30)
  
  # Community detection via Modularity maximization
  mod_groups = cluster_fast_greedy(ig)
  ig = ig %>% set_vertex_attr('ig_cluster', index = V(ig), value = factor(mod_groups$membership))
  
  # Clustering assignments
  vertex_clusterid = data.frame(barcode = as_ids(V(ig)), ig_cluster = mod_groups$membership)
  meta = obj@meta.data %>% 
    mutate(barcode = rownames(.), ig_cluster = NULL) %>%
    left_join(vertex_clusterid, by = 'barcode') %>% 
    suppressMessages()
  
  # Calculate total nCount per cluster
  cl.df = meta %>% 
    reframe(.by = ig_cluster, totalcounts = sum(nCount_Spatial)) %>% 
    arrange(totalcounts)
  
  # Calculate nCount threshold
  # - log the total counts per cluster to smooth out the density
  # - take the mean of the top 5 clusters as the upper limit of the
  #   optimization window; this seems to work pretty well.
  den = density(log(cl.df$totalcounts))
  interval.max = cl.df$totalcounts %>% tail(5) %>% mean %>% log
  thres = exp(optimize(approxfun(den$x,den$y),interval=c(0,interval.max))$minimum)
  cl.df = cl.df %>% mutate(thres.pass = totalcounts > thres)
  meta = meta %>% left_join(cl.df) %>% suppressMessages()
  
  # Add results to Seurat object
  # - Return all results, let the user determine if the filter is
  #   appropriate for their analysis. If not, individual clusters can
  #   be manually re-added to "pass" the threshold
  obj = AddMetaData(obj, meta$thres.pass, 'threshold')
  obj = AddMetaData(obj, factor(meta$ig_cluster %>% ifelse(is.na(.), 'singlet', .)), 'ig_cluster')
  obj = AddMetaData(obj, meta$totalcounts, 'cluster_nCount')
  
  # Store igraph object in Seurat object
  obj@misc[['igraph_qc']] = ig
  
  return(obj)
}

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

# Create a function to plot on x,y coordinates some vertex
# attribute from a corresponding igraph object
# - Coordinate matrix should have two columns (x,y) with rownames
#   matching all of the vertices in the igraph object
# - function will assume first column is x coordinates and 
#   second column is y coordinates
# - group.by should be a vertex attribute in the igraph object
SpatialPlotGraph = function(igraph_object, coord, group.by, label = T) {
  ig = igraph_object
  v_names = as_ids(V(ig))
  if(!all(rownames(coord) %in% v_names)) {
    stop('igraph vertices do not match coordinates')
  }
  
  # Join coordinates together with vertex information
  colnames(coord) = c('x', 'y')
  coord = coord %>% mutate(barcode = rownames(.))
  df = data.frame(barcode = as_ids(V(ig)), groups = vertex_attr(ig, group.by))
  df = df %>% left_join(coord, by = 'barcode')
  
  # Create plot
  plt = ggplot(df, aes(x = y, y = -x)) +
    geom_point(aes(color = groups))
  
  # Label groups if desired
  if (label) {
    label.df = df %>% 
      reframe(.by = groups, x = mean(x), y = mean(y))
    plt = plt + 
      geom_label(data = label.df, 
                 aes(x = y, y = -x,
                     label = groups, fill = groups))
  }
  return(plt)
}

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