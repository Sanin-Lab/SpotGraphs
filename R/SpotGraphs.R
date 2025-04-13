# Remove spots that are not part of the major tissue mass on each
# slide by identifying which spots are connected by:
# 1. Identify neighboring spots based on euclidean distance (on x/y coordinates)
# 2. Performing modularity maximization to identify clusters
# 3. Calculating total number of transcripts detected in each cluster
# 4. Setting a threshold to identify which clusters should be filtered out

CleanSlide = function(obj) {
  # Get coordinates and calculate euclidean distance
  coord = GetTissueCoordinates(obj)
  
  # Create igraph object with coordinates
  ig = SpotGraph(coord[,1:2], max.dist = 30)
  
  # Community detection via Modularity maximization
  mod_groups = cluster_fast_greedy(ig)
  ig = ig %>% set_vertex_attr('cluster', index = V(ig), value = factor(mod_groups$membership))
  
  # Clustering assignments
  vertex_clusterid = data.frame(barcode = as_ids(V(ig)), cluster = mod_groups$membership)
  meta = obj@meta.data %>% 
    mutate(barcode = rownames(.)) %>%
    left_join(vertex_clusterid) %>% 
    suppressMessages()
  
  # Calculate total nCount per cluster
  cl.df = meta %>% 
    reframe(.by = cluster, totalcounts = sum(nCount_Spatial)) %>% 
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
  obj = AddMetaData(obj, factor(meta$cluster %>% ifelse(is.na(.), 'singlet', .)), 'cluster')
  obj = AddMetaData(obj, meta$totalcounts, 'cluster_nCount')
  
  # Store igraph object in Seurat object
  obj@misc[['igraph_qc']] = ig
  
  return(obj)
}

# Define a function to create an igraph object, given tissue coordinates
# - input should be a two-column dataframe or matrix
# - each column corresponds to x and y coordinates

SpotGraph = function(coord, max.dist = 30) {
  # Get coordinates and calculate euclidean distance
  d = dist(coord, method = 'euclidean')
  m = as.matrix(d)
  
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
  return(ig)
}

