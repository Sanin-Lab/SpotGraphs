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
