#' Identify spots isolated from the largest groups of immediately adjacent spots
#' @description
#' Identify spots in a 10X Visium Seurat object that are potentially low
#' quality. Mostly useful for slides where there are large areas without tissue
#' and standard spot-filtering based on nCounts is insufficient.
#' 1. Identify neighboring spots on x,y coordinates with SpotGraph()
#' 2. Performing modularity maximization to identify clusters
#' 3. Calculate total number of transcripts detected in each cluster (nCount)
#' 4. Set a threshold to identify which clusters should be filtered out
#'
#' @param obj A Seurat object with 10X Visium data.
#'
#' @return A Seurat object with several columns of meta-data added:
#' 1. igraph modularity maximization results stored in 'ig_cluster'
#' 2. summed transcript counts within each 'ig_cluster' group in 'cluster_nCount'
#' 3. whether the respective 'ig_cluster' group passed the automatic threshold
#' @export
#'
#' @examples
#' library(Seurat)
#' obj = CleanSlide(obj)
#' SpatialDimPlot(obj, group.by = 'ig_cluster')
CleanSlide = function(obj) {
  # Get coordinates, calculate euclidean distance, create igraph object
  coord = GetTissueCoordinates(obj)
  ig = SpotGraph(coord[,1:2])

  # Community detection via Modularity maximization
  mod_groups = cluster_fast_greedy(ig)
  ig = set_vertex_attr(ig, 'ig_cluster', index = V(ig),
                       value = factor(mod_groups$membership))

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
  obj = AddMetaData(obj, factor(meta$ig_cluster %>% ifelse(is.na(.), 'singlet', .)), 'ig_cluster')
  obj = AddMetaData(obj, meta$totalcounts, 'cluster_nCount')
  obj = AddMetaData(obj, meta$thres.pass, 'threshold')

  # Store igraph object in Seurat object
  obj@misc[['igraph_qc']] = ig

  return(obj)
}
