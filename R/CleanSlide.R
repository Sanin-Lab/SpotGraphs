#' Identify spots isolated from the largest groups of immediately adjacent spots
#'
#' @description Identify spots in a 10X Visium Seurat object that are potentially low quality. Mostly useful for slides where there are large areas without tissue and standard spot-filtering based on nCounts is insufficient.
#' 1. Identify neighboring spots on x,y coordinates with SpotGraph()
#' 2. Performing modularity maximization to identify clusters
#' 3. Calculate total number of transcripts detected in each cluster (nCount)
#' 4. Set a threshold to identify which clusters should be filtered out
#'
#' @param coord A two-column data.frame or matrix, where each column contains x or y coordinates.
#' @param nCount A named vector of transcript counts per spot, where names are
#' spot barcodes/ids, and each element should match each row of `coord`.
#'
#' @return A data.frame with the following columns:
#' 1. ig_cluster: igraph modularity maximization results
#' 2. cluster_nCount: summed transcript counts within each 'ig_cluster' group
#' 3. threshold: whether the respective 'ig_cluster' group passed the
#' automatically determiend threshold
#' 4. barcode: the spot id or barcode sequence
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a coordinate data frame with two isolated groups of points
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#' nspots = dim(df)[1]
#' rownames(df) = paste0('spot', 1:nspots)
#'
#' # Create some vector of transcript counts
#' spotcounts = rbinom(25, 10000, 0.3)
#' spotcounts = c(spotcounts, rbinom(9, 100, 0.3))
#' names(spotcounts) = rownames(df)
#'
#' # The resulting data.frame should indicate that the smaller
#' # cluster of 3x3 points did not pass the threshold, while the
#' # larger group of 5x5 points should have passed the threshold.
#' res = CleanSlide(coord = df, nCount = spotcounts)
#' }

CleanSlide = function(coord, nCount) {
  # Get coordinates, calculate euclidean distance, create igraph object
  ig = SpotGraph(coord)

  # Community detection via Modularity maximization
  cl = igraph::cluster_fast_greedy(ig)
  igraph::V(ig)$ig_cluster <- factor(cl$membership)

  # Clustering assignments
  vertex_clusterid = data.frame(barcode = names(V(ig)), ig_cluster = igraph::V(ig)$ig_cluster)
  meta = data.frame(barcode = names(nCount), nCount = nCount) %>%
    dplyr::left_join(vertex_clusterid, by = 'barcode')

  # Calculate total nCount per cluster
  cl.df = meta %>%
    dplyr::reframe(.by = ig_cluster, cluster_nCount = sum(nCount)) %>%
    dplyr::arrange(cluster_nCount)

  # Calculate nCount threshold
  thres = get_threshold(cl.df$cluster_nCount)
  cl.df = cl.df %>% dplyr::mutate(thres.pass = cluster_nCount > thres)
  meta = meta %>% dplyr::left_join(cl.df, by = 'ig_cluster')

  # Add results to Seurat object
  # - Return all results, let the user determine if the filter is
  #   appropriate for their analysis. If not, individual clusters can
  #   be manually re-added to "pass" the threshold
  meta = meta %>%
    dplyr::mutate(ig_cluster = factor(ifelse(is.na(ig_cluster), 'singlet', ig_cluster))) %>%
    dplyr::select(ig_cluster, cluster_nCount, threshold = thres.pass, barcode)

  rownames(meta) = meta$barcode

  return(meta)
}

#' Calculate a threshold given some set of statistics
#'
#' @param stats a numeric vector.
#' @param percentile1 the percentile to set as the lower bound for calculating a local minimum.
#' @param percentile2 the percentile to set as the upper bound for calculating a local minimum.
#' @param log logical (default `TRUE`), whether to use the log of the input statistics before
#'  calculating a threshold; useful for calculating threshold for filtering spots based on
#'  number of transcripts detected per-spot.
#' @param print.plot logical (default `TRUE`), whether to print a plot of approximated function
#' with intervals (red) and threshold (blue) plotted.
#'
#' @returns The local minimum of the approximated function between the percentiles
#'  set by `percentile1` and `percentile2`. A plot of the approximated function will
#'  also be printed if `print.plot == TRUE` with the returned value and intervals plotted.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # To set a threshold for number of transcripts per spot to filter low quality spots:
#' seu_obj = scc_s1
#' nCount = seu_obj$nCount_Spatial
#' thres = get_threshold(nCount)
#' seu_obj_filtered = subset(seu_obj, nCount_Spatial > thres)
#' }
get_threshold = function(stats,
                         percentile1 = 0.05,
                         percentile2 = 0.95,
                         log = T,
                         print.plot = T) {
  # log the total counts per cluster to smooth out the density
  if (log) stats = log10(stats+1)
  den = density(stats)

  # define an interval to check for local minima
  int.min = spatstat.univar::quantile.density(den, 0.05)
  int.max = spatstat.univar::quantile.density(den, 0.95)

  interval = c(int.min, int.max)

  # calculate local minima
  thres = approxfun(den$x, den$y) %>% optimize(interval=interval)
  if (log) thres = 10^(thres$minimum) else thres = thres$minimum

  # print a plot of thresholding with stats and estimated function
  if (print.plot) {
    df = data.frame(x = den$x, y = den$y)

    plot.labels = c(
      paste('threshold:', round(log10(thres), 2)),
      paste('interval max:', round(int.max, 2)),
      paste('interval min:', round(int.min, 2))
    )

    plt = ggplot(df) +
      geom_line(aes(x = x, y = y)) +
      geom_vline(xintercept = interval, color = 'red') +
      geom_vline(xintercept = log10(thres), color = 'blue') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      annotate(
        geom = 'text',
        x = 0.5,
        y = c(0.7, 0.75, 0.8),
        label = plot.labels
      )
    print(plt)
  }

  return(thres)
}


