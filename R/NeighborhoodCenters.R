#' Identify center of some neighborhood of spots, e.g., the center of a tumor
#'
#' @param coord A two-column data.frame or matrix, where each column contains x or y coordinates,
#' where the rownames are spot barcodes/ids.
#' @param is_neighborhood A boolean vector (`TRUE`/`FALSE`) of equal length as the number of rows
#' in `coord`. Centers will be identified for each neighborhood of spots labeled `TRUE`.
#'
#' @returns a list with the following elements:
#' 1. `eigen.scores`: a numeric vector with the output from igraph::centr_eigen, with a
#' score for each spot, where `max(eigen.scores)` within each neighborhood will
#' identify its center.
#' 2. `centers`: a character vector with barcode ids of spots that were identified as the
#' center of each neighborhood of spots.
#' 3. `boundaries`: a character vector with barcode ids of spots that were identified be on
#' the boundary of each neighborhood of spots.
#' @export
#'
#' @examples
#' \dontrun{}
NeighborhoodCenters = function(coord = NULL, is_neighborhood) {

  # make sure barcodes in coordinates and label vector are aligned
  coord = coord[match(rownames(coord), names(is_neighborhood)),]
  ig = SpotGraph(coord)

  # Remove edges between tumor and non-tumor spots
  is_neighborhood = is_neighborhood[match(names(igraph::V(ig)), names(is_neighborhood))]
  igraph::V(ig)$is_neighborhood = is_neighborhood
  ig = CutEdges(igraph_object = ig,
                cluster_pairs = list(unique(is_neighborhood)),
                cluster.col = 'is_neighborhood')

  # Clustering to label each connected community
  cl = factor(igraph::cluster_leiden(ig, resolution = 0)$membership)
  igraph::V(ig)$clusters = cl

  # Calculate the centralization eigenvector for each cluster
  # cl = igraph::V(ig)$clusters
  v.scores = lapply(levels(cl), function(xx) {
    ig.sub = igraph::subgraph(ig, which(cl==xx))
    res = igraph::centr_eigen(ig.sub)$vector
    names(res) = names(igraph::V(ig.sub))
    return(res)
  }) %>% unlist
  v.scores = v.scores[match(names(igraph::V(ig)), names(v.scores))]
  igraph::V(ig)$center_eigen = v.scores

  # Identify center vertices of tumor clusters
  centers = as.data.frame(igraph::vertex_attr(ig)) %>%
    filter(.by = clusters,
           is_neighborhood &
             center_eigen == max(center_eigen)) %>%
    .$name

  # Set eigenvalue to 0 on spots where is_neighborhood == F
  finalscores = as.data.frame(igraph::vertex_attr(ig)) %>%
    mutate(center_eigen = ifelse(is_neighborhood, center_eigen, 0)) %>%
    select(barcode = name, center_eigen)

  # Ensure that the output is in the same order as the input
  finalscores = finalscores[match(names(is_neighborhood), finalscores$barcode),]

  # output = list(eigen.scores = finalscores,
  #               centers = centers,
  #               boundaries = boundaries)
  output = list(eigen.scores = finalscores,
                centers = centers)
  return(output)
}
