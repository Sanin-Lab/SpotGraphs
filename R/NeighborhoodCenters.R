#' Identify center of some neighborhood of spots, e.g., the center of a tumor
#'
#' @description
#' Score spots using igraph::centr_eigen to identify spots in the center of some region/neighborhood.
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
#' \dontrun{
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#' nspots = dim(df)[1]
#' rownames(df) = paste0('spot', 1:nspots)
#'
#' igraph::V(ig)$cluster = igraph::cluster_leiden(ig, resolution = 0)$membership
#' is_cluster = igraph::V(ig)$cluster=='1'
#' names(is_cluster) = names(igraph::V(ig))
#'
#' res = NeighborhoodCenters(coord = df, is_neighborhood = is_cluster)
#' }
NeighborhoodCenters = function(coord = NULL, is_neighborhood) {
  # check if is_neighborhood vector is boolean
  if (!is.logical(is_neighborhood)) {
    stop('is_neighborhood should be logical (TRUE/FALSE) indicating neighborhood membership')
  }

  # Ensure spots in is_neighborhood and coord are aligned with each other
  # - if is_neighborhood is not named with spot ids (i.e., barcodes), assume
  #   that they are properly ordered in the same way as the coordinate matrix
  if (!is.null(names(is_neighborhood))) {
    coord = coord[match(rownames(coord), names(is_neighborhood)),]
  } else {
    names(is_neighborhood) = rownames(coord)
  }

  # Create igraph object
  ig = SpotGraph(coord)

  # Remove edges between neighborhood and non-neighborhood spots
  is_neighborhood = is_neighborhood[match(names(igraph::V(ig)), names(is_neighborhood))]
  igraph::V(ig)$is_neighborhood = is_neighborhood
  ig = CutEdges(igraph_object = ig,
                cluster_pairs = list(unique(is_neighborhood)),
                cluster.col = 'is_neighborhood')

  # Clustering to label each connected community
  cl = factor(igraph::components(ig)$membership)
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

  # Identify boundary vertices
  nbh_boundary = igraph::degree(ig) < max(igraph::degree(ig))
  boundary = names(igraph::V(ig))[nbh_boundary & is_neighborhood]

  # Identify center vertices
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

  # Format and return the output
  output = list(eigen.scores = finalscores,
                centers = centers,
                boundary = boundary)
  return(output)
}
