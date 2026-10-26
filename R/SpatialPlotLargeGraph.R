#' Visualize a large graph on tissue x,y coordinates
#'
#' @param igraph_object An igraph object
#' @param group.by Some vertex attribute to plot onto tissue coordinates. Must be present in `vertex_attr(ig)`.
#' @param aspect.ratio Numeric input to constrain aspect ratio of the generated plot.
#' @param only_edges Boolean to determine if only the edges of the graph should be returned. Default is `TRUE`.
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a coordinate data frame and create and igraph object
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#' ig = SpotGraph(df)
#'
#' # Plot cluster results
#' SpatialPlotLargeGraph(igraph_object = ig)
#' }
SpatialPlotLargeGraph <- function(
    igraph_object = NULL,
    only_edges = T,
    aspec.ratio = 1.6,
    group.by = NULL){

  if (is.null(igraph_object)) {
    stop("Must provide a SpotGraphs igraph_object")
  }

  # Retrieved coordinates from the igraph object after running SpotGraph()
  coord = data.frame(coord_x = igraph::V(igraph_object)$coord_x, coord_y = igraph::V(igraph_object)$coord_y)

  # Create data.frame with all igraph vertex attributes and
  # coordinates for geom_segment using ggnetwork
  df = ggnetwork::ggnetwork(igraph_object, layout = as.matrix(coord))

  # Base plot
  plt = ggplot(df, aes(x = x, y = y)) +
    ggplot2::geom_segment(
      aes(xend = xend, yend = yend),
      linewidth = 0.05,
      alpha = 0.2) +
    ggplot2::theme_void() +
    ggplot2::theme(aspect.ratio = aspec.ratio)

  if(!only_edges){ # add vertex
    group.by = ifelse(is.null(group.by), NA, group.by)
    if (group.by %in% colnames(df)) {
      plt = plt +
        ggplot2::geom_point(
          aes(color = .data[[group.by]]),
          shape = ".", alpha = 0.5)

      if (is.numeric(df[, group.by])) {
        plt = plt + ggplot2::guides(color = ggplot2::guide_colourbar())
      }
    }
    else {
      plt = plt + ggplot2::geom_point(shape = ".", alpha = 0.5)
    }
  }

  # Return plot
  return(plt)
}
