#' Visualize igraph vertex attribute on tissue x,y coordinates
#'
#' @param igraph_object An igraph object
#' @param coord (optional) A two-column data.frame or matrix, where each column contains x or y coordinates.
#' if not provided, will look for coord_x and coord_y attributes in the provided igraph object.
#' @param group.by Some vertex attribute to plot onto tissue coordinates. Must be
#' present in vertex_attr(ig).
#' @param label Logical, whether to label groups on the plot based on group.by
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
#' SpatialPlotGraph(igraph_object = ig, group.by = 'is_boundary')
#' }
SpatialPlotGraph = function(igraph_object, coord = NULL, group.by, label = T) {
  ig = igraph_object

  # Check if coordinates are provided, if not assume they are
  # stored in the igraph object from running SpotGraph()
  if (is.null(coord)) {
    coord = data.frame(coord_x = V(ig)$coord_x, coord_y = V(ig)$coord_y)
  } else {
    colnames(coord) = c('coord_x', 'coord_y')
    if(!all(rownames(coord) %in% names(V(ig)))) {
      stop('igraph vertices do not match coordinates')
    }
  }

  # Create data.frame with all igraph vertex attributes and
  # coordinates for geom_segment using ggnetwork
  df = ggnetwork::ggnetwork(ig, layout = as.matrix(coord))
  df$groups = df[,group.by]

  # Create plot
  plt = ggplot(df, aes(x = x, y = y)) +
    geom_segment(aes( xend = xend, yend = yend, alpha = weight)) +
    geom_point(aes(color = groups)) +
    guides(color = guide_legend(title = group.by))

  # Label groups if desired
  if (label) {
    label.df = df %>%
      reframe(.by = groups, x = mean(x), y = mean(y))
    plt = plt +
      geom_label(data = label.df,
                 aes(x = x, y = y, label = groups, fill = groups),
                 show.legend = F)
  }
  return(plt)
}
