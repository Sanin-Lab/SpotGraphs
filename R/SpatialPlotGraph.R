#' Visualize igraph vertex attribute on tissue x,y coordinates
#'
#' @param igraph_object An igraph object
#' @param coord (optional) A two-column data.frame or matrix, where each column contains x or y coordinates.
#' if not provided, will look for `coord_x` and `coord_y` attributes in the provided igraph object.
#' @param group.by Some vertex attribute to plot onto tissue coordinates. Must be
#' present in `vertex_attr(ig)`.
#' @param label Logical, whether to label groups on the plot based on `group.by`.
#' @param flip.axes Default is TRUE. Invert the x and y axes, typically needed with coordinates from Visium Seurat objects to align the spots with the tissue image.
#' @param linewidth Numeric input to determine width of segments between nodes. Used in `geom_segment().`
#' @param pt.size Numeric input to determine size of each point on the plot. Used as the size parameter in `geom_point()`.
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
SpatialPlotGraph = function(igraph_object = NULL,
                            coord = NULL,
                            group.by = 'is_boundary',
                            label = FALSE,
                            flip.axes = TRUE,
                            linewidth = 0.5,
                            pt.size = 1.6) {
  # Create igraph object from coordinates if igraph object isn't provided
  if (is.null(igraph_object)) {
    ig = SpotGraph(coord)
  } else {
    ig = igraph_object
  }

  # Check if coordinates are provided, if not assume they are
  # stored in the igraph object from running SpotGraph()
  if (is.null(coord)) {
    coord = data.frame(
      coord_x = igraph::V(ig)$coord_x,
      coord_y = igraph::V(ig)$coord_y
      )
  } else {
    colnames(coord) = c('coord_x', 'coord_y')
    if(!all(rownames(coord) %in% names(V(ig)))) {
      stop('igraph vertices do not match coordinates')
    }
  }

  if (is.null(igraph_object) & is.null(coord)) {
    stop('must provide either igraph_object or coord')
  }

  if (flip.axes) {
    coord$x.new = coord$coord_y
    coord$y.new = -coord$coord_x
    coord[,c('coord_x', 'coord_y')] = NULL
    colnames(coord) = c('coord_x', 'coord_y')
  }

  # Create data.frame with all igraph vertex attributes and
  # coordinates for geom_segment using ggnetwork
  df = ggnetwork::ggnetwork(ig, layout = as.matrix(coord))

  # Create base plot
  plt = ggplot(df, aes(x = x, y = y)) +
    geom_segment(aes(xend = xend, yend = yend),
                 linewidth = linewidth) +
    theme_classic()

  # Add vertices based on whether group.by is present or not
  if (group.by %in% colnames(df)) {
    plt = plt +
      geom_point(aes_string(color = group.by), size = pt.size) +
      guides(color = guide_legend(title = group.by))
  } else {
    plt = plt +
      geom_point(size = pt.size)
  }

  # Label groups if desired
  if (label) {
    label.df = df %>%
      dplyr::reframe(.by = group.by, x = mean(x), y = mean(y))
    plt = plt +
      geom_label(data = label.df,
                 aes_string(x = 'x', y = 'y', label = group.by, fill = group.by),
                 show.legend = F)
  }
  return(plt)
}
