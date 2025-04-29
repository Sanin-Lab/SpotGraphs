#' Visualize igraph vertex attribute on tissue x,y coordinates
#'
#' @param igraph_object An igraph object
#' @param coord A two-column data.frame or matrix, where each column contains x or y coordinates.
#' @param group.by Some vertex attribute to plot onto tissue coordinates. Must be
#' present in vertex_attr(ig).
#' @param label Logical, whether to label groups on the plot based on group.by
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' # Create a coordinate data frame and create and igraph object
#' df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
#' colnames(df) = c('x', 'y')
#' ig = SpotGraph(df, cluster = TRUE)
#'
#' # Plot cluster results
#' SpatialPlotGraph(igraph_object = ig, coord = df, group.by = 'iglouvain_cluster')
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
                     label = groups, fill = groups)) + #get rid of all of that
  }
  return(plt)
}
