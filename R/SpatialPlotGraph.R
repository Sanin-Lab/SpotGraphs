# Create a function to plot on x,y coordinates some vertex
# attribute from a corresponding igraph object
# - Coordinate matrix should have two columns (x,y) with rownames
#   matching all of the vertices in the igraph object
# - function will assume first column is x coordinates and
#   second column is y coordinates
# - group.by should be a vertex attribute in the igraph object
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
                     label = groups, fill = groups))
  }
  return(plt)
}
