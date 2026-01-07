#' @keywords internal
"_PACKAGE"

utils::globalVariables(c('barcode',
                         'cluster_nCount',
                         'edge_cl',
                         'groups',
                         'ig_cluster',
                         'node1',
                         'node2',
                         'thres.pass',
                         'weight',
                         'x',
                         'y',
                         'xend',
                         'yend',
                         'clusters',
                         'name',
                         'center_eigen',
                         '.'))

## usethis namespace: start
#' @importFrom dplyr arrange
#' @importFrom dplyr filter
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr reframe
#' @importFrom dplyr select
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 annotate
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_vline
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 guide_colourbar
#' @importFrom ggplot2 guide_legend
#' @importFrom ggplot2 guides
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 theme_classic
#' @importFrom igraph as_ids
#' @importFrom igraph cluster_fast_greedy
#' @importFrom igraph cluster_louvain
#' @importFrom igraph components
#' @importFrom igraph degree
#' @importFrom igraph E
#' @importFrom igraph graph_from_edgelist
#' @importFrom igraph set_vertex_attr
#' @importFrom igraph V
#' @importFrom igraph vertex_attr
#' @importFrom interp tri.mesh
#' @importFrom rlang .data
#' @importFrom stats approxfun
#' @importFrom stats density
#' @importFrom stats dist
#' @importFrom stats optimize
#' @importFrom stats quantile
#' @importFrom stringr str_flatten
#' @importFrom stringr str_split
#' @importFrom tidyselect matches
#' @importFrom utils tail
## usethis namespace: end
NULL
