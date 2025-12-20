# SpatialPlotGraph produces expected output with coordinates

    Code
      waldo::compare(res_expected, res)
    Output
      `attr(old@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff1b8856250>
      `attr(new@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff198a89110>
      
      attr(old@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$myid vs attr(new@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$myid
      - "94977dcc-b409-40e7-94c0-6e4f82e526e2"
      + "d42c5609-fa3b-48e1-810b-ea81fcc5fb18"

# SpatialPlotGraph produces expected output with igraph

    Code
      waldo::compare(res_expected, res)
    Output
      `attr(old@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff1b8823d40>
      `attr(new@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff1a8915870>

