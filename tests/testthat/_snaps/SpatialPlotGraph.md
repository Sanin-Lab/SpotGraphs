# SpatialPlotGraph produces expected output with coordinates

    Code
      waldo::compare(res_expected, res)
    Output
      `attr(old@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff1bf2af140>
      `attr(new@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff1b8b293a0>
      
      attr(old@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$myid vs attr(new@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$myid
      - "94977dcc-b409-40e7-94c0-6e4f82e526e2"
      + "34c3c7ab-6913-4969-87c2-eda94732be9d"

# SpatialPlotGraph produces expected output with igraph

    Code
      waldo::compare(res_expected, res)
    Output
      `attr(old@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff1b899b610>
      `attr(new@layers$geom_segment$mapping$xend, '.Environment')$ig[[10]]$igraph` is <pointer: 0x7ff19c0edd50>

