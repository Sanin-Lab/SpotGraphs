test_that("CutEdges produces expected output", {
  ig = readRDS(test_path('fixtures', 'sccs1_igraph.rds'))
  clusters = readRDS(test_path('fixtures', 'cutedges_clusters.rds'))
  igraph::V(ig)$cluster = factor(clusters$membership)

  res = CutEdges(
    igraph_object = ig,
    cluster_pairs = list(
      c('5', '2'),
      c('5', '4')
    )
  )

  res_expected = readRDS(test_path('fixtures', 'cutedges_expected.rds'))

  # compare adjacency matrices
  res.adj = igraph::as_adjacency_matrix(res)
  res_expected.adj = igraph::as_adjacency_matrix(res_expected)

  expect_snapshot(
    waldo::compare(res_expected.adj, res.adj)
  )
})
