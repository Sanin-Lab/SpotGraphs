test_that("CutEdges produces expected output", {
  ig = readRDS(test_path('fixtures', 'sccs1_igraph.rds'))
  clusters = readRDS(test_path('fixtures', 'cutedges_clusters.rds'))
  igraph::V(ig)$cluster = factor(clusters$membership)

  res = CutEdges(
    igraph_object = ig,
    cluster_pairs = list(
      c('1', '12'),
      c('1', '4')
    )
  )

  res_expected = readRDS(test_path('fixtures', 'cutedges_expected.rds'))

  expect_snapshot(
    waldo::compare(res_expected, res)
  )
})
