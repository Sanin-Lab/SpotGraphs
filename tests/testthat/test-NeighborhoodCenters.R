test_that("NeighborhoodCenters produces expected output", {
  coord = readRDS(test_path('fixtures', 'sccs1_coord.rds'))
  ig = SpotGraph(coord)
  igraph::V(ig)$name = rownames(coord)
  is_neighborhood = readRDS(test_path('fixtures', 'neighborhoodcenters_isneighborhood.rds'))

  res = NeighborhoodCenters(igraph_object = ig, is_neighborhood = is_neighborhood)
  res_expected = readRDS(test_path('fixtures', 'neighborhoodcenters_expected.rds'))
  expect_snapshot(
    waldo::compare(res_expected, res, tolerance = testthat::testthat_tolerance())
  )
})
