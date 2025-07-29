test_that("NeighborhoodCenters produces expected output", {
  coord = readRDS(test_path('fixtures', 'sccs1_coord.rds'))
  is_neighborhood = readRDS(test_path('fixtures', 'neighborhoodcenters_isneighborhood.rds'))

  res = NeighborhoodCenters(coord = coord, is_neighborhood = is_neighborhood)
  res_expected = readRDS(test_path('fixtures', 'neighborhoodcenters_expected.rds'))
  expect_snapshot(
    waldo::compare(res_expected, res)
  )
})
