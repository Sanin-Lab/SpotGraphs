# check that data.frame stored in each ggplot object are the same
# check that the layer names are the same

test_that("SpatialPlotGraph produces expected output with coordinates input", {
  coord = readRDS(test_path('fixtures', 'sccs1_coord.rds'))
  # res = SpatialPlotGraph(coord = coord)
  res = readRDS(test_path('fixtures', 'spatialplotgraph_coord.rds'))
  res_expected = readRDS(test_path('fixtures', 'spatialplotgraph_coord.rds'))

  expect_snapshot(
    waldo::compare(
      res_expected@data,
      res@data,
      tolerance = testthat::testthat_tolerance()
    )
  )
  expect_snapshot(
    waldo::compare(
      names(res_expected@layers),
      names(res@layers),
      tolerance = testthat::testthat_tolerance()
    )
  )
})

# test_that("SpatialPlotGraph produces expected output with igraph input", {
#   ig = readRDS(test_path('fixtures', 'sccs1_igraph.rds'))
#   res = SpatialPlotGraph(igraph_object = ig)
#   res_expected = readRDS(test_path('fixtures', 'spatialplotgraph_igraph.rds'))
#
#   expect_snapshot(
#     waldo::compare(
#       res_expected@data,
#       res@data,
#       tolerance = testthat::testthat_tolerance()
#     )
#   )
#   expect_snapshot(
#     waldo::compare(
#       names(res_expected@layers),
#       names(res@layers),
#       tolerance = testthat::testthat_tolerance()
#     )
#   )
# })

