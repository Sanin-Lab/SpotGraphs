test_that("SpotGraph returns an igraph object", {
  df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
  colnames(df) = c('x', 'y')
  expect_equal(class(SpotGraph(df)), 'igraph')
})

test_that("min.dist is less than max.dist", {
  dist.buffer = 1.05
  min.dist = 10
  max.dist = sqrt(2*min.dist^2)*dist.buffer
  expect_gt(max.dist, min.dist)
})

test_that("SpotGraph produces expected output", {
  coord = readRDS(test_path('fixtures', 'sccs1_coord.rds'))

  res = SpotGraph(coord = coord)
  res_expected = readRDS(test_path('fixtures', 'sccs1_igraph.rds'))

  expect_snapshot(
    waldo::compare(res_expected, res)
  )
})
