test_that("min.dist is less than max.dist", {
  dist.buffer = 1.05
  min.dist = 10
  max.dist = sqrt(2*min.dist^2)*dist.buffer
  expect_gt(max.dist, min.dist)
})

test_that("coordinates are properly stored", {
  # expect_{}()
})

test_that("edge weights are expected values (.25, .5, 1)", {
  # expect_{}()
})
