test_that("SpotGraph returns an igraph object", {
  df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
  colnames(df) = c('x', 'y')
  expect_equal(class(SpotGraph(df)), 'igraph')
})

test_that("density calculation results in a number in the range of nCount", {
  nCount = rbinom(20, 100, 0.5)
  thres = get_threshold(nCount)
  expect_true(thres > 0 & thres < max(nCount-min(nCount)))
})
