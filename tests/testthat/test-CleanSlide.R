test_that("density calculation returns a number within the range of nCount", {
  nCount = rbinom(20, 100, 0.5)
  thres = get_threshold(nCount)
  expect_true(thres > 0 & (thres-min(nCount) < max(nCount)-min(nCount)))
})

test_that("CleanSlide produces expected output", {
  coord = readRDS(test_path('fixtures', 'sccs1_coord.rds'))
  nCount = readRDS(test_path('fixtures', 'sccs1_nCount.rds'))
  res = CleanSlide(coord = coord, nCount = nCount)

  res_expected = readRDS(test_path('fixtures', 'cleanslide_expected.rds'))

  expect_snapshot(
    waldo::compare(res_expected, res)
  )
})
