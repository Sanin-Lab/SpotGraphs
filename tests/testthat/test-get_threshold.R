test_that("get_threshold produces expected output", {
  x = readRDS(test_path('fixtures', 'getthres_vector.rds'))
  res = get_threshold(x)
  res_expected = readRDS(test_path('fixtures', 'getthres_result.rds'))

  expect_snapshot(
    waldo::compare(res_expected, res)
  )
})
