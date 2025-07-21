test_that("SpotGraph returns an igraph object", {
  df = rbind(expand.grid(1:5, 1:5), expand.grid(9:11, 9:11))
  colnames(df) = c('x', 'y')
  expect_equal(class(SpotGraph(df)), 'igraph')
})

test_that("density calculation returns a number within the range of nCount", {
  nCount = rbinom(20, 100, 0.5)
  thres = get_threshold(nCount)
  expect_true(thres > 0 & (thres-min(nCount) < max(nCount)-min(nCount)))
})

test_that("CleanSlide produces expected output", {
  coord = Seurat::GetTissueCoordinates(scc_s1)[,c('x', 'y')]
  nCount = scc_s1$nCount_Spatial
  res = CleanSlide(coord = coord, nCount = nCount)

  res_expected = readRDS(test_path('fixtures', 'cleanslide_expected.rds'))

  expect_snapshot(
    waldo::compare(res_expected, res)
  )
})
