# SpatialPlotGraph produces expected output with coordinates input

    Code
      waldo::compare(res_expected@data, res@data, tolerance = testthat::testthat_tolerance())
    Output
      v No differences

---

    Code
      waldo::compare(names(res_expected@layers), names(res@layers), tolerance = testthat::testthat_tolerance())
    Output
      v No differences

