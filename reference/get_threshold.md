# Calculate a threshold given some set of statistics

Calculate a threshold given some set of statistics

## Usage

``` r
get_threshold(
  stats,
  percentile1 = 0.05,
  percentile2 = 0.95,
  log = T,
  print.plot = T
)
```

## Arguments

- stats:

  a numeric vector.

- percentile1:

  the percentile to set as the lower bound for calculating a local
  minimum.

- percentile2:

  the percentile to set as the upper bound for calculating a local
  minimum.

- log:

  logical (default `TRUE`), whether to use the log of the input
  statistics before calculating a threshold; useful for calculating
  threshold for filtering spots based on number of transcripts detected
  per-spot.

- print.plot:

  logical (default `TRUE`), whether to print a plot of approximated
  function with intervals (red) and threshold (blue) plotted.

## Value

The local minimum of the approximated function between the percentiles
set by `percentile1` and `percentile2`. A plot of the approximated
function will also be printed if `print.plot == TRUE` with the returned
value and intervals plotted.

## Examples

``` r
if (FALSE) { # \dontrun{
# To set a threshold for number of transcripts per spot to filter low quality spots:
seu_obj = scc_s1
nCount = seu_obj$nCount_Spatial
thres = get_threshold(nCount)
seu_obj_filtered = subset(seu_obj, nCount_Spatial > thres)
} # }
```
