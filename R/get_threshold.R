#' Calculate a threshold given some set of statistics
#'
#' @param stats a numeric vector.
#' @param percentile1 the percentile to set as the lower bound for calculating a local minimum.
#' @param percentile2 the percentile to set as the upper bound for calculating a local minimum.
#' @param log logical (default `TRUE`), whether to use the log of the input statistics before
#'  calculating a threshold; useful for calculating threshold for filtering spots based on
#'  number of transcripts detected per-spot.
#' @param print.plot logical (default `TRUE`), whether to print a plot of approximated function
#' with intervals (red) and threshold (blue) plotted.
#'
#' @returns The local minimum of the approximated function between the percentiles
#'  set by `percentile1` and `percentile2`. A plot of the approximated function will
#'  also be printed if `print.plot == TRUE` with the returned value and intervals plotted.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # To set a threshold for number of transcripts per spot to filter low quality spots:
#' seu_obj = scc_s1
#' nCount = seu_obj$nCount_Spatial
#' thres = get_threshold(nCount)
#' seu_obj_filtered = subset(seu_obj, nCount_Spatial > thres)
#' }
get_threshold = function(stats,
                         percentile1 = 0.05,
                         percentile2 = 0.95,
                         log = T,
                         print.plot = T) {
  # log the total counts per cluster to smooth out the density
  if (log) stats = log10(stats+1)
  den = density(stats)

  # define an interval to check for local minima
  int.min = spatstat.univar::quantile.density(den, percentile1)
  int.max = spatstat.univar::quantile.density(den, percentile2)

  interval = c(int.min, int.max)

  # calculate local minima
  thres = approxfun(den$x, den$y) %>% optimize(interval=interval)
  if (log) thres = 10^(thres$minimum) else thres = thres$minimum

  # print a plot of thresholding with stats and estimated function
  if (print.plot) {
    df = data.frame(x = den$x, y = den$y)

    plot.labels = c(
      paste('threshold:', round(log10(thres), 2)),
      paste('interval max:', round(int.max, 2)),
      paste('interval min:', round(int.min, 2))
    )

    plt = ggplot(df) +
      geom_line(aes(x = x, y = y)) +
      geom_vline(xintercept = interval, color = 'red') +
      geom_vline(xintercept = log10(thres), color = 'blue') +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      annotate(
        geom = 'text',
        x = 0.5,
        y = c(0.7, 0.75, 0.8),
        label = plot.labels
      )
    print(plt)
  }

  return(thres)
}


