# Xenium data (BETA)

``` r

library(SpotGraphs)
library(dplyr)
library(ggplot2)
```

## Load example dataset

SpotGraph was designed for spot-based spatial transcripmotic data, but
we provide the function
[`CellGraph()`](https://potential-adventure-or7z9q9.pages.github.io/reference/CellGraph.md)
to use with cell coordinates generated with Xenium or other image based
technologies. Our goal with this package is support analysis by
providing a function that constructs an igraph object based on the
number of adjacent cells on a tissue section This opens up the analysis
pipeline to various igraph functions to perform certain tasks such as:

1.  identifying spots lying on tissue debris to be removed from
    downstream analysis
2.  distance calculations between spots by measuring the shortest path
    through the network between each spot
3.  easily identifying edges of pre-defined regions of tissue

To construct an igraph object, all we need are the x,y coordinates of
the spots on our slide. To demonstrate the general usage of our package
with this type of data, we collected coordinates from a Xenium Breast
Cancer dataset.

``` r

data("coord", package = "SpotGraphs")
class(coord)
#> [1] "data.frame"
dim(coord)
#> [1] 36913     2
head(coord)
#>                   x        y
#> bjfiibfn-1 681.4629 7513.613
#> bjfijchg-1 673.1639 7521.474
#> bjfikalj-1 678.1020 7518.308
#> bjfikgph-1 694.6487 7664.667
#> bjfimpba-1 704.1942 7543.554
#> bjfinfkh-1 702.9787 7535.547
ig = CellGraph(coord = coord)
ig
#> IGRAPH d4fb3bd U--- 36913 844453 -- 
#> + attr: coord_x (v/n), coord_y (v/n)
#> + edges from d4fb3bd:
#>  [1] 1--    2 1--    3 1--    5 1--    6 1--    7 1--    8 1--   10 1--   11
#>  [9] 1--   14 1--   15 1-- 6851 1-- 6854 1-- 6877 1-- 7023 1-- 7446 1-- 7447
#> [17] 1-- 7448 1-- 7449 1-- 7463 1-- 7466 1-- 7479 1-- 7480 1-- 7481 1-- 7485
#> [25] 1-- 7487 1-- 7488 1-- 7489 1-- 7490 1-- 7491 1-- 7856 1-- 7857 1-- 7858
#> [33] 1-- 7859 1-- 7860 1-- 7861 1-- 7862 1-- 7863 1-- 7864 1-- 7865 1-- 7866
#> [41] 1-- 7867 1-- 7871 1-- 7872 1-- 7873 1-- 7874 1-- 7875 1-- 7876 1-- 7877
#> [49] 1-- 7878 1-- 7880 1-- 7881 1-- 7883 1-- 7884 1-- 7885 1--36300 1--36301
#> [57] 2--    3 2--    5 2--    6 2--    7 2--    8 2--   10 2--   11 2--   14
#> + ... omitted several edges
```

We can plot this igraph object with the
[`SpatialPlotLargeGraph()`](https://potential-adventure-or7z9q9.pages.github.io/reference/SpatialPlotLargeGraph.md)
function provided in this package to visualize the edges drawn between
each node in the original x,y coordinates stored in the igraph object.
**As these technologies generate much larger graphs, it is often better
to visualize only the edges**, though the function can also plot the
nodes (colored by some attribute).

``` r

plt <- SpatialPlotLargeGraph(igraph_object = ig, only_edges = T, aspec.ratio = 1, group.by = NULL)
ggrastr::rasterise(plt, scale = 1, dpi = 300)
```

With this result, you may calculate any other property of a graph that
may be informative. For example, cell degree (number of connections for
each cell) varies greatly by cell type and may reveal interesting
biology.
