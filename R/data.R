#' Oral SCC Visium spatial transcriptomics data
#'
#' A subset of data from GSE208253
#'
#' @format Visium data from a fresh-frozen oral squamous cell carcinoma sample:
#' An object of class Seurat
#' 36601 features across 1185 samples within 1 assay
#' Active assay: Spatial (36601 features, 0 variable features)
#'  1 layer present: counts
#'  1 spatial field of view present: slice1
#' \describe{
#'   \item{Get x,y-coordinates of all spots}{Seurat::GetTissueCoordinates(scc_s1)}
#'   \item{Access Seurat object meta data}{scc_s1@meta.data}
#' }
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208253>
'scc_s1'
