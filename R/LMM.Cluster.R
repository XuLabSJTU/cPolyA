###########
#' poly(A) sites clustering based on computational algorithm
#' @name LMM.Cluster
#' @usage LMM.Cluster(polyA, cutoff, method)
#' @param polyA A data.frame containing seqname, strand, poly(A) sites, and gene id.
#' @param cutoff cutoff from LMMcutoff function.
#' @param method polyA cluster methods including polyAclust,peaklu, QuantifyR and simple cluster.
#' @return The polyA data.frame will be return cluster of PolyA site.
#' @importFrom dplyr %>% mutate ungroup filter group_by arrange filter summarize summarise n tibble bind_rows
#' @import data.table
#' @importFrom IRanges IRanges extractList
#' @importFrom S4Vectors countQueryHits aggregate
#' @importFrom BiocGenerics width start
#' @import ggplot2
#' @import VGAM
#' @importFrom GenomicRanges GRanges findOverlaps follow distance reduce
#' @importFrom GenomeInfoDb seqlevels seqlevels<-
#' @import Rgb
#' @import pbmcapply
#' @import outliers
#' @import stringr
#' @importFrom stats dpois kmeans p.adjust density median quantile var
#' @importFrom utils globalVariables
utils::globalVariables(c("sum.wts","pCount","adPvalue","seqname","gene_id","coord","result","cluster","x","size","variance","tags","strand"))

#'@export
LMM.Cluster <- function(polyA, cutoff, method) {
  if (method == "polyAclust") {
    cluster.data <- polyAclust(polyA, cutoff)
  } else if (method == "peaklu") {
    tss <- polyA
    colnames(tss)[c(1,3,4)] <- c("chr","pos", "tags")
    cluster.data <- clusterByPeak(data.table(tss), peakDistance = 100, localThreshold = 0.02, extensionDistance = cutoff)
  } else if (method == "QuantifyR") {
    cluster.data <- Cluster.PolyA(polyA, max.gapwidth = cutoff, mc.cores = 1)
  } else {
    cluster.data <- simpleCluster(polyA, max.gapwidth = cutoff)
  }
  return(cluster.data)
}
