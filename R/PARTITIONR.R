#' PARTITIONR: A package for partitioning diversity across spatial scales
#'
#' The PARTITIONR package provides tools for partitioning diversity across spatial scales. It has
#'     three types of functions.
#'
#' @section Partition Function:
#' The partition function partitions diversity across nested scales.
#'
#' @section Support Functions:
#' The support functions provide an opportunity to visualize the data and output of
#'     the \code{\link{partition}} function. These functions are good.
#'
#' @section Community Analysis Functions:
#' The community analysis functions are NOT FUNCTIONAL YET.
#'
#' @importFrom plyr rbind.fill
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise_all
#' @importFrom dplyr count_
#' @importFrom dplyr sample_n
#' @importFrom dplyr distinct_
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom parallel clusterExport
#' @importFrom parallel stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach
#' @import ggplot2
#'
#' @docType package
#' @name PARTITIONR
NULL
