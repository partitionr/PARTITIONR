#' PARTITIONR: A package for partitioning diversity across spatial scales
#'
#' The PARTITIONR package provides tools for partitioning diversity across spatial scales. It has
#'     three types of functions.
#'
#' @section Partition Function:
#' \code{\link{partition}} allows for both individual- and sample-based randomizations
#'     of the provided data. Individual-based randomizations can be run with data
#'     from single level and non-nested experimental designs. Sample-based randomizations
#'     can ONLY be run with data from nested (hierarchical) experimental designs.
#'
#' @section Support Functions:
#' The support functions provide an opportunity to visualize the data and output of
#'     the \code{\link{partition}} function. \code{\link{partiplot}} visualizes the
#'     output of \code{\link{partition}}. \code{\link{summary.partition}} provides
#'     observed diversity, expected diversity, and associated p-values (based on
#'     randomizations). 
#'
#' @section Community Analysis Functions:
#' The community analysis functions are NOT functional yet. Working through these.
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
