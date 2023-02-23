#' cicero
#'
#' @import monocle
#' @import VGAM
#' @import data.table
#' @import ggplot2
#' @importFrom Biobase exprs pData fData ExpressionSet annotatedDataFrameFrom
#'   multiassign assayDataNew "fData<-" "pData<-"
#' @importFrom grDevices col2rgb dev.cur dev.off palette rainbow
#' @importFrom methods as callNextMethod is new
#' @importFrom stats as.formula cov dist filter median
#' @importFrom utils combn read.table
#' @importFrom BiocGenerics estimateDispersions estimateSizeFactors
"_PACKAGE"

## quiets concerns of R CMD check re: the .'s that appear in pipelines
#utils::globalVariables(c("."))

## temporary until i figure out a fix
#utils::globalVariables(c("val", "value", "CCAN", "V1", "f_id"))
