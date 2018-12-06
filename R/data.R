#' Example single-cell chromatin accessibility data
#'
#' A dataset containing a subset of a single-cell ATAC-seq
#' dataset collected on Human Skeletal Muscle Myoblasts.
#' Only includes data from chromosome 18.
#'
#' @format A data frame with 35137 rows and 3 variables:
#' \describe{
#'   \item{Peak}{Peak information}
#'   \item{Cell}{Cell ID}
#'   \item{Count}{Reads per cell per peak}
#' }
"cicero_data"

#' Chromosome lengths from human genome hg19
#'
#' A list of the chromosomes in hg19 and their lengths
#' in base pairs.
#'
#' @format A data frame with 93 rows and 2 variables:
#' \describe{
#'   \item{V1}{Chromosome}
#'   \item{V2}{Chromosome length, base pairs}
#' }
"human.hg19.genome"

#' Example gene annotation information
#'
#' Gencode gene annotation data from chromosome 18 of the
#' human genome (hg19).
#'
#' @format A data frame with 15129 rows and 8 variables:
#' \describe{
#'   \item{chromosome}{Chromosome}
#'   \item{start}{Exon starting base}
#'   \item{end}{Exon ending base}
#'   \item{strand}{Exon mapping direction}
#'   \item{feature}{Feature type}
#'   \item{gene}{Gene ID}
#'   \item{transcript}{Transcript ID}
#'   \item{symbol}{Gene symbol}
#' }
"gene_annotation_sample"

#' Metadata for example cells in cicero_data
#'
#' Metadata information for cicero_data
#'
#' @format A data frame with 200 rows and 2 variables:
#' \describe{
#'   \item{timepoint}{Time at cell collection}
#'   \item{cell}{Cell barcode}
#' }
"cell_data"
