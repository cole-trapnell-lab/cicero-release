#' Make an aggregate count cds by collapsing nearby peaks
#'
#' @param cds A CellDataSet (CDS) object. For example, output of
#' \code{\link{make_atac_cds}}
#' @param distance The distance within which peaks should be collapsed
#'
#' @return A CDS object with aggregated peaks.
#'
#' @details This function takes an input CDS object and collapses features
#'   within a given distance by summing the values for the collapsed features.
#'   Ranges of features are determined by their feature name, so the feature
#'   names must be in the form "chr1:1039013-2309023".
#'
#' @export
#'
#' @examples
#'   data("cicero_data")
#'   input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
#'   agg_cds <- aggregate_nearby_peaks(input_cds, distance = 10000)
#'
aggregate_nearby_peaks <- function(cds,
                                   distance = 1000) {
    assertthat::assert_that(assertthat::is.number(distance))
    assertthat::assert_that(is(cds, "CellDataSet"))

    fData(cds)$bin <- make_bin_col(cds, distance)
    cds <- cds[!is.na(fData(cds)$bin),]

    exprs_dt <- sparse_to_datatable(Matrix::Matrix(exprs(cds), sparse = TRUE))
    bin_info <- data.table::data.table(site = row.names(fData(cds)),
                                       bin = fData(cds)$bin)
    data.table::setkey(bin_info, "site")
    data.table::setkey(exprs_dt, "site")
    exprs_dt <- merge(exprs_dt, bin_info)

    data.table::setkey(exprs_dt, "cell", "bin")
    genomic_bins <- exprs_dt[,sum(val), by="cell,bin"]
    out <- Matrix::sparseMatrix(j=as.numeric(factor(genomic_bins$cell)),
                                i=as.numeric(factor(genomic_bins$bin)),
                                x=genomic_bins$V1)

    match_table <-
        data.table::data.table(num = as.numeric(factor(genomic_bins$bin)),
                               name = genomic_bins$bin)
    match_table <- unique(match_table)

    match_table2 <-
        data.table::data.table(num = as.numeric(factor(genomic_bins$cell)),
                               name = genomic_bins$cell)
    match_table2 <- unique(match_table2)

    fdf <- data.frame(site_name = levels(factor(genomic_bins$bin)),
                      row.names = levels(factor(genomic_bins$bin)))
    pdf <- data.frame(cells = levels(factor(genomic_bins$cell)),
                      row.names = levels(factor(genomic_bins$cell)))
    fdf$bin <- NULL
    pdf <- pdf[row.names(pData(cds)),]
    pdf <- cbind(pdf, pData(cds))
    pdf$pdf <- NULL

    data.table::setorder(match_table, "num")
    row.names(out) <- match_table$name

    data.table::setorder(match_table2, "num")
    colnames(out) <- match_table2$name

    out <- out[row.names(fdf), row.names(pdf)]

    fd <- new("AnnotatedDataFrame", data = fdf)
    pd <- new("AnnotatedDataFrame", data = pdf)

    if (is(exprs(cds), "dgCMatrix")) {
        compart_cds <-  suppressWarnings(newCellDataSet(as(out, "sparseMatrix"),
                                         phenoData = pd,
                                         featureData = fd,
                                         expressionFamily=negbinomial.size(),
                                         lowerDetectionLimit=0))
    } else {
        compart_cds <-  suppressWarnings(newCellDataSet(as.matrix(out),
                                         phenoData = pd,
                                         featureData = fd,
                                         expressionFamily=negbinomial.size(),
                                         lowerDetectionLimit=0))
    }

    return(compart_cds)
}


make_bin_col <- function(cds, distance) {
    coords_string_df <- df_for_coords(row.names(exprs(cds)))
    names(coords_string_df)[2:3] <- c("start", "stop")
    coords_ranges <- GenomicRanges::makeGRangesFromDataFrame(coords_string_df)
    coords_range_merge <- GenomicRanges::reduce(coords_ranges,
                                                min.gapwidth = distance)

    merge_df <- data.frame(seqnames=GenomicRanges::seqnames(coords_range_merge),
                           starts=GenomicRanges::start(coords_range_merge),
                           ends=GenomicRanges::end(coords_range_merge))
    merge_df$name <- paste(merge_df$seqnames,
                           merge_df$starts,
                           merge_df$ends, sep="_")

    overlaps <- GenomicRanges::findOverlaps(coords_ranges,
                                            coords_range_merge,
                                            select="first")
    overlaps <- as.data.frame(overlaps)

    merge_df <- merge_df[overlaps$overlaps,]
    merge_df$name
}

sparse_to_datatable <- function(sparse) {
    dgt_mat <- as(Matrix::t(sparse), "dgTMatrix")
    dt <- data.table(cell = dgt_mat@Dimnames[[1]][dgt_mat@i+1],
                     site=dgt_mat@Dimnames[[2]][dgt_mat@j+1],
                     val = dgt_mat@x)
    setkey(dt, "site", "cell")
    dt
}

#' Aggregate count CDS by groups of cells
#'
#' Aggregates a CDS based on an indicator column in the \code{pData} table
#'
#' @importFrom dplyr %>%
#' @importFrom plyr .
#' @param cds A CDS object to be aggregated
#' @param group_col The name of the column in the \code{pData} table that
#'   indicates the cells assignment to its aggregate bin.
#'
#' @details This function takes an input CDS object and collapses cells based
#'   on a column in the \code{pData} table by summing the values within the
#'   cell group.
#'
#' @return A count cds aggregated by group_col
#' @export
#'
#' @examples
#'   data("cicero_data")
#'   #input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
#'   #pData(input_cds)$cell_subtype <- rep(1:10, times=20)
#'   #binned_input_lin <-aggregate_by_cell_bin(input_cds, "cell_subtype")
#'
aggregate_by_cell_bin <- function(cds, group_col) {
    assertthat::assert_that(is(cds, "CellDataSet"))
    assertthat::assert_that(is.character(group_col))
    assertthat::assert_that(group_col %in% names(pData(cds)),
                            msg = "group_col is missing from your pData table")

    pData_grouping <- pData(cds) %>%
        tibble::rownames_to_column() %>%
        dplyr::group_by_at(group_col)

    cell_bins <- pData_grouping %>% dplyr::do(agg_cells(exprs(cds)[,.$rowname]))
    var_cols <- setdiff(colnames(cell_bins), c("site", "compartment_count"))

    agg_counts <- reshape2::dcast(cell_bins,
                                  as.formula(paste("site", "~",
                                      paste(var_cols, collapse="+"))),
                                  value.var="compartment_count")

    pData_cols <- as.data.frame(pData_grouping %>%
                                dplyr::group_by_at(group_col) %>%
                                dplyr::add_tally() %>%
                                dplyr::summarise_if(is.numeric,
                                                    mean,
                                                    na.rm = TRUE))

    rownames(pData_cols) <- colnames(agg_counts)[-1]

    fdf <- data.frame(site_name = agg_counts$site, row.names = agg_counts$site)

    bin_names <- colnames(agg_counts)[-1]

    pdf <- pData_cols

    fd <- new("AnnotatedDataFrame", data = fdf)
    pd <- new("AnnotatedDataFrame", data = pdf)
    out <- agg_counts[,bin_names]

    compart_cds <-  suppressWarnings(newCellDataSet(as.matrix(out),
                                     phenoData = pd,
                                     featureData = fd,
                                     expressionFamily=negbinomial.size(),
                                     lowerDetectionLimit=0))

    compart_cds <- detectGenes(compart_cds, min_expr=0.1)
    compart_cds <- estimateSizeFactorsSimp(compart_cds)
    compart_cds <- estimateDispersionsSimp(compart_cds)

    fData(compart_cds)$use_for_ordering <- FALSE

    compart_cds
}

agg_cells <- function(exprs_mat){
    cell_bins <- data.frame(compartment_count=Matrix::rowSums(exprs_mat))
    cell_bins$site <- row.names(exprs_mat)
    return (cell_bins)
}
