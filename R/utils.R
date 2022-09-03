#' Make ATAC CDS object
#'
#' This function takes as input a data frame or a path to a file in a sparse
#' matrix format and returns a properly formatted \code{CellDataSet} (CDS)
#' object.
#'
#' @param input Either a data frame or a path to input data. If a file, it
#'   should be a tab-delimited text file with three columns and no header. For
#'   either a file or a data frame, the first column is the peak coordinates in
#'   the form "chr10_100013372_100013596", the second column is the cell name,
#'   and the third column is an integer that represents the number of reads
#'   from that cell overlapping that peak. Zero values do not need to be
#'   included (sparse matrix format).
#' @param binarize Logical. Should the count matrix be converted to binary?
#'
#' @return A CDS object containing your ATAC data in proper format.
#' @export
#'
#' @examples
#'   data("cicero_data")
#'   input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
#'
make_atac_cds <- function(input, binarize = FALSE) {
  #check input:
  if(is(input, "character")) {
    assertthat::is.readable(input)
    intersect_lean <- as.data.frame(data.table::fread(input, header=FALSE))
  } else if (class(input) %in% c("matrix", "data.frame")) {
    intersect_lean <- input
  } else {
    stop("Input must be file path, matrix, or data.frame")
  }
  assertthat::assert_that(assertthat::are_equal(ncol(intersect_lean), 3))
  assertthat::assert_that(is.logical(binarize))
  
  names(intersect_lean) <- c("site_name", "cell_name", "read_count")
  
  assertthat::assert_that(is.numeric(intersect_lean$read_count))
  
  intersect_lean$site_name <- as.factor(intersect_lean$site_name)
  intersect_lean$cell_name <- as.factor(intersect_lean$cell_name)
  cellinfo <- data.frame(cells=levels(intersect_lean$cell_name))
  row.names(cellinfo) <- cellinfo$cells
  cellinfo$temp <- seq_len(nrow(cellinfo))
  
  dhsinfo <- data.frame(site_name = levels(intersect_lean$site_name))
  dhsinfo <- cbind(dhsinfo, split_peak_names(dhsinfo$site_name))
  row.names(dhsinfo) <- dhsinfo$site_name
  names(dhsinfo) <- c("site_name", "chr", "bp1", "bp2")
  dhsinfo$chr <- gsub("chr","", dhsinfo$chr)
  dhsinfo <- dhsinfo[order(as.character(dhsinfo$chr),
                           as.numeric(as.character(dhsinfo$bp2))),]
  
  intersect_lean_ord <- intersect_lean[order(intersect_lean$site_name,
                                             intersect_lean$cell_name),]
  dhsinfo <- dhsinfo[order(dhsinfo$site_name),]
  cellinfo <- cellinfo[order(cellinfo$cells),]
  intersect_lean_ord$site_name <- factor(intersect_lean_ord$site_name)
  intersect_lean_ord$cell_name <- factor(intersect_lean_ord$cell_name)
  intersect_lean_ord$site_name_num <- as.numeric(intersect_lean_ord$site_name)
  intersect_lean_ord$cell_name_num <- as.numeric(intersect_lean_ord$cell_name)
  
  if(binarize) intersect_lean_ord$read_count <-
    as.numeric(intersect_lean_ord$read_count > 0)
  sparse_intersect <- Matrix::sparseMatrix(i=intersect_lean_ord$site_name_num,
                                           j=intersect_lean_ord$cell_name_num,
                                           x=intersect_lean_ord$read_count)
  
  fd <- methods::new("AnnotatedDataFrame", data = dhsinfo)
  pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
  atac_cds <-  suppressWarnings(newCellDataSet(methods::as(sparse_intersect,
                                                           "sparseMatrix"),
                                               phenoData = pd,
                                               featureData = fd,
                                               expressionFamily=negbinomial.size(),
                                               lowerDetectionLimit=0))
  if(binarize) {
    atac_cds@expressionFamily <- binomialff()
    atac_cds@expressionFamily@vfamily <- "binomialff"
  }
  pData(atac_cds)$temp <- NULL
  fData(atac_cds)$chr <- as.character(fData(atac_cds)$chr)
  fData(atac_cds)$bp1 <- as.numeric(as.character(fData(atac_cds)$bp1))
  fData(atac_cds)$bp2 <- as.numeric(as.character(fData(atac_cds)$bp2))
  atac_cds <- atac_cds[order(fData(atac_cds)$chr, fData(atac_cds)$bp1),]
  atac_cds <- monocle::detectGenes(atac_cds)
  atac_cds
}

#' Construct GRanges objects from coordinate strings
#'
#' @param coord_strings A list of coordinate strings (in the form
#'   "chr1:500000-1000000")
#' @param with_names logical - should meta data include coordinate string
#'   (field coord_string)?
#' @param meta_data_df A data frame with any meta data columns you want
#'   included with the ranges. Must be in the same order as coord_strings.
#'
#' @details Coordinate strings consist of three pieces of information:
#'   chromosome, start, and stop. These pieces of information can be separated
#'   by the characters ":", "_", or "-". Commas will be removed, not used as
#'   separators (ex: "chr18:8,575,097-8,839,855" is ok).
#'
#' @return GRanges object of the input strings
#'
#' @examples
#'   ran1 <- ranges_for_coords("chr1:2039-30239", with_names = TRUE)
#'   ran2 <- ranges_for_coords(c("chr1:2049-203902", "chrX:489249-1389389"),
#'                             meta_data_df = data.frame(dat = c("1", "X")))
#'   ran3 <- ranges_for_coords(c("chr1:2049-203902", "chrX:489249-1389389"),
#'                             with_names = TRUE,
#'                             meta_data_df = data.frame(dat = c("1", "X"),
#'                                            stringsAsFactors = FALSE))
#'
#' @seealso \code{\link[GenomicRanges]{GRanges-class}}
#'
#' @export
ranges_for_coords <- function(coord_strings,
                              meta_data_df = NULL,
                              with_names = FALSE) {
  assertthat::assert_that(is.logical(with_names))
  if (!is.null(meta_data_df)) {
    assertthat::assert_that(is.data.frame(meta_data_df))
    assertthat::assert_that(assertthat::are_equal(length(coord_strings),
                                                  nrow(meta_data_df)))
  }
  
  coord_strings <- gsub(",", "", coord_strings)
  coord_cols <- split_peak_names(coord_strings)
  gr <- GenomicRanges::GRanges(coord_cols[, 1],
                               ranges = IRanges::IRanges(as.numeric(coord_cols[,2]),
                                                         as.numeric(coord_cols[, 3])),
                               mcols = meta_data_df)
  if (!is.null(meta_data_df)) {
    for (n in names(meta_data_df)) {
      newname <- paste0("mcols.", n)
      names(GenomicRanges::mcols(gr))[which(names(GenomicRanges::mcols(gr)) ==
                                              newname)] <- n
    }
  }
  if (with_names) {
    gr$coord_string <- coord_strings
  }
  gr
}

#' Construct a data frame of coordinate info from coordinate strings
#'
#' @param coord_strings A list of coordinate strings (each like
#'   "chr1:500000-1000000")
#'
#' @details Coordinate strings consist of three pieces of information:
#'   chromosome, start, and stop. These pieces of information can be separated
#'   by the characters ":", "_", or "-". Commas will be removed, not used as
#'   separators (ex: "chr18:8,575,097-8,839,855" is ok).
#'
#' @return data.frame with three columns, chromosome, starting base pair and
#'   ending base pair
#'
#' @examples
#'   df_for_coords(c("chr1:2,039-30,239", "chrX:28884:101293"))
#'
#' @export
df_for_coords <- function(coord_strings) {
  coord_strings <- gsub(",", "", coord_strings)
  coord_cols <- as.data.frame(split_peak_names(coord_strings),
                              stringsAsFactors = FALSE )
  names(coord_cols) <- c("chr", "bp1", "bp2")
  coord_cols$Peak <- coord_strings
  coord_cols$bp1 <- as.numeric(coord_cols$bp1)
  coord_cols$bp2 <- as.numeric(coord_cols$bp2)
  coord_cols
}

#' Add feature data columns to fData
#'
#' Annotate the sites of your CDS with feature data based on coordinate overlap.
#'
#' @param cds A CDS object.
#' @param feature_data Data frame, or a character path to a file of
#'   feature data. If a path, the file should be tab separated. Default assumes
#'   no header, if your file has a header, set \code{header = FALSE}. For
#'   either a data frame or a path, the file should be in bed-like format, with
#'   the first 3 columns containing chromosome, start and stop respectively.
#'   The remaining columns will be added to the \code{fData} table as feature
#'   data.
#' @param verbose Logical, should progress messages be printed?
#' @param maxgap The maximum number of base pairs allowed between the peak and
#'   the feature for the feature and peak to be considered overlapping.
#'   Default = 0 (overlapping). Details in
#'   \code{\link[IRanges]{findOverlaps-methods}}. If \code{maxgap}
#'   is set to "nearest" then the nearest feature will be assigned regardless
#'   of distance.
#' @param all Logical, should all overlapping intervals be reported? If all is
#'   FALSE, the largest overlap is reported.
#' @param header Logical, if reading a file, is there a header?
#'
#' @details \code{annotate_cds_by_site} will add columns to the \code{fData}
#'   table of a CDS object based on the overlap of peaks with features in a
#'   data frame or file. An "overlap" column will be added, along with any
#'   columns beyond the three required columns in the feature data. The
#'   "overlap" column is the number of base pairs overlapping the \code{fData}
#'   site. When maxgap is used, the true overlap is still calculated (overlap
#'   will be 0 if the two features only overlap because of maxgap) \code{NA}
#'   means that there was no overlapping feature. If a peak overlaps multiple
#'   data intervals and all is FALSE, the largest overlapping interval will be
#'   chosen (in a tie, the first entry is taken), otherwise all intervals will
#'   be chosen and annotations will be collapsed using a comma as a separator.
#'
#' @return A CDS object with updated \code{fData} table.
#'
#' @examples
#'   data("cicero_data")
#'   input_cds <- make_atac_cds(cicero_data, binarize = TRUE)
#'   feat <- data.frame(chr = c("chr18", "chr18", "chr18", "chr18"),
#'                      bp1 = c(10000, 10800, 50000, 100000),
#'                      bp2 = c(10700, 11000, 60000, 110000),
#'                      type = c("Acetylated", "Methylated", "Acetylated",
#'                      "Methylated"))
#'   input_cds <- annotate_cds_by_site(input_cds, feat)
#'
#' @importFrom IRanges findOverlaps
#'
#' @export
annotate_cds_by_site <- function(cds,
                                 feature_data,
                                 verbose = FALSE,
                                 maxgap = 0,
                                 all = FALSE,
                                 header = FALSE) {
  
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(is.character(feature_data) |
                            is.data.frame(feature_data))
  assertthat::assert_that(assertthat::is.number(maxgap) | maxgap == "nearest")
  assertthat::assert_that(is.logical(verbose))
  assertthat::assert_that(is.logical(all))
  assertthat::assert_that(is.logical(header))
  
  if (verbose) print("Generating fData ranges")
  granges <- ranges_for_coords(rownames(fData(cds)), with_names=TRUE)
  
  if (is(feature_data, "character")) {
    if (verbose) print("Reading data file")
    ch <- read.table(feature_data, header=header, stringsAsFactors = FALSE)
    if (verbose) print("Generating feature data ranges")
    names(ch)[c(1,2,3)] <- c("chr", "start", "stop")
    dtt <- GenomicRanges::makeGRangesFromDataFrame(ch,
                                                   keep.extra.columns = TRUE)
  } else {
    if (verbose) print("Generating feature data ranges")
    names(feature_data)[c(1,2,3)] <- c("chr", "start", "stop")
    dtt <- GenomicRanges::makeGRangesFromDataFrame(feature_data,
                                                   keep.extra.columns = TRUE)
  }
  
  if (verbose) print("Determining overlaps")
  
  if(maxgap == "nearest") {
    ol <- GenomicRanges::nearest(granges, dtt, select = "all")
  } else {
    ol <- GenomicRanges::findOverlaps(granges, dtt, select = "all",
                                      maxgap = maxgap)
  }
  
  olaps <- data.frame(
    row_name = GenomicRanges::mcols(granges[
      S4Vectors::queryHits(ol)])@listData$coord_string,
    width = GenomicRanges::width(
      IRanges::pintersect(granges[S4Vectors::queryHits(ol)],
                          dtt[S4Vectors::subjectHits(ol)]))
  )
  
  olaps <- cbind(olaps,
                 as.data.frame(GenomicRanges::mcols(
                   dtt[S4Vectors::subjectHits(ol)])))
  if (verbose) print("Assigning labels")
  
  if (all) {
    olaps <- olaps %>%
      dplyr::rename(overlap = width) %>%
      dplyr::group_by(row_name) %>%
      dplyr::summarise_all(paste, collapse=",")
  } else {
    olaps <- olaps[order(olaps$width, decreasing = TRUE),]
    olaps <- olaps[!duplicated(olaps$row_name),]
    olaps <- olaps %>% dplyr::rename(overlap = width)
  }
  
  if (verbose) print("Merging to fData table")
  fd <- fData(cds)
  fd$row_name <- row.names(fd)
  fd <- data.table::as.data.table(fd)
  data.table::setkey(fd, "row_name")
  
  fd <- merge(fd, olaps, by="row_name", all.x=TRUE)
  
  fData(cds) <- as.data.frame(fd)
  row.names(fData(cds)) <- fData(cds)$row_name
  fData(cds)$row_name <- NULL
  fData(cds) <- fData(cds)[row.names(exprs(cds)),]
  
  cds
}

#' Make a symmetric square sparse matrix from data frame
#'
#' Convert a data frame into a square sparse matrix (all versus all)
#'
#' @param data data frame
#' @param i.name name of i column
#' @param j.name name of j column
#' @param x.name name of value column
#'
#' @return sparse matrix
#'
#'
make_sparse_matrix <- function(data,
                               i.name = "Peak1",
                               j.name = "Peak2",
                               x.name = "value") {
  if(!i.name %in% names(data) |
     !j.name %in% names(data) |
     !x.name %in% names(data)) {
    stop('i.name, j.name, and x.name must be columns in data')
  }
  
  data$i <- as.character(data[,i.name])
  data$j <- as.character(data[,j.name])
  data$x <- data[,x.name]
  
  if(!class(data$x) %in%  c("numeric", "integer"))
    stop('x.name column must be numeric')
  
  peaks <- data.frame(Peak = unique(c(data$i, data$j)),
                      index = seq_len(length(unique(c(data$i, data$j)))))
  
  data <- data[,c("i", "j", "x")]
  
  data <- rbind(data, data.frame(i=peaks$Peak, j = peaks$Peak, x = 0))
  data <- data[!duplicated(data[,c("i", "j", "x")]),]
  data <- data.table::as.data.table(data)
  peaks <- data.table::as.data.table(peaks)
  data.table::setkey(data, "i")
  data.table::setkey(peaks, "Peak")
  data <- data[peaks]
  data.table::setkey(data, "j")
  data <- data[peaks]
  data <- as.data.frame(data)
  
  data <- data[,c("index", "i.index", "x")]
  data2 <- data
  names(data2) <- c("i.index", "index", "x")
  
  data <- rbind(data, data2)
  
  data <- data[!duplicated(data[,c("index", "i.index")]),]
  data <- data[data$index >= data$i.index,]
  
  sp_mat <- Matrix::sparseMatrix(i=as.numeric(data$index),
                                 j=as.numeric(data$i.index),
                                 x=data$x,
                                 symmetric = TRUE)
  
  colnames(sp_mat) <- peaks[order(peaks$index),]$Peak
  row.names(sp_mat) <- peaks[order(peaks$index),]$Peak
  return(sp_mat)
}

#' Compare Cicero connections to other datasets
#'
#' Compare two sets of connections and return a vector of logicals for whether
#' connections in one are present in the other.
#'
#' @param conns1 A data frame of Cicero connections, like those output from
#'   \code{assemble_connections}. The first two columns must be the coordinates
#'   of peaks that are connected.
#' @param conns2 A data frame of connections to be searched for overlap. The
#'   first two columns must be coordinates of genomic sites that are connected.
#' @param maxgap The number of base pairs between peaks allowed to be called
#'   overlapping. See \code{\link[IRanges]{findOverlaps-methods}} in the IRanges
#'   package for further description.
#'
#'
#' @return A vector of logicals of whether the Cicero pair is present in the
#'   alternate dataset.
#' @export
#'
#' @examples
#' \dontrun{
#' cons$in_dataset <- compare_connections(conns, alt_data)
#' }
#'
compare_connections <- function(conns1,
                                conns2,
                                maxgap = 0) {
  
  assertthat::assert_that(is(conns1, "data.frame"))
  assertthat::assert_that(is(conns2, "data.frame"))
  assertthat::assert_that(assertthat::is.number(maxgap))
  
  conns2 <- conns2[,c(1,2)]
  names(conns2) <- c("Peak1", "Peak2")
  conns22 <- conns2[,c(2,1)]
  names(conns22) <- c("Peak1", "Peak2")
  conns2 <- rbind(conns2, conns22)
  conns2 <- conns2[!duplicated(conns2),]
  
  alt1 <- ranges_for_coords(conns2$Peak1)
  alt2 <- ranges_for_coords(conns2$Peak2)
  
  conns11 <- ranges_for_coords(conns1$Peak1)
  conns12 <- ranges_for_coords(conns1$Peak2)
  
  
  ol1 <- GenomicRanges::findOverlaps(conns11,
                                     alt1,
                                     maxgap = maxgap,
                                     select="all")
  ol2 <- GenomicRanges::findOverlaps(conns12,
                                     alt2,
                                     maxgap = maxgap,
                                     select="all")
  ol1 <- as.list(ol1)
  ol2 <- as.list(ol2)
  
  olap <- mapply(function(o1, o2) {
    if (length(o1) == 0) out <- FALSE
    else if (length(intersect(o1,o2) > 0)) out <- TRUE
    else out <- FALSE
    return(out)
  }, ol1, ol2)
  
  return(olap)
}


#' Find peaks that overlap a specific genomic location
#'
#' @param coord_list A list of coordinates to be searched for overlap in the
#'   form chr_100_2000.
#' @param coord The coordinates that you want to find in the form chr1_100_2000.
#' @param maxgap The maximum distance in base pairs between coord and the
#'  coord_list that should count as overlapping. Default is 0.
#'
#' @return A character vector of the peaks that overlap coord.
#' @export
#'
#' @examples
#'   test_coords <- c("chr18_10025_10225", "chr18_10603_11103",
#'                    "chr18_11604_13986",
#'                    "chr18_157883_158536", "chr18_217477_218555",
#'                    "chr18_245734_246234")
#'   find_overlapping_coordinates(test_coords, "chr18:10,100-1246234")
#'
#'
find_overlapping_coordinates <- function(coord_list,
                                         coord,
                                         maxgap = 0) {
  coord <- gsub(",", "", coord)
  cons_gr <- ranges_for_coords(coord_list)
  if(length(coord) == 1) {
    ol1 <- GenomicRanges::findOverlaps(ranges_for_coords(coord),
                                       cons_gr,
                                       maxgap = maxgap,
                                       select="all")
    ol1 <- as.list(ol1)
    return(as.character(coord_list[unlist(ol1)]))
  } else {
    ol1 <- lapply(coord, function(x) {
      y <- suppressWarnings(unlist(as.list(
        GenomicRanges::findOverlaps(ranges_for_coords(x),
                                    cons_gr,
                                    maxgap = maxgap,
                                    select="all"))))
      if(length(y) == 0) return(NA)
      return(coord_list[y])
    })
    
    return(as.character(unlist(ol1)))
  }
}

is_chr <- function(x) {
  assertthat::assert_that(is.character(x))
  grepl("chr", x)
}

assertthat::on_failure(is_chr) <- function(call, env) {
  paste0(deparse(call$x), " must be of format 'chr1'")
}

is_color <- function(x, df=NULL) {
  if (!is.null(df)) {
    if (all(vapply(x, function(y) y %in% names(df), TRUE))) return(TRUE)
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
  } else {
    tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE)
  }
}

assertthat::on_failure(is_color) <- function(call, env) {
  paste0(deparse(call$x),
         " must be a valid color or a column in input data frame")
}

split_peak_names <- function(inp) {
  out <- stringr::str_split_fixed(stringi::stri_reverse(inp), 
                                  ":|-|_", 3)
  out[,1] <- stringi::stri_reverse(out[,1])
  out[,2] <- stringi::stri_reverse(out[,2])
  out[,3] <- stringi::stri_reverse(out[,3])
  out[,c(3,2,1), drop=FALSE]
}

########### Borrowed from monocle to avoid bug ##############

estimateSizeFactorsSimp <- function(cds) {
  if (any(class(exprs(cds)) %in% c("dgCMatrix", "dgTMatrix"))) {
    sizeFactors(cds) <- estimate_sf_sparse(exprs(cds))
  }else{
    sizeFactors(cds) <- estimate_sf_dense(exprs(cds))
  }
  return(cds)
}

estimate_sf_sparse <- function(counts){
  counts <- round(counts)
  cell_total <- Matrix::colSums(counts)
  sfs <- cell_total / exp(mean(log(cell_total)))
  sfs[is.na(sfs)] <- 1
  sfs
}

estimate_sf_dense <- function(counts){
  CM <- round(counts)
  cell_total <- apply(CM, 2, sum)
  sfs <- cell_total / exp(mean(log(cell_total)))
  sfs[is.na(sfs)] <- 1
  sfs
}

estimateDispersionsSimp <- function (object, modelFormulaStr = "~ 1", 
                                     relative_expr = TRUE, 
                                     min_cells_detected = 1, 
                                     remove_outliers = TRUE, cores = 1, 
                                     ...) {
  dispModelName = "blind"
  object@dispFitInfo = new.env(hash = TRUE)
  dfi <- estimateDispersionsForCellDataSet(object, modelFormulaStr, 
                                           relative_expr, min_cells_detected, 
                                           remove_outliers, 
                                           cores)
  object@dispFitInfo[[dispModelName]] <- dfi
  object
}


estimateDispersionsForCellDataSet <- function (cds, modelFormulaStr, 
                                               relative_expr, 
                                               min_cells_detected,
                                               removeOutliers, verbose = FALSE) 
{
  mu <- NA
  model_terms <- unlist(lapply(stringr::str_split(modelFormulaStr, "~|\\+|\\*"), 
                               stringr::str_trim))
  model_terms <- model_terms[model_terms != ""]
  cds_pdata <- dplyr::group_by(dplyr::select(tibble::rownames_to_column(pData(cds)), 
                                              "rowname"))
  disp_table <- as.data.frame(cds_pdata %>% 
                                dplyr::do(disp_calc_helper_NB(cds[, .$rowname],
                                                              cds@expressionFamily, 
                                                              min_cells_detected)))
  
  disp_table <- subset(disp_table, is.na(mu) == FALSE)
  res <- parametricDispersionFit(disp_table, verbose)
  fit <- res[[1]]
  coefs <- res[[2]]
  CD <- cooks.distance(fit)
  cooksCutoff <- 4/nrow(disp_table)
  outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), 
                                                         names(CD)))
  res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% 
                                              outliers == FALSE, ], verbose)
  fit <- res[[1]]
  coefs <- res[[2]]
  names(coefs) <- c("asymptDisp", "extraPois")
  ans <- function(q) coefs[1] + coefs[2]/q
  attr(ans, "coefficients") <- coefs
  res <- list(disp_table = disp_table, disp_func = ans)
  return(res)
}

disp_calc_helper_NB <- function (cds, expressionFamily, min_cells_detected) 
{
  rounded <- round(exprs(cds))
  nzGenes <- Matrix::rowSums(rounded > cds@lowerDetectionLimit)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])
  x <- t(t(rounded[nzGenes, ])/pData(cds[nzGenes, ])$Size_Factor)
  xim <- mean(1/pData(cds[nzGenes, ])$Size_Factor)
  f_expression_mean <- Matrix::rowMeans(x)
  f_expression_var <- Matrix::rowMeans((x - f_expression_mean)^2)
  disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
  disp_guess_meth_moments <- disp_guess_meth_moments/(f_expression_mean^2)
  res <- data.frame(mu = as.vector(f_expression_mean), 
                    disp = as.vector(disp_guess_meth_moments))
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA
  res$disp[res$disp < 0] <- 0
  res <- cbind(gene_id = row.names(fData(cds[nzGenes, ])), 
               res)
  res
}

parametricDispersionFit <- function (disp_table, verbose = FALSE, 
                                     initial_coefs = c(1e-06, 1)) {
  coefs <- initial_coefs
  iter <- 0
  while (TRUE) {
    residuals <- disp_table$disp/(coefs[1] + coefs[2]/disp_table$mu)
    good <- disp_table[which((residuals > initial_coefs[1]) & 
                               (residuals < 10000)), ]
    if (verbose) 
      fit <- glm(disp ~ I(1/mu), data = good, family = Gamma(link = "identity"), 
                 start = coefs)
    else suppressWarnings(fit <- glm(disp ~ I(1/mu), data = good, 
                                     family = Gamma(link = "identity"), 
                                     start = coefs))
    oldcoefs <- coefs
    coefs <- coefficients(fit)
    if (coefs[1] < initial_coefs[1]) {
      coefs[1] <- initial_coefs[1]
    }
    if (coefs[2] < 0) {
      stop("Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')")
    }
    if (sum(log(coefs/oldcoefs)^2) < initial_coefs[1]) 
      break
    iter <- iter + 1
    if (iter > 10) {
      warning("Dispersion fit did not converge.")
      break
    }
  }
  if (!all(coefs > 0)) {
    stop("Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')")
  }
  list(fit, coefs)
}







