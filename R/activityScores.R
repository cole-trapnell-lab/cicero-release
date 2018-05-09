

#' Calculate initial Cicero gene activity matrix
#'
#' This function calculates the initial Cicero gene activity matrix. After this
#' function, the activity matrix should be normalized with any comparison
#' matrices using the function \code{\link{normalize_gene_activities}}.
#'
#' @param input_cds Binary sci-ATAC-seq input CDS. The input CDS must have a
#'   column in the fData table called "gene" which is the gene name if the site
#'   is a promoter, and \code{NA} if the site is distal.
#' @param cicero_cons_info Cicero connections table, generally the output of
#'   \code{\link{run_cicero}}. This table is a data frame with three required
#'   columns named "Peak1", "Peak2", and "coaccess". Peak1 and Peak2 contain
#'  coordinates for the two compared elements, and coaccess contains their
#'  Cicero co-accessiblity score.
#' @param site_weights NULL or an individual weight for each site in input_cds.
#' @param dist_thresh The maximum distance in base pairs between pairs of sites
#'   to included in the gene activity calculation.
#' @param coaccess_cutoff The minimum Cicero co-accessibility score that should
#'   be considered connected.
#'
#' @return Unnormalized gene activity matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' unnorm_ga <- build_gene_activity_matrix(input_cds, conns)
#' }
#'
build_gene_activity_matrix <- function(input_cds,
                                      cicero_cons_info,
                                      site_weights=NULL,
                                      dist_thresh=250000,
                                      coaccess_cutoff=0.25){
  assertthat::assert_that(class(input_cds) == "CellDataSet")
  assertthat::assert_that(is.data.frame(cicero_cons_info))
  assertthat::assert_that(assertthat::has_name(cicero_cons_info, "Peak1"),
                          assertthat::has_name(cicero_cons_info, "Peak2"),
                          assertthat::has_name(cicero_cons_info, "coaccess"))

  assertthat::assert_that(assertthat::has_name(fData(input_cds), "gene"),
                          msg = strwrap("The fData table of the input CDS must
                                         have a column called 'gene'. See
                                         documentation for details."))

  accessibility_mat <- exprs(input_cds)
  if (is.null(site_weights)){
    site_weights = Matrix::rowMeans(accessibility_mat) /
      Matrix::rowMeans(accessibility_mat)
    site_weights[names(site_weights)] = 1
  }

  gene_promoter_activity <-
    build_composite_gene_activity_matrix(input_cds,
                                         site_weights,
                                         cicero_cons_info,
                                         dist_thresh=dist_thresh,
                                         coaccess_cutoff=coaccess_cutoff)


  gene_activity_scores = gene_promoter_activity

  return(gene_activity_scores)
}

build_composite_gene_activity_matrix <- function(input_cds,
                                                 site_weights,
                                                 cicero_cons_info,
                                                 dist_thresh=250000,
                                                 coaccess_cutoff=0.25){
  accessibility_mat <- exprs(input_cds)
  promoter_peak_table <- fData(input_cds)
  promoter_peak_table$peak <- as.character(row.names(promoter_peak_table))
  promoter_peak_table <- promoter_peak_table[!is.na(promoter_peak_table$gene),]
  promoter_peak_table <- promoter_peak_table[,c("peak", "gene")]
  promoter_peak_table$gene = as.character(promoter_peak_table$gene)

  # Make site_weight matrix
  site_names = names(site_weights)
  site_weights = as(Matrix::Diagonal(x=as.numeric(site_weights)), "sparseMatrix")
  row.names(site_weights) = site_names
  colnames(site_weights) = site_names

  # Find distance between cicero peaks. If distance already calculated, skip
  if ("dist" %in% colnames(cicero_cons_info) == FALSE){
    Peak1_cols = stringr::str_split_fixed(cicero_cons_info$Peak1, "_", 3)
    Peak2_cols = stringr::str_split_fixed(cicero_cons_info$Peak2, "_", 3)
    Peak1_bp = round((as.integer(Peak1_cols[,3]) +
                        as.integer(Peak1_cols[,2])) / 2)
    Peak2_bp = round((as.integer(Peak2_cols[,3]) +
                        as.integer(Peak2_cols[,2])) / 2)
    cicero_cons_info$dist = abs(Peak2_bp - Peak1_bp)
  }

  # Get connections between promoters and distal sites above coaccess threshold
  nonneg_cons <- cicero_cons_info[(cicero_cons_info$Peak1 %in%
                                     promoter_peak_table$peak |
                                     cicero_cons_info$Peak2 %in%
                                     promoter_peak_table$peak) &
                                  cicero_cons_info$coaccess >= coaccess_cutoff &
                                  cicero_cons_info$dist < dist_thresh,]
  nonneg_cons <- nonneg_cons[,c("Peak1", "Peak2", "coaccess")]
  nonneg_cons <- nonneg_cons[!duplicated(nonneg_cons),]

  nonneg_cons$Peak1 = as.character(nonneg_cons$Peak1)
  nonneg_cons$Peak2 = as.character(nonneg_cons$Peak2)

  nonneg_cons = rbind(nonneg_cons,
                      data.frame(Peak1=unique(promoter_peak_table$peak),
                                 Peak2=unique(promoter_peak_table$peak),
                                 coaccess=0))

  # Make square matrix of connections from distal to proximal
  distal_connectivity_matrix = make_sparse_matrix(nonneg_cons,
                                                           x.name="coaccess")

  # Make connectivity matrix of promoters versus all
  promoter_conn_matrix =
    distal_connectivity_matrix[unique(promoter_peak_table$peak),]

  # Get list of promoter and distal sites in accessibility mat
  promoter_safe_sites <- intersect(rownames(promoter_conn_matrix),
                                   row.names(accessibility_mat))
  distal_safe_sites  = intersect(colnames(promoter_conn_matrix),
                                 row.names(accessibility_mat))
  distal_safe_sites = setdiff(distal_safe_sites, promoter_safe_sites)

  # Get accessibility info for promoters
  promoter_access_mat_in_cicero_map = accessibility_mat[promoter_safe_sites,]

  # Get accessibility for distal sites
  distal_activity_scores = accessibility_mat[distal_safe_sites,]

  # Scale connectivity matrix by site_weights
  scaled_site_weights = site_weights[distal_safe_sites,distal_safe_sites]
  total_linked_site_weights = promoter_conn_matrix[,distal_safe_sites] %*%
    scaled_site_weights
  total_linked_site_weights = 1/Matrix::rowSums(total_linked_site_weights,
                                                na.rm=T)
  total_linked_site_weights[is.finite(total_linked_site_weights) == FALSE] = 0
  total_linked_site_weights[is.na(total_linked_site_weights)] = 0
  total_linked_site_weights[is.nan(total_linked_site_weights)] = 0
  total_linked_site_weights = Diagonal(x=total_linked_site_weights)
  scaled_site_weights = total_linked_site_weights %*%
    promoter_conn_matrix[,distal_safe_sites] %*%
    scaled_site_weights
  scaled_site_weights@x[scaled_site_weights@x > 1] = 1

  # Multiply distal accessibility by site weights
  distal_activity_scores =  scaled_site_weights %*%
    distal_activity_scores

  distal_activity_scores =
    distal_activity_scores[row.names(promoter_access_mat_in_cicero_map),]

  # Sum distal and promoter scores
  promoter_activity_scores = distal_activity_scores +
    promoter_access_mat_in_cicero_map

  # Make and populate final matrix
  promoter_gene_mat <-
    Matrix::sparseMatrix(j=as.numeric(factor(promoter_peak_table$peak)),
                         i=as.numeric(factor(promoter_peak_table$gene)),
                         x=1)
  colnames(promoter_gene_mat) = levels(factor(promoter_peak_table$peak))
  row.names(promoter_gene_mat) = levels(factor(promoter_peak_table$gene))
  promoter_gene_mat = promoter_gene_mat[,row.names(promoter_activity_scores)]
  gene_activity_scores = promoter_gene_mat %*% promoter_activity_scores

  return(gene_activity_scores)
}



#' Normalize gene activities
#'
#' Normalize the output of \code{\link{build_gene_activity_matrix}}. Input is
#' either one or multiple gene activity matrices. Any gene activties to be
#' compared amongst each other should be normalized together.
#'
#'
#' @param activity_matrices A gene activity matix, output from
#'   \code{\link{build_gene_activity_matrix}}, or a list of gene activity
#'   matrices to be normalized together.
#' @param cell_num_genes A named vector of the total number of accessible sites
#'   per cell. Names should correspond to the cell names in the activity
#'   matrices. These values can be found in the "num_genes_expressed" column of
#'   the pData table of the CDS used to calculate the gene activity matrix.
#'
#' @return Normalized activity matrix or matrices.
#' @export
#'
#' @examples
#' \dontrun{
#' num_genes <- pData(input_cds)$num_genes_expressed
#' names(num_genes) <- row.names(pData(input_cds))
#' cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
#' }
#'
#'
normalize_gene_activities = function(activity_matrices,
                                     cell_num_genes){
  if (!is.list(activity_matrices)) {
    scores <- activity_matrices
    normalization_df <- data.frame(cell = colnames(activity_matrices),
                                            cell_group=1)
  } else {
    scores <- do.call(Matrix::cBind, activity_matrices)

    normalization_df <- do.call(rbind,
                               lapply(seq_along(activity_matrices), function(x) {
      data.frame(cell = colnames(activity_matrices[[x]]),
                 cell_group=rep(x, ncol(activity_matrices[[x]])))
    }))
  }

  normalization_df$cell_group <- factor(normalization_df$cell_group)
  normalization_df$total_activity <- Matrix::colSums(scores)
  normalization_df$total_sites <-
    cell_num_genes[as.character(normalization_df$cell)]

  if (!is.list(activity_matrices)) {
    activity_model <- stats::lm(log(total_activity) ~ log(total_sites),
                         data=normalization_df)
  } else {
    activity_model <- stats::lm(log(total_activity) ~ log(total_sites) * cell_group,
                         data=normalization_df)
  }

  normalization_df$fitted_curve <- exp(as.vector(predict(activity_model,
                                                         type="response")))

  size_factors <- log(normalization_df$fitted_curve) /
    mean(log(normalization_df$fitted_curve))

  size_factors <- Matrix::Diagonal(x=1/size_factors)
  row.names(size_factors) <- normalization_df$cell
  colnames(size_factors) <- row.names(size_factors)

  # Adjust the scores by the size factors
  scores <- Matrix::t(size_factors %*% Matrix::t(scores))

  scores@x <- pmin(1e9, exp(scores@x) - 1)

  sum_activity_scores <- Matrix::colSums(scores)

  scale_factors <- Diagonal(x=1/sum_activity_scores)
  row.names(scale_factors) <- normalization_df$cell
  colnames(scale_factors) <- row.names(scale_factors)

  scores <- Matrix::t(scale_factors %*% Matrix::t(scores))

  if (!is.list(activity_matrices)) {
    ret <- scores[row.names(activity_matrices), colnames(activity_matrices)]
  } else {
    ret <- lapply(activity_matrices, function(x) {
      scores[row.names(x), colnames(x)]
    })
  }
  return(ret)
}


