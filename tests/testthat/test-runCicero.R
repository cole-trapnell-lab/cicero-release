context("runCicero")

#### make_cicero_cds ####

data(cicero_data)
load("../tsne_coord.Rdata")
data("human.hg19.genome")

sample_genome <- subset(human.hg19.genome, V1 == "chr18")
input_cds <- make_atac_cds(cicero_data)

set.seed(2017)
input_cds <- detectGenes(input_cds, min_expr = .1)
input_cds <- estimateSizeFactors(input_cds)
input_cds <- suppressWarnings(suppressMessages(estimateDispersions(input_cds)))

set.seed(2018)
cicero_cds <- make_cicero_cds(input_cds,
                              reduced_coordinates = tsne_coords,
                              silent = TRUE,
                              summary_stats = c("num_genes_expressed"))

test_that("make_cicero_cds aggregates correctly", {
  #skip_on_bioc()
  expect_is(cicero_cds, "CellDataSet")
  expect_equal(nrow(fData(cicero_cds)), nrow(fData(input_cds)))
  expect_named(pData(cicero_cds),c("agg_cell", "mean_num_genes_expressed",
                                   "Size_Factor", "num_genes_expressed"))
  expect_equal(nrow(exprs(cicero_cds)), 6146)
  expect_equal(ncol(exprs(cicero_cds)), 34)

  set.seed(2018)
  expect_warning(cicero_cds <- make_cicero_cds(input_cds,
                                 reduced_coordinates = tsne_coords,
                                 silent = FALSE,
                                 summary_stats = c("num_genes_expressed")))
  input_cds2 <- input_cds
  fData(input_cds2)$bp1 <- NULL
  set.seed(2018)
  cicero_cds <- make_cicero_cds(input_cds2,
                                reduced_coordinates = tsne_coords,
                                silent = TRUE,
                                size_factor_normalize = FALSE,
                                summary_stats = c("num_genes_expressed"))
  expect_is(cicero_cds, "CellDataSet")
  expect_equal(nrow(fData(cicero_cds)), nrow(fData(input_cds)))
  expect_named(pData(cicero_cds),c("agg_cell", "mean_num_genes_expressed",
                                   "Size_Factor", "num_genes_expressed"))
  expect_equal(nrow(exprs(cicero_cds)), 6146)
  expect_equal(ncol(exprs(cicero_cds)), 34)
})

set.seed(2018)
cicero_cds_temp <- make_cicero_cds(input_cds,
                                   reduced_coordinates = tsne_coords,
                                   silent = TRUE,
                                   summary_stats = c("num_genes_expressed"),
                                   return_agg_info = TRUE,
                                   size_factor_normalize = FALSE)
cicero_cds2 <- cicero_cds_temp[[1]]
agg_info <- cicero_cds_temp[[2]]

test_that("make_cicero_cds returns agg_info", {
  expect_is(cicero_cds2, "CellDataSet")
  expect_equal(nrow(fData(cicero_cds2)), nrow(fData(input_cds)))
  expect_named(pData(cicero_cds2),c("agg_cell", "mean_num_genes_expressed",
                                   "Size_Factor", "num_genes_expressed"))
  expect_equal(nrow(exprs(cicero_cds2)), 6146)
  expect_equal(ncol(exprs(cicero_cds2)), 34)
  
  expect_is(agg_info, "data.frame")

  agg_test_cell <- agg_info$agg_cell[[1]]
  test_agg <- as.character(agg_info[agg_info$agg_cell == agg_test_cell,]$cell)
  temp_exprs <- exprs(input_cds)[,test_agg]
  test_agg <- Matrix::rowSums(temp_exprs)
  expect_equal(sum(exprs(cicero_cds2)[,agg_test_cell] != test_agg), 0)
  
})

#### estimate_distance_parameter ####

test_that("estimate_distance_parameter gives correct mean", {
  #skip_on_bioc()
  set.seed(200)
  alphas <- estimate_distance_parameter(cicero_cds, window=500000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome,
                                        max_sample_windows = 6)
  mean_alpha <- mean(alphas)
  expect_equal(length(alphas), 2)
  expect_equal(mean_alpha, 2.25, tolerance = 1e-2)
  set.seed(200)
  alphas <- estimate_distance_parameter(cicero_cds, window=500000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = "../human.hg19.genome_sub.txt",
                                        max_sample_windows = 6)
  mean_alpha <- mean(alphas)
  expect_equal(length(alphas), 2)
  expect_equal(mean_alpha, 2.25, tolerance = 1e-2)
  set.seed(200)
  expect_error(expect_warning(alphas <- estimate_distance_parameter(cicero_cds,
                                        window=500000,
                                        maxit=100, sample_num = 2,
                                        max_elements = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome,
                                        max_sample_windows = 6)))
  testthat::expect_error(alphas <- estimate_distance_parameter(cicero_cds,
                                        window=10000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome,
                                        max_sample_windows = 6),
                         "distance_constraint not less than window")
  set.seed(205)
  testthat::expect_warning(alphas <- estimate_distance_parameter(cicero_cds,
                                        window=10000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 5000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome,
                                        max_sample_windows = 6))
})

#### generate_cicero_models ####
set.seed(203)
mean_alpha <- 2.030655
con_list <- generate_cicero_models(cicero_cds,
                                  mean_alpha,
                                  s=.75,
                                  genomic_coords = sample_genome)

test_that("generate_cicero_models gives output", { #slow
  skip_on_bioc()
  expect_is(con_list, "list")
  expect_equal(length(con_list), 313)
  expect_equal(con_list[[1]]$w[1,2], 0.866, tolerance = 1e-3)
  set.seed(203)
  con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     s=.75,
                                     genomic_coords = "../human.hg19.genome_sub.txt")
  expect_is(con_list, "list")
  expect_equal(length(con_list), 313)
  expect_equal(con_list[[1]]$w[1,2], 0.866, tolerance = 1e-3)
  set.seed(203)
  con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     window = 5000000,
                                     s=0.75,
                                     genomic_coords = sample_genome)
  expect_equal(length(con_list), 32)
  expect_equal(con_list[[1]], "Too many elements in range")

  set.seed(203)
  expect_error(con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     window = 5000000,
                                     s=1,
                                     genomic_coords = sample_genome),
               "s not less than 1")

  set.seed(203)
  con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     window = 500000,
                                     s=0.1,
                                     genomic_coords = sample_genome)
  expect_equal(con_list[[1]]$w[1,3], -3.7, tolerance = 1e-2)
})

#### assemble_connections ####

test_that("assemble_connections gives output", {
  #skip_on_bioc()
  expect_is(con_list, "list")
  expect_equal(length(con_list), 313)
  expect_equal(con_list[[1]]$w[1,2], 0.866, tolerance = 1e-3)
  set.seed(203)
  con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     s=.75,
                                     genomic_coords = "../human.hg19.genome_sub.txt")

  cons <- assemble_connections(con_list, silent = FALSE)
  expect_equal(cons[cons$Peak1 == "chr18_10025_10225" & 
                      cons$Peak2 == "chr18_10603_11103",]$coaccess, 0.877, 
               tolerance = 1e-3)
  expect_equal(ncol(cons), 3)
  expect_equal(nrow(cons), 543286)
})

#### run_cicero ####
set.seed(2000)
cons <- run_cicero(cicero_cds, sample_genome, window = 500000, silent=TRUE,
                   sample_num = 2)

test_that("run_cicero gives output", {
  #skip_on_bioc()
  expect_equal(cons[cons$Peak1 == "chr18_10025_10225" & 
                      cons$Peak2 == "chr18_10603_11103",]$coaccess, 0.877, 
               tolerance = 1e-3)
  expect_equal(ncol(cons), 3)
  expect_equal(nrow(cons), 543286)
  cons <- run_cicero(cicero_cds, window = 500000, silent=TRUE, sample_num = 2,
                     genomic_coords = "../human.hg19.genome_sub.txt")
  expect_equal(cons[cons$Peak1 == "chr18_10025_10225" & 
                      cons$Peak2 == "chr18_10603_11103",]$coaccess, 0.877, 
               tolerance = 1e-3)
  expect_equal(ncol(cons), 3)
  expect_equal(nrow(cons), 543286)
})

test_that("run_cicero gives output bad chromosomes", {
  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  input_cds <- make_atac_cds(cicero_data)
  
  fdata <- fData(input_cds)
  mtx <- exprs(input_cds)
  pdata <- pData(input_cds)
  row.names(fdata) <- gsub("chr", "A0", row.names(fdata))
  fdata$site_name <- row.names(fdata)
  row.names(mtx) <- row.names(fdata)
  pdata <- new("AnnotatedDataFrame", data = pdata)
  fdata <- new("AnnotatedDataFrame", data = fdata)
  new_inp <- suppressWarnings(newCellDataSet(mtx, pdata, fdata))
  
  set.seed(2017)
  new_inp <- detectGenes(new_inp, min_expr = .1)
  new_inp <- estimateSizeFactors(new_inp)
  
  set.seed(2018)
  cicero_cds <- make_cicero_cds(new_inp,
                                reduced_coordinates = tsne_coords,
                                silent = TRUE,
                                summary_stats = c("num_genes_expressed"))
  
  cons <- run_cicero(cicero_cds, window = 500000, silent=TRUE, sample_num = 2,
                     genomic_coords = "../human.hg19.genome_sub.txt")
  
  #skip_on_bioc()
  expect_equal(cons[cons$Peak1 == "A018_10025_10225" & 
                      cons$Peak2 == "A018_10603_11103",]$coaccess, 0.877, 
               tolerance = 1e-3)
  expect_equal(ncol(cons), 3)
  expect_equal(nrow(cons), 543286)
})

#### generate_ccans ####

test_that("generate_ccans gives output", { #slow
  skip_on_bioc()
  expect_output(CCAN_assigns <- generate_ccans(cons),
                "Coaccessibility cutoff used: 0.47")
  #expect_equal(CCAN_assigns["chr18_217477_218555",]$CCAN, 3, tolerance = 1e-7)
  expect_equal(ncol(CCAN_assigns), 2)
  expect_equal(nrow(CCAN_assigns), 1905)
  expect_equal(length(unique(CCAN_assigns$CCAN)), 116)

  expect_output(CCAN_assigns <- generate_ccans(cons,
                                               coaccess_cutoff_override = 0.25),
                "Coaccessibility cutoff used: 0.25")
  expect_output(CCAN_assigns <- generate_ccans(cons, tolerance_digits = 1),
                "Coaccessibility cutoff used: 0.5")
  expect_output(CCAN_assigns <- generate_ccans(cons, tolerance_digits = 1,
                                               coaccess_cutoff_override = .25),
                "Coaccessibility cutoff used: 0.25")
})

#### compare_connections ####

test_that("compare_connections works", {
  #skip_on_bioc()
  chia_conns <-  data.frame(Peak1 = c("chr18_10000_10200", "chr18_10000_10200",
                                      "chr18_49500_49600"),
                            Peak2 = c("chr18_10600_10700", "chr18_111700_111800",
                                      "chr18_10600_10700"))
  cons$in_dataset <- compare_connections(cons, chia_conns)
  cons$in_dataset2 <- compare_connections(cons, chia_conns, maxgap=1000)

  expect_is(cons, "data.frame")
  expect_equal(sum(cons$in_dataset), 4)
  expect_equal(sum(cons$in_dataset2), 22)
  expect_equal(cons[cons$Peak1 == "chr18_10025_10225" & 
                      cons$Peak2 == "chr18_10603_11103",]$in_dataset[1], TRUE)
})

#### find_overlapping_ccans ####

test_that("find_overlapping_ccans works", {
  CCAN_assigns <- generate_ccans(cons, coaccess_cutoff_override = 0.25)
  over <- find_overlapping_ccans(CCAN_assigns)
  expect_is(over, "data.frame")
  expect_equal(ncol(over), 2)
  skip_on_bioc()
  expect_equal(nrow(over), 98)
  over <- find_overlapping_ccans(CCAN_assigns, min_overlap = 3000000)
  expect_equal(nrow(over), 2)
})

#### activity scores ####

input_cds <- make_atac_cds(cicero_data, binarize=TRUE)
input_cds <- detectGenes(input_cds, min_expr = .1)

data(gene_annotation_sample)
gene_annotation_sub <- gene_annotation_sample[,c(1:3, 8)]
names(gene_annotation_sub)[4] <- "gene"

input_cds <- suppressWarnings(annotate_cds_by_site(input_cds, 
                                                   gene_annotation_sub))
unnorm_ga <- build_gene_activity_matrix(input_cds, cons)
expect_equal(nrow(unnorm_ga), 626)
expect_equal(ncol(unnorm_ga), 200)
expect_equal(unnorm_ga[1,1], 1.19, tolerance = 1e-2)

exprs(input_cds) <- as.matrix(exprs(input_cds))
unnorm_ga <- build_gene_activity_matrix(input_cds, cons)

test_that("build_gene_activity_matrix works", {
  #skip_on_bioc()
  expect_equal(nrow(unnorm_ga), 626)
  expect_equal(ncol(unnorm_ga), 200)
  expect_equal(unnorm_ga[1,1], 1.19, tolerance = 1e-2)
})

test_that("normalize_gene_activities works", {
  #skip_on_bioc()

  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))

  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  expect_equal(nrow(cicero_gene_activities), 626)
  expect_equal(ncol(cicero_gene_activities), 200)
  expect_equal(cicero_gene_activities[1,1], 0.0086, tolerance = 1e-4)

  cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga,
                                                           unnorm_ga),
                                                      num_genes)
  expect_is(cicero_gene_activities, "list")
  expect_equal(length(cicero_gene_activities), 2)
  cicero_gene_activities1 <- cicero_gene_activities[[1]]
  cicero_gene_activities2 <- cicero_gene_activities[[2]]
  expect_equal(nrow(cicero_gene_activities1), 626)
  expect_equal(ncol(cicero_gene_activities1), 200)
  expect_equal(cicero_gene_activities1[1,1], 0.0086, tolerance = 1e-4)

  expect_equal(nrow(cicero_gene_activities2), 626)
  expect_equal(ncol(cicero_gene_activities2), 200)
  expect_equal(cicero_gene_activities2[1,1], 0.0086, tolerance = 1e-4)


  unnorm_ga <- as.matrix(unnorm_ga)
  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  expect_equal(nrow(cicero_gene_activities), 626)
  expect_equal(ncol(cicero_gene_activities), 200)
  expect_equal(cicero_gene_activities[1,1], 0.0086, tolerance = 1e-4)

  cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga,
                                                           unnorm_ga),
                                                      num_genes)
  expect_is(cicero_gene_activities, "list")
  expect_equal(length(cicero_gene_activities), 2)
  cicero_gene_activities1 <- cicero_gene_activities[[1]]
  cicero_gene_activities2 <- cicero_gene_activities[[2]]
  expect_equal(nrow(cicero_gene_activities1), 626)
  expect_equal(ncol(cicero_gene_activities1), 200)
  expect_equal(cicero_gene_activities1[1,1], 0.0086, tolerance = 1e-4)

  expect_equal(nrow(cicero_gene_activities2), 626)
  expect_equal(ncol(cicero_gene_activities2), 200)
  expect_equal(cicero_gene_activities2[1,1], 0.0086, tolerance = 1e-4)

})

