context("runCicero")

#### make_cicero_cds ####

test_that("make_cicero_cds aggregates correctly", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
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
  expect_is(cicero_cds, "CellDataSet")
  expect_equal(nrow(fData(cicero_cds)), nrow(fData(input_cds)))
  expect_named(pData(cicero_cds),c("agg_cell", "mean_num_genes_expressed",
                                   "Size_Factor", "num_genes_expressed"))
  expect_equal(nrow(exprs(cicero_cds)), 6146)
  expect_equal(ncol(exprs(cicero_cds)), 36)

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
  expect_equal(ncol(exprs(cicero_cds)), 36)
})

#### estimate_distance_parameter ####

test_that("estimate_distance_parameter gives correct mean", {
  skip_on_bioc()
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
  set.seed(200)
  alphas <- estimate_distance_parameter(cicero_cds, window=500000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome)
  mean_alpha <- mean(alphas)
  expect_equal(length(alphas), 2)
  expect_equal(mean_alpha, 2.030655, tolerance = 1e-6)
  set.seed(200)
  alphas <- estimate_distance_parameter(cicero_cds, window=500000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = "../human.hg19.genome_sub.txt")
  mean_alpha <- mean(alphas)
  expect_equal(length(alphas), 2)
  expect_equal(mean_alpha, 2.030655, tolerance = 1e-6)
  set.seed(200)
  expect_error(expect_warning(alphas <- estimate_distance_parameter(cicero_cds,
                                        window=500000,
                                        maxit=100, sample_num = 2,
                                        max_elements = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome)))
  testthat::expect_error(alphas <- estimate_distance_parameter(cicero_cds,
                                        window=10000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 250000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome),
                         "distance_constraint not less than window")
  set.seed(202)
  testthat::expect_warning(alphas <- estimate_distance_parameter(cicero_cds,
                                        window=10000,
                                        maxit=100, sample_num = 2,
                                        distance_constraint = 5000,
                                        distance_parameter_convergence = 1e-22,
                                        genomic_coords = sample_genome))
})


#### generate_cicero_models ####


test_that("generate_cicero_models gives output", {
  skip_on_bioc()
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
  set.seed(203)
  mean_alpha <- 2.030655
  con_list <- generate_cicero_models(cicero_cds,
                                    mean_alpha,
                                    s=.75,
                                    genomic_coords = sample_genome)
  expect_is(con_list, "list")
  expect_equal(length(con_list), 313)
  expect_equal(con_list[[1]]$w[1,2], 1.028059, tolerance = 1e-6)
  set.seed(203)
  con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     s=.75,
                                     genomic_coords = "../human.hg19.genome_sub.txt")
  expect_is(con_list, "list")
  expect_equal(length(con_list), 313)
  expect_equal(con_list[[1]]$w[1,2], 1.028059, tolerance = 1e-6)
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
  expect_equal(con_list[[1]]$w[1,3], -4.225899, tolerance = 1e-5)
})


#### assemble_connections ####



test_that("assemble_connections gives output", {
  skip_on_bioc()
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
  set.seed(203)
  mean_alpha <- 2.030655
  con_list <- generate_cicero_models(cicero_cds,
                                    mean_alpha,
                                    s=.75,
                                    genomic_coords = sample_genome)
  expect_is(con_list, "list")
  expect_equal(length(con_list), 313)
  expect_equal(con_list[[1]]$w[1,2], 1.028059, tolerance = 1e-6)
  set.seed(203)
  con_list <- generate_cicero_models(cicero_cds,
                                     mean_alpha,
                                     s=.75,
                                     genomic_coords = "../human.hg19.genome_sub.txt")

  cons <- assemble_connections(con_list, silent = FALSE)
  expect_equal(cons$coaccess[1], 0.908402, tolerance = 1e-7)
  expect_equal(ncol(cons), 3)
  expect_equal(nrow(cons), 543482)
})

#### run_cicero ####

test_that("run_cicero gives output", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
  data("human.hg19.genome")

  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  set.seed(2000)
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
  cons <- run_cicero(cicero_cds, sample_genome, window = 500000, silent=TRUE,
                     sample_num = 2)
  expect_equal(cons$coaccess[1], 0.9083997, tolerance = 1e-6)
  expect_equal(ncol(cons), 3)
  expect_equal(nrow(cons), 543482)
  cons <- run_cicero(cicero_cds, window = 500000, silent=TRUE, sample_num = 2,
                     genomic_coords = "../human.hg19.genome_sub.txt")
})

#### generate_ccans ####

test_that("generate_ccans gives output", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
  data("human.hg19.genome")

  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  set.seed(2000)
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
  cons <- run_cicero(cicero_cds, sample_genome, window = 500000, silent=TRUE,
                     sample_num = 2)
  expect_output(CCAN_assigns <- generate_ccans(cons),
                "Coaccessibility cutoff used: 0.7")
  expect_equal(CCAN_assigns$CCAN[3], 4, tolerance = 1e-7)
  expect_equal(ncol(CCAN_assigns), 2)
  expect_equal(nrow(CCAN_assigns), 1524)
  expect_equal(length(unique(CCAN_assigns$CCAN)), 170)

  expect_output(CCAN_assigns <- generate_ccans(cons,
                                               coaccess_cutoff_override = 0.25),
                "Coaccessibility cutoff used: 0.25")
  expect_output(CCAN_assigns <- generate_ccans(cons, tolerance_digits = 1),
                "Coaccessibility cutoff used: 0.7")
  expect_output(CCAN_assigns <- generate_ccans(cons, tolerance_digits = 1,
                                               coaccess_cutoff_override = .25),
                "Coaccessibility cutoff used: 0.25")
})

#### compare_connections ####


test_that("compare_connections works", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
  data("human.hg19.genome")

  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  set.seed(2000)
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
  cons <- run_cicero(cicero_cds, sample_genome, window = 500000, silent=TRUE,
                     sample_num = 2)
  chia_conns <-  data.frame(Peak1 = c("chr18_10000_10200", "chr18_10000_10200",
                                      "chr18_49500_49600"),
                            Peak2 = c("chr18_10600_10700", "chr18_111700_111800",
                                      "chr18_10600_10700"))
  cons$in_dataset <- compare_connections(cons, chia_conns)
  cons$in_dataset2 <- compare_connections(cons, chia_conns, maxgap=1000)

  expect_is(cons, "data.frame")
  expect_equal(sum(cons$in_dataset), 4)
  expect_equal(sum(cons$in_dataset2), 22)
  expect_equal(cons$in_dataset[1], TRUE)
})

#### find_overlapping_ccans ####

test_that("find_overlapping_ccans works", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
  data("human.hg19.genome")

  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  set.seed(2000)
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
  cons <- run_cicero(cicero_cds, sample_genome, window = 500000, silent=TRUE,
                     sample_num = 2)
  CCAN_assigns <- generate_ccans(cons, coaccess_cutoff_override = 0.25)
  over <- find_overlapping_ccans(CCAN_assigns)
  expect_is(over, "data.frame")
  expect_equal(ncol(over), 2)
  expect_equal(nrow(over), 40)
  over <- find_overlapping_ccans(CCAN_assigns, min_overlap = 3000000)
  expect_equal(nrow(over), 2)
})


#### acivity scores ####


test_that("build_gene_activity_matrix works", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
  data("human.hg19.genome")

  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  input_cds <- make_atac_cds(cicero_data, binarize=TRUE)
  input_cds <- detectGenes(input_cds, min_expr = .1)
  gene_annotation_sub <- gene_annotation_sample[,c(1:3, 8)]
  names(gene_annotation_sub)[4] <- "gene"

  input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
  unnorm_ga <- build_gene_activity_matrix(input_cds, cons)
  expect_equal(nrow(unnorm_ga), 626)
  expect_equal(ncol(unnorm_ga), 200)
  expect_equal(unnorm_ga[1,1], 1.183498, tolerance = 1e-6)

  exprs(input_cds) <- as.matrix(exprs(input_cds))
  unnorm_ga <- build_gene_activity_matrix(input_cds, cons)
  expect_equal(nrow(unnorm_ga), 626)
  expect_equal(ncol(unnorm_ga), 200)
  expect_equal(unnorm_ga[1,1], 1.183498, tolerance = 1e-6)
})


test_that("normalize_gene_activities works", {
  skip_on_bioc()
  data(cicero_data)
  load("../tsne_coord.Rdata")
  data("human.hg19.genome")

  sample_genome <- subset(human.hg19.genome, V1 == "chr18")
  input_cds <- make_atac_cds(cicero_data, binarize=TRUE)
  input_cds <- detectGenes(input_cds, min_expr = .1)
  gene_annotation_sub <- gene_annotation_sample[,c(1:3, 8)]
  names(gene_annotation_sub)[4] <- "gene"

  input_cds <- annotate_cds_by_site(input_cds, gene_annotation_sub)
  unnorm_ga <- build_gene_activity_matrix(input_cds, cons)
  num_genes <- pData(input_cds)$num_genes_expressed
  names(num_genes) <- row.names(pData(input_cds))

  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  expect_equal(nrow(cicero_gene_activities), 626)
  expect_equal(ncol(cicero_gene_activities), 200)
  expect_equal(cicero_gene_activities[1,1], 0.008026735, tolerance = 1e-6)

  cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga,
                                                           unnorm_ga),
                                                      num_genes)
  expect_is(cicero_gene_activities, "list")
  expect_equal(length(cicero_gene_activities), 2)
  cicero_gene_activities1 <- cicero_gene_activities[[1]]
  cicero_gene_activities2 <- cicero_gene_activities[[2]]
  expect_equal(nrow(cicero_gene_activities1), 626)
  expect_equal(ncol(cicero_gene_activities1), 200)
  expect_equal(cicero_gene_activities1[1,1], 0.008026735, tolerance = 1e-6)

  expect_equal(nrow(cicero_gene_activities2), 626)
  expect_equal(ncol(cicero_gene_activities2), 200)
  expect_equal(cicero_gene_activities2[1,1], 0.008026735, tolerance = 1e-6)


  unnorm_ga <- as.matrix(unnorm_ga)
  cicero_gene_activities <- normalize_gene_activities(unnorm_ga, num_genes)
  expect_equal(nrow(cicero_gene_activities), 626)
  expect_equal(ncol(cicero_gene_activities), 200)
  expect_equal(cicero_gene_activities[1,1], 0.008026735, tolerance = 1e-6)

  cicero_gene_activities <- normalize_gene_activities(list(unnorm_ga,
                                                           unnorm_ga),
                                                      num_genes)
  expect_is(cicero_gene_activities, "list")
  expect_equal(length(cicero_gene_activities), 2)
  cicero_gene_activities1 <- cicero_gene_activities[[1]]
  cicero_gene_activities2 <- cicero_gene_activities[[2]]
  expect_equal(nrow(cicero_gene_activities1), 626)
  expect_equal(ncol(cicero_gene_activities1), 200)
  expect_equal(cicero_gene_activities1[1,1], 0.008026735, tolerance = 1e-6)

  expect_equal(nrow(cicero_gene_activities2), 626)
  expect_equal(ncol(cicero_gene_activities2), 200)
  expect_equal(cicero_gene_activities2[1,1], 0.008026735, tolerance = 1e-6)

})

