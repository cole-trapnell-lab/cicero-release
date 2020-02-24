context("plotting")

df <- data.frame(Peak1 = c("chr18_10034652_10034983", "chr18_10034652_10034983",
                           "chr18_10034652_10034983", "chr18_10034652_10034983",
                           "chr18_10087586_10087901", "chr18_10120685_10127115",
                           "chr18_10097718_10097934", "chr18_10087586_10087901",
                           "chr18_10154818_10155215", "chr18_10238762_10238983",
                           "chr18_10198959_10199183", "chr18_10250985_10251585"),
                 Peak2 = c("chr18_10097718_10097934", "chr18_10087586_10087901",
                           "chr18_10154818_10155215", "chr18_10238762_10238983",
                           "chr18_10198959_10199183", "chr18_10250985_10251585",
                           "chr18_10034652_10034983", "chr18_10034652_10034983",
                           "chr18_10034652_10034983", "chr18_10034652_10034983",
                           "chr18_10087586_10087901", "chr18_10120685_10127115"),
                 coaccess = c(0.0051121787, 0.0016698617, 0.0006570246,
                              0.0013466927, 0.0737935011, 0.3264019452,
                              0.0051121787, 0.0016698617, 0.0006570246,
                              0.0013466927, 0.0737935011, 0.3264019452),
                 peak1_type = c("a","a","a","a","a","b","a","a","b","b","b","b"),
                 peak1_log = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE,
                               FALSE, FALSE, FALSE, FALSE),
                 peak_colors = c("#66C2A5", "#66C2A5", "#66C2A5", "#66C2A5",
                                 "#66C2A5", "#FC8D62", "#66C2A5", "#66C2A5",
                                 "#FC8D62", "#FC8D62", "#FC8D62", "#FC8D62" ),
                 bad_peak_colors = c("#FC8D62", "#66C2A5", "#66C2A5", "#66C2A5",
                                     "#FC8D62", "#FC8D62", "#FC8D62", "#FC8D62",
                                     "#66C2A5", "#66C2A5", "#66C2A5", "#FC8D62"),
                 connection_type = c("h", "h", "h", "h", "h", "p", "h", "h",
                                     "h", "h", "h", "p"),
                 colors = c("#66C2A5", "#66C2A5", "#66C2A5", "#66C2A5",
                            "#66C2A5", "#FC8D62", "#66C2A5", "#66C2A5",
                            "#66C2A5", "#66C2A5", "#66C2A5", "#FC8D62"))

library(data.table)
dt <- as.data.table(df)
test_that("plot_connections with data.table", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections plot dt",
                              plot_connections(dt, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585))
  })
test_that("plot_connections with bad chromosomes", {
  df$Peak1 <- gsub("chr", "A0", df$Peak1)
  df$Peak2 <- gsub("chr", "A0", df$Peak2)
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections plot bad chr",
                              plot_connections(df, chr = "A018",
                                               minbp = 10034652,
                                               maxbp = 10251585))
})


test_that("plot_connections with coaccess_cutoff", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections plot",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585))
  vdiffr::expect_doppelganger("basic connections plot cutoff",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               coaccess_cutoff = .1))
  testthat::expect_warning(vdiffr::expect_doppelganger("basic connections high cutoff",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               coaccess_cutoff = 5)))
  testthat::expect_error(plot_connections(df, chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          coaccess_cutoff = -1))
  })

test_that("plot_connections with peak_color", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("peak_color",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               peak_color = "purple"))
  vdiffr::expect_doppelganger("peak_color color column",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               peak_color = "peak_colors"))
  vdiffr::expect_doppelganger("peak_color type column",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               peak_color = "peak1_type"))
  vdiffr::expect_doppelganger("peak_color logical column",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               peak_color = "peak1_log"))
  testthat::expect_error(plot_connections(df, chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          peak_color = "bad_peak_colors"))
})

test_that("plot_connections with connection_color and connection_width and alpha_by_coaccess", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("connection_color",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_color = "purple"))
  vdiffr::expect_doppelganger("connection_color color column",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_color = "colors"))
  vdiffr::expect_doppelganger("connection_color type column",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_color = "connection_type"))
  vdiffr::expect_doppelganger("connection_color type column coaccess",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_color = "connection_type",
                                               alpha_by_coaccess = TRUE))
  vdiffr::expect_doppelganger("connection_color type column coaccess no legend",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_color = "connection_type",
                                               connection_color_legend = FALSE,
                                               alpha_by_coaccess = TRUE))
  vdiffr::expect_doppelganger("connection_color connection_width",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_color = "connection_type",
                                               connection_width = 7))
  testthat::expect_error(plot_connections(df, chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          connection_color = "connection_type",
                                          connection_width = 0))
})

test_that("plot_connections with connection_ymax", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("connection_ymax",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_ymax = 7))
  testthat::expect_warning(vdiffr::expect_doppelganger("connection_ymax plus cutoff",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               connection_ymax = 2,
                                               coaccess_cutoff = 1)))
  testthat::expect_error(plot_connections(df, chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          connection_ymax = 0))
})

test_that("plot_connections makes plots with gene model", {
  skip_on_bioc()
  data(gene_annotation_sample)
  vdiffr::expect_doppelganger("connections plot with gene model",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample))
  vdiffr::expect_doppelganger("connections plot with gene model color",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample,
                                               gene_model_color = "purple"))
  vdiffr::expect_doppelganger("connections plot with gene model no genes",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10114652,
                                               gene_model = gene_annotation_sample,
                                               gene_model_color = "purple",
                                               collapseTranscripts = "longest"))
  vdiffr::expect_doppelganger("connections plot with collapseTranscripts true",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample,
                                               collapseTranscripts = TRUE))
  vdiffr::expect_doppelganger("connections plot with collapseTranscripts longest",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample,
                                               collapseTranscripts = "longest"))
  vdiffr::expect_doppelganger("connections plot with collapseTranscripts gene",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample,
                                               collapseTranscripts = "gene"))
  vdiffr::expect_doppelganger("connections plot with collapseTranscripts shortest",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample,
                                               collapseTranscripts = "shortest"))
  vdiffr::expect_doppelganger("connections plot with collapseTranscripts meta",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample,
                                               collapseTranscripts = "meta"))
})

test_that("plot_connections with comparison_coaccess_cutoff", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections comparison plot",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585))
  vdiffr::expect_doppelganger("basic connections plot comparison cutoff",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_coaccess_cutoff = .1))
  testthat::expect_warning(vdiffr::expect_doppelganger("basic connections high comparison cutoff",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_coaccess_cutoff = 5)))
  testthat::expect_error(plot_connections(df, chr = "chr18",
                                          comparison_track = df,
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          comparison_coaccess_cutoff = -1))
})

test_that("plot_connections with comparison_peak_color", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("comparison_peak_color",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_peak_color = "purple"))
  vdiffr::expect_doppelganger("comparison_peak_color color column",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_peak_color = "peak_colors"))
  vdiffr::expect_doppelganger("comparison_peak_color type column",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_peak_color = "peak1_type"))
  vdiffr::expect_doppelganger("comparison_peak_color logical column",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_peak_color = "peak1_log"))
  testthat::expect_error(plot_connections(df, comparison_track = df,
                                          chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          comparison_peak_color = "bad_peak_colors"))
})

test_that("plot_connections with comparison_connection_color and comparison_connection_width and alpha_by_coaccess", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("comparison_connection_color",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_connection_color = "purple"))
  vdiffr::expect_doppelganger("comparison_connection_color color column",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_connection_color = "colors"))
  vdiffr::expect_doppelganger("comparison_connection_color type column",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_connection_color = "connection_type"))
  vdiffr::expect_doppelganger("comparison_connection_color type column coaccess",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_connection_color = "connection_type",
                                               alpha_by_coaccess = TRUE))
  vdiffr::expect_doppelganger("comparison_connection_color type column coaccess no legend",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_connection_color = "connection_type",
                                               comparison_connection_color_legend = FALSE,
                                               alpha_by_coaccess = TRUE))
  vdiffr::expect_doppelganger("comparison_connection_color comparison_connection_width",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_connection_color = "connection_type",
                                               comparison_connection_width = 7))
  testthat::expect_error(plot_connections(df, comparison_track = df,
                                          chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          comparison_connection_color = "connection_type",
                                          comparison_connection_width = 0))
})

test_that("plot_connections with comparison_ymax", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("comparison_ymax",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_ymax = 7))
  testthat::expect_warning(vdiffr::expect_doppelganger("comparison_ymax plus cutoff",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_ymax = 2,
                                               comparison_coaccess_cutoff = 1)))
  testthat::expect_error(plot_connections(df, comparison_track = df,
                                          chr = "chr18",
                                          minbp = 10034652,
                                          maxbp = 10251585,
                                          comparison_ymax = 0))
})

test_that("plot_connections makes plots with gene model with comparison", {
  skip_on_bioc()
  data(gene_annotation_sample)
  vdiffr::expect_doppelganger("connections plot with gene model with comparison",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               gene_model = gene_annotation_sample))
})

test_that("plot_connections makes plots with comparison dataset", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("connections plot with comparison",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_track = df,
                                               include_axis_track = FALSE))
  vdiffr::expect_doppelganger("connections plot with comparison color",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_track = df,
                                               comparison_connection_color = "connection_type"))
  vdiffr::expect_doppelganger("connections plot with comparison peak color",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_track = df,
                                               comparison_peak_color = "peak1_type"))
  vdiffr::expect_doppelganger("connections plot with comparison peak color hex",
                              plot_connections(df, chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               comparison_track = df,
                                               comparison_connection_color = "colors"))
})

test_that("plot_connections makes plots with viewpoint", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections plot with viewpoint",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               viewpoint = "chr18_10119685_10128115"))
  vdiffr::expect_doppelganger("basic connections plot with viewpoint no comp",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               viewpoint = "chr18_10119685_10128115",
                                               comparison_viewpoint=FALSE))
  vdiffr::expect_doppelganger("basic connections plot with viewpoint change colors",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               viewpoint = "chr18:10,119,685-10,128,115",
                                               viewpoint_color = "purple",
                                               viewpoint_fill = "orange",
                                               viewpoint_alpha = .1))
})

test_that("plot_connections no axis", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections include_axis_track",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               include_axis_track = FALSE))
})

df$chr_1 <- "chr18"
df$chr_2 <- 18

test_that("plot_connections given chr info", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections chr",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               include_axis_track = FALSE))
})

df$bp1_1 <- c(10034652, 10034652, 10034652, 10034652, 10087586, 10120685,
              10097718, 10087586, 10154818, 10238762, 10198959, 10250985)

df$bp2_1 <- c(10034983, 10034983, 10034983, 10034983, 10087901, 10127115,
              10097934, 10087901, 10155215, 10238983, 10199183, 10251585)

test_that("plot_connections given chr partial info", {
  skip_on_bioc()
  vdiffr::expect_doppelganger("basic connections chr bp1",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               include_axis_track = FALSE))
})

test_that("plot_connections given all partial info", {
  skip_on_bioc()
  df$chr_1 <- "18"

  df$bp1_2 <- c(10097718, 10087586, 10154818, 10238762, 10198959, 10250985,
                10034652, 10034652, 10034652, 10034652, 10087586, 10120685)

  df$bp2_2 <- c(10097934, 10087901, 10155215, 10238983, 10199183, 10251585,
                10034983, 10034983, 10034983, 10034983, 10087901, 10127115)
  vdiffr::expect_doppelganger("basic connections all bp",
                              plot_connections(df, comparison_track = df,
                                               chr = "chr18",
                                               minbp = 10034652,
                                               maxbp = 10251585,
                                               include_axis_track = FALSE))
})

test_that("plot_accessibility_in_pseudotime works", {
  skip_on_bioc()
  input_cds <- make_atac_cds(cicero_data)

  pData(input_cds)$Pseudotime <- c(3.48, 2.19, 0.70, 0.95, 3.19, 2.31, 2.26,
                                  1.54, 4.00, 2.87, 1.11, 2.36, 2.11, 0.58,
                                  3.58, 1.02, 0.15, 4.30, 3.59, 3.00, 2.71,
                                  0.18, 3.63, 3.02, 3.31, 3.61, 2.70, 0.85,
                                  0.11, 0.13, 3.54, 2.75, 2.51, 2.27, 0.24,
                                  3.14, 2.80, 2.21, 0.26, 3.60, 3.07, 2.45,
                                  0.23, 0.21, 0, 4.84, 2.03, 3.13, 3.40, 2.03,
                                  2.38, 2.65, 0.75, 4.65, 3.59, 2.79, 4.36,
                                  3.66, 0.15, 0.23, 3.00, 4.46, 2.95, 2.81,
                                  3.10, 0.94, 3.60, 3.68, 2.46, 2.87, 1.94,
                                  0.10, 2.29, 0.60, 2.23, 0.15, 4.17, 2.74,
                                  4.44, 2.65, 2.22, 2.93, 3.23, 3.57, 1.75,
                                  0.21, 1.92, 0.24, 3.54, 2.77, 2.57, 2.42,
                                  3.72, 4.57, 3.57, 3.56, 3.66, 3.55, 2.02,
                                  3.08, 2.13, 2.21, 4.43, 0.14, 2.92, 0.17,
                                  3.02, 3.61, 2.52, 3.12, 2.2, 1.60, 4.19,
                                  3.76, 1.35, 2.40, 0.72, 3.30, 3.23, 4.40,
                                  3.23, 3.55, 4.39, 2.20, 0.17, 0.28, 2.93,
                                  3.63, 4.51, 2.22, 2.91, 4.15, 1.87, 3.30,
                                  2.96, 4.18, 0.44, 3.59, 4.52, 0.46, 3.55,
                                  3.13, 4.42, 3.59, 1.12, 2.78, 0.81, 3.34,
                                  2.37, 3.00, 1.89, 0.18, 3.68, 1.04, 2.81,
                                  3.00, 2.95, 0.26, 3.24, 4.52, 2.90, 4.35,
                                  2.79, 2.19, 3.13, 1.10, 2.13, 4.12, 2.80,
                                  0.06, 0.11, 2.54, 0.12, 3.66, 1.30, 2.93,
                                  0.26, 0.05, 4.33, 4.26, 2.84, 3.04, 1.02,
                                  2.81, 4.17, 1.60, 2.76, 2.77, 2.24, 0.19,
                                  3.35, 3.05, 0.19, 3.00, 0.14, 2.15, 3.93,
                                  1.16, 0.10, 2.01)



  vdiffr::expect_doppelganger("basic bar",
                              plot_accessibility_in_pseudotime(
                                input_cds[c("chr18_38156577_38158261",
                                            "chr18_48373358_48374180",
                                            "chr18_60457956_60459080"),]))
  vdiffr::expect_doppelganger("basic bar high breaks",
                              plot_accessibility_in_pseudotime(
                                input_cds[c("chr18_38156577_38158261",
                                            "chr18_48373358_48374180",
                                            "chr18_60457956_60459080"),],
                                breaks = 20))
  vdiffr::expect_doppelganger("basic bar one",
                              plot_accessibility_in_pseudotime(
                                input_cds[c("chr18_38156577_38158261"),],
                                breaks = 4))
  testthat::expect_error(plot_accessibility_in_pseudotime(
                           input_cds[c("chr18_38156577_38158261",
                                       "chr18_48373358_48374180",
                                       "chr18_60457956_60459080"),],
                           breaks = 1),
                         "breaks not greater than or equal to 2")
  testthat::expect_error(plot_accessibility_in_pseudotime(input_cds,
                                                          breaks = 10),
                         paste("Too many sites to plot. Be sure you are",
                               "passing only a subset of your CDS.",
                               collapse = " "))
  testthat::expect_output(plot_accessibility_in_pseudotime(input_cds[2,],
                                breaks = 4), "Error! Curve fit failed!")
})

