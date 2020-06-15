context("test-utils.R")


test_that("make_atac_cds makes a valid cds object", {
  #skip_on_bioc()
  data("cicero_data")
  #### make_atac_cds ####
  test_cds <- make_atac_cds(cicero_data)
  expect_is(test_cds, "CellDataSet")
  expect_equal(nrow(exprs(test_cds)), 6146)
  expect_equal(ncol(exprs(test_cds)), 200)
  expect_match(row.names(test_cds)[1], "chr18_10025_10225")
  expect_match(colnames(test_cds)[1], "AGCGATAGAACGAATTCGGCGCAATGACCCTATCCT")
  expect_is(exprs(test_cds), "dgCMatrix")
  test_cds <-make_atac_cds(cicero_data, binarize=TRUE)
  expect_is(test_cds, "CellDataSet")
  expect_equal(nrow(exprs(test_cds)), 6146)
  expect_equal(ncol(exprs(test_cds)), 200)
  expect_match(row.names(test_cds)[1], "chr18_10025_10225")
  expect_match(colnames(test_cds)[1], "AGCGATAGAACGAATTCGGCGCAATGACCCTATCCT")
  expect_is(exprs(test_cds), "dgCMatrix")
  expect_error(test_cds <- make_atac_cds(3),
               "Input must be file path, matrix, or data.frame")
  test_cds <-make_atac_cds("../cicero_data_sub.txt", binarize=TRUE)
  expect_is(test_cds, "CellDataSet")
  expect_equal(nrow(exprs(test_cds)), 2149)
  expect_equal(ncol(exprs(test_cds)), 7)
  expect_match(row.names(test_cds)[1], "chr18_10025_10225")
  expect_match(colnames(test_cds)[1], "AGCGATAGGCGCTATGGTGGAATTCAGTCAGGACGT")
  expect_is(exprs(test_cds), "dgCMatrix")
})

#### ranges_for_coords ####

test_that("ranges_for_coords works", {
  #skip_on_bioc()

  wn <- ranges_for_coords("chr1:2039-30239", with_names = TRUE)
  wmd <- ranges_for_coords(c("chr1:2049-203902", "chrX:489249-1389389"),
                         meta_data_df = data.frame(dat = c("1", "X")))
  wmdn <- ranges_for_coords(c("chr1:2049-203902", "chrX:489249-1389389"),
                          with_names = TRUE,
                          meta_data_df = data.frame(dat = c("1", "X"),
                                                    stringsAsFactors = FALSE))


  expect_is(ranges_for_coords("chr1_2039_30239"), "GRanges")
  expect_is(ranges_for_coords("chr1_random_2039_30239"), "GRanges")
  expect_is(ranges_for_coords("chr1:2039:30239"), "GRanges")
  expect_is(ranges_for_coords("chr1-2039-30239"), "GRanges")
  expect_is(ranges_for_coords("chr1:2,039-30,239"), "GRanges")
  expect_is(ranges_for_coords(c("chr1:2,039-30,239", "chrX:28884:101293")),
            "GRanges")
  expect_is(ranges_for_coords(c("chr1:2,039-30,239", "chrX:28884:101293"),
                              with_names = TRUE), "GRanges")
  expect_is(wn, "GRanges")
  expect_is(wmd, "GRanges")
  expect_match(wn$coord_string, "chr1:2039-30239")
  expect_match(as.character(wmd$dat[2]), "X")
  expect_match(wmdn$coord_string[1], "chr1:2049-203902")
  expect_match(as.character(wmdn$dat[2]), "X")
})

#### df_for_coords ####

test_that("df_for_coords works", {
  #skip_on_bioc()
  expect_is(df_for_coords(c("chr1:2,039-30,239", "chrX:28884:101293")),
            "data.frame")
  expect_equal(df_for_coords(c("chr1:2,039-30,239",
                               "chrX:28884:101293"))$bp2[1], 30239)
  
  expect_is(df_for_coords(c("chr1:2,039-30,239", "chrX:28884:101293", 
                            "chr1_random_2039_30239")),
            "data.frame")
  expect_equal(df_for_coords(c("chr1:2,039-30,238", "chrX:28884:101293", 
                               "chr1_random_2039_30239"))$bp2[3], 30239)
})

#### annotate_cds_by_site ####



test_that("annotate_cds_by_site works", {
  #skip_on_bioc()
  data("cicero_data")
  #### make_atac_cds ####
  test_cds <- make_atac_cds(cicero_data)

  feat <- data.frame(chr = c("chr18", "chr18", "chr18", "chr18"),
                   bp1 = c(10000, 10800, 50000, 100000),
                   bp2 = c(10700, 11000, 60000, 110000),
                   type = c("Acetylated", "Methylated",
                            "Acetylated", "Methylated"),
                   stringsAsFactors = FALSE)

  test_cds2 <- annotate_cds_by_site(test_cds, feat, verbose = TRUE)
  test_cds3 <- annotate_cds_by_site(test_cds, feat, all=TRUE, verbose = TRUE)

  expect_is(test_cds2, "CellDataSet")
  expect_is(test_cds3, "CellDataSet")
  expect_equal(nrow(fData(test_cds2)), nrow(fData(test_cds)))
  expect_equal(nrow(fData(test_cds3)), nrow(fData(test_cds)))
  expect_equal(ncol(fData(test_cds2)), ncol(fData(test_cds)) + 2)
  expect_equal(ncol(fData(test_cds3)), ncol(fData(test_cds)) + 2)

  expect_equal(fData(test_cds2)$overlap[2], 201)
  expect_equal(fData(test_cds3)$overlap[2], "98,201")
  expect_equal(fData(test_cds2)$type[2], "Methylated")
  expect_equal(fData(test_cds3)$type[2], "Acetylated,Methylated")

  expect_true(is.na(fData(test_cds2)$overlap[3]))
  expect_true(is.na(fData(test_cds3)$overlap[3]))
  expect_true(is.na(fData(test_cds2)$type[3]))
  expect_true(is.na(fData(test_cds3)$type[3]))

  test_cds2 <- annotate_cds_by_site(test_cds, feat)
  test_cds3 <- annotate_cds_by_site(test_cds, feat, all=TRUE)

  expect_is(test_cds2, "CellDataSet")
  expect_is(test_cds3, "CellDataSet")
  expect_equal(nrow(fData(test_cds2)), nrow(fData(test_cds)))
  expect_equal(nrow(fData(test_cds3)), nrow(fData(test_cds)))
  expect_equal(ncol(fData(test_cds2)), ncol(fData(test_cds)) + 2)
  expect_equal(ncol(fData(test_cds3)), ncol(fData(test_cds)) + 2)

  expect_equal(fData(test_cds2)$overlap[2], 201)
  expect_equal(fData(test_cds3)$overlap[2], "98,201")
  expect_equal(fData(test_cds2)$type[2], "Methylated")
  expect_equal(fData(test_cds3)$type[2], "Acetylated,Methylated")

  expect_true(is.na(fData(test_cds2)$overlap[3]))
  expect_true(is.na(fData(test_cds3)$overlap[3]))
  expect_true(is.na(fData(test_cds2)$type[3]))
  expect_true(is.na(fData(test_cds3)$type[3]))

  test_cds2 <- annotate_cds_by_site(test_cds, "../feat.txt", verbose =TRUE)
  test_cds3 <- annotate_cds_by_site(test_cds, "../feat.txt", all=TRUE)

  expect_is(test_cds2, "CellDataSet")
  expect_is(test_cds3, "CellDataSet")
  expect_equal(nrow(fData(test_cds2)), nrow(fData(test_cds)))
  expect_equal(nrow(fData(test_cds3)), nrow(fData(test_cds)))
  expect_equal(ncol(fData(test_cds2)), ncol(fData(test_cds)) + 2)
  expect_equal(ncol(fData(test_cds3)), ncol(fData(test_cds)) + 2)

  expect_equal(fData(test_cds2)$overlap[2], 201)
  expect_equal(fData(test_cds3)$overlap[2], "98,201")
  expect_equal(fData(test_cds2)$V4[2], "Methylated")
  expect_equal(fData(test_cds3)$V4[2], "Acetylated,Methylated")

  expect_true(is.na(fData(test_cds2)$overlap[3]))
  expect_true(is.na(fData(test_cds3)$overlap[3]))
  expect_true(is.na(fData(test_cds2)$V4[3]))
  expect_true(is.na(fData(test_cds3)$V4[3]))

  test_cds2 <- annotate_cds_by_site(test_cds, "../feat_head.txt",
                                    header = TRUE)
  test_cds3 <- annotate_cds_by_site(test_cds, "../feat_head.txt",
                                    header = TRUE, all=TRUE)

  expect_is(test_cds2, "CellDataSet")
  expect_is(test_cds3, "CellDataSet")
  expect_equal(nrow(fData(test_cds2)), nrow(fData(test_cds)))
  expect_equal(nrow(fData(test_cds3)), nrow(fData(test_cds)))
  expect_equal(ncol(fData(test_cds2)), ncol(fData(test_cds)) + 2)
  expect_equal(ncol(fData(test_cds3)), ncol(fData(test_cds)) + 2)

  expect_equal(fData(test_cds2)$overlap[2], 201)
  expect_equal(fData(test_cds3)$overlap[2], "98,201")
  expect_equal(fData(test_cds2)$type[2], "Methylated")
  expect_equal(fData(test_cds3)$type[2], "Acetylated,Methylated")

  expect_true(is.na(fData(test_cds2)$overlap[3]))
  expect_true(is.na(fData(test_cds3)$overlap[3]))
  expect_true(is.na(fData(test_cds2)$type[3]))
  expect_true(is.na(fData(test_cds3)$type[3]))

  # check tie
  feat2 <- data.frame(chr = c("chr18", "chr18", "chr18", "chr18", "chr18_GL456216_random"),
                      bp1 = c(10125, 10125, 50000, 100000, 32820116),
                      bp2 = c(10703, 10703, 60000, 110000, 32820118),
                     type = c("Acetylated", "Methylated",
                              "Acetylated", "Methylated", "Other"),
                     stringsAsFactors = FALSE)
  test_cds2 <- annotate_cds_by_site(test_cds, feat2, all=FALSE)
  expect_equal(fData(test_cds2)$type[2], "Acetylated")
  test_cds2 <- annotate_cds_by_site(test_cds, feat2, all=FALSE, maxgap = 901)
  expect_equal(fData(test_cds2)$type[3], "Acetylated")

  # check maxgap = "nearest"
  test_cds2 <- annotate_cds_by_site(test_cds, feat2, all=FALSE, maxgap = "nearest")
  expect_equal(sum(is.na(fData(test_cds2)$type)), 0)

})


#### make_sparse_matrix ####

test_that("make_sparse_matrix works", {
  #skip_on_bioc()
  df <- data.frame(icol = c("chr18_30209631_30210783",
                            "chr18_45820294_45821666",
                            "chr18_32820116_32820994"),
                   jcol = c("chr18_41888433_41890138",
                            "chr18_33038287_33039444",
                            "chr18_random_25533921_25534483"),
                   xcol = c(1,2,3))
  sm <- make_sparse_matrix(df, "icol", "jcol", "xcol")
  expect_equal(sm["chr18_30209631_30210783", "chr18_41888433_41890138"], 1)
  expect_equal(sm["chr18_45820294_45821666", "chr18_33038287_33039444"], 2)
  expect_equal(sm["chr18_random_25533921_25534483", "chr18_32820116_32820994"], 3)
  expect_equal(sm["chr18_random_25533921_25534483", "chr18_30209631_30210783"], 0)
  expect_error(make_sparse_matrix(df, "icol", "xcol", "jcol"),
               "x.name column must be numeric")
  expect_error(make_sparse_matrix(df, "icol", "hannah", "jcol"),
               "i.name, j.name, and x.name must be columns in data")
})

#### compare_connections ####
# IN test-runCicero.R

#### find_overlapping_coordinates ####



test_that("find_overlapping_coordinates works", {
  #skip_on_bioc()
  test_coords <- c("chr18_10025_10225", "chr18_10603_11103", "chr18_11604_13986",
                 "chr18_157883_158536", "chr18_217477_218555",
                 "chr18_245734_246234", "chr18_random_245734_246234")
  expect_equal(length(find_overlapping_coordinates(test_coords,
                                                   "chr18:10,100-1246234")), 6)
  expect_equal(length(find_overlapping_coordinates(test_coords,
                                                   "chr18_10227_10601")), 0)
  expect_equal(length(find_overlapping_coordinates(test_coords,
                                                   "chr18_10227_10601",
                                                   maxgap = 1)), 2)
  expect_equal(length(find_overlapping_coordinates(test_coords,
                                                   "chr18_random_10227_245736",
                                                   maxgap = 1)), 1)
  expect_equal(length(find_overlapping_coordinates(test_coords,
                                                   c("chr18_10227_10602",
                                                     "chr18:11604-246234"))), 5)
  expect_equal(length(find_overlapping_coordinates(test_coords,
                                                   c("chr18_10226_10602",
                                                     "chr18:11604-246234"),
                                                   maxgap = 1)), 6)
  expect_true(all(is.na(find_overlapping_coordinates(test_coords,
                                                   c("chr19_10226_10602",
                                                     "chr19:11604-246234"),
                                                   maxgap = 1))))
  expect_true(all(is.na(find_overlapping_coordinates(test_coords,
                                                c("chr18_1022600_1060200",
                                                  "chr18:1160400-24623400"),
                                                maxgap = 1))))
})






