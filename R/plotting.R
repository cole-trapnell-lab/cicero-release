monocle_theme_opts <- function()
{
  theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank()) +
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

#' Plot connections
#'
#' Plotting function for Cicero connections. Uses \code{\link[Gviz]{plotTracks}}
#' as its basis
#'
#' @param connection_df Data frame of connections, which must include the
#'   columns 'Peak1', 'Peak2', and 'coaccess'. Generally, the output of
#'   run_cicero or assemble_connections.
#' @param chr The chromosome of the region you would like to plot in the form
#'   'chr10'.
#' @param minbp The base pair coordinate of the start of the region to be
#'   plotted.
#' @param maxbp The base pair coordinate of the end of the region to be plotted.
#' @param coaccess_cutoff The minimum cicero co-accessibility score you would
#'   like to be plotted. Default is 0.
#' @param peak_color Color for peak annotations - a single color, the name of a
#'   column containing color values that correspond to Peak1, or the name of
#'   column containing a character or factor to base peak colors on.
#' @param connection_color Color for connection lines. A single color, the name
#'   of a column containing color values, or the name of a column containing a
#'   character or factor to base connection colors on.
#' @param alpha_by_coaccess Logical, should the transparency of connection
#'   lines be scaled based on co-accessibility score?
#' @param connection_width Width of connection lines.
#' @param connection_ymax Connection y-axis height. If \code{NULL}, chosen
#'   automatically.
#' @param gene_model Either \code{NULL} or a data.frame. The data.frame should
#'   be in a form compatible with the Gviz function
#'   \code{\link[Gviz]{GeneRegionTrack-class}} (cannot have NA as column names).
#' @param gene_model_color Color for gene annotations.
#' @param gene_model_shape Shape for gene models, passed to
#'   \code{\link[Gviz]{GeneRegionTrack-class}}. Options described at
#'   \code{\link[Gviz]{GeneRegionTrack-class}}.
#' @param comparison_track Either \code{NULL} or a data frame. If a data frame,
#'   a second track of connections will be plotted based on this data. This
#'   data frame has the same requirements as connection_df (Peak1, Peak2 and
#'   coaccess columns).
#' @param comparison_coaccess_cutoff The minimum cicero co-accessibility score
#'   you would like to be plotted for the comparison dataset. Default = 0.
#' @param comparison_peak_color Color for comparison peak annotations - a
#'   single color, the name of a column containing color values that correspond
#'   to Peak1, or the name of a column containing a character or factor to base
#'   peak colors on.
#' @param comparison_connection_color Color for comparison connection lines. A
#'   single color, the name of a column containing color values, or the name of
#'   a column containing a character or factor to base connection colors on.
#' @param comparison_connection_width Width of comparison connection lines.
#' @param comparison_ymax Connection y-axis height for comparison track. If
#'   \code{NULL}, chosen automatically.
#' @param collapseTranscripts Logical or character scalar. Can be one in
#'   \code{gene}, \code{longest}, \code{shortest} or \code{meta}. Variable is
#'   passed to the \code{\link[Gviz]{GeneRegionTrack-class}} function of Gviz.
#'   Determines whether and how to collapse related transcripts. See Gviz
#'   documentation for details.
#' @param include_axis_track Logical, should a genomic axis be plotted?
#' @param return_as_list Logical, if TRUE, the function will not plot, but will
#'   return the plot components as a list. Allows user to add/customize Gviz
#'   components and plot them separately using \code{\link[Gviz]{plotTracks}}.
#' @param viewpoint \code{NULL} or Coordinates in form "chr1_10000_10020". Use
#'   viewpoint if you would like to plot cicero connections "4C-seq style".
#'   Only connections originating in the viewpoint will be shown. Ideal for
#'   comparisons with 4C-seq data. If comparison_viewpoint is \code{TRUE}, any
#'   comparison track will be subsetted as well.
#' @param comparison_viewpoint Logical, should viewpoint apply to comparison
#'   track as well?
#' @param viewpoint_color Color for the highlight border.
#' @param viewpoint_fill Color for the highlight fill.
#' @param viewpoint_alpha Alpha value for the highlight fill.
#' @param connection_color_legend Logical, should connection color legend be
#'   shown?
#' @param comparison_connection_color_legend Logical, should comparison
#'   connection color legend be shown?
#'
#' @return A gene region plot, or list of components if return_as_list is
#'   \code{TRUE}.
#' @export
#' @import Gviz
#'
#' @examples
#'   cicero_cons <- data.frame(
#'              Peak1 = c("chr18_10034652_10034983", "chr18_10034652_10034983",
#'                        "chr18_10034652_10034983", "chr18_10034652_10034983",
#'                        "chr18_10087586_10087901", "chr18_10120685_10127115",
#'                        "chr18_10097718_10097934", "chr18_10087586_10087901",
#'                        "chr18_10154818_10155215", "chr18_10238762_10238983",
#'                        "chr18_10198959_10199183", "chr18_10250985_10251585"),
#'              Peak2 = c("chr18_10097718_10097934", "chr18_10087586_10087901",
#'                        "chr18_10154818_10155215", "chr18_10238762_10238983",
#'                        "chr18_10198959_10199183", "chr18_10250985_10251585",
#'                        "chr18_10034652_10034983", "chr18_10034652_10034983",
#'                        "chr18_10034652_10034983", "chr18_10034652_10034983",
#'                        "chr18_10087586_10087901", "chr18_10120685_10127115"),
#'              coaccess = c(0.0051121787, 0.0016698617, 0.0006570246,
#'                           0.0013466927, 0.0737935011, 0.3264019452,
#'                           0.0051121787, 0.0016698617, 0.0006570246,
#'                           0.0013466927, 0.0737935011, 0.3264019452))
#'   plot_connections(cicero_cons, chr = "chr18",
#'                    minbp = 10034652,
#'                    maxbp = 10251585,
#'                    peak_color = "purple")
#'
plot_connections <- function(connection_df,
                             chr,
                             minbp,
                             maxbp,
                             coaccess_cutoff = 0,
                             peak_color = "#B4656F",
                             connection_color = "#7F7CAF",
                             connection_color_legend = TRUE,
                             alpha_by_coaccess = FALSE,
                             connection_width = 2,
                             connection_ymax = NULL,
                             gene_model = NULL,
                             gene_model_color = "#81D2C7",
                             gene_model_shape = c("smallArrow", "box"),
                             comparison_track = NULL,
                             comparison_coaccess_cutoff = 0,
                             comparison_peak_color = "#B4656F",
                             comparison_connection_color = "#7F7CAF",
                             comparison_connection_color_legend = TRUE,
                             comparison_connection_width = 2,
                             comparison_ymax = NULL,
                             collapseTranscripts = FALSE,
                             include_axis_track = TRUE,
                             return_as_list = FALSE,
                             viewpoint = NULL,
                             comparison_viewpoint = TRUE,
                             viewpoint_color = "#F0544F",
                             viewpoint_fill = "#EFD8D7",
                             viewpoint_alpha = 0.5) {
  # Check inputs:

  assertthat::assert_that(is.data.frame(connection_df))
  assertthat::assert_that(assertthat::has_name(connection_df, "Peak1"),
                          assertthat::has_name(connection_df, "Peak2"),
                          assertthat::has_name(connection_df, "coaccess"))

  #assertthat::assert_that(is.character(chr), is_chr(chr))
  assertthat::assert_that(assertthat::is.number(minbp),
                          assertthat::is.number(maxbp))

  assertthat::assert_that(assertthat::is.number(coaccess_cutoff))
  assertthat::assert_that(coaccess_cutoff >=0)

  assertthat::assert_that(is_color(peak_color, connection_df))

  assertthat::assert_that(is_color(connection_color, connection_df))

  assertthat::assert_that(is.logical(alpha_by_coaccess))

  assertthat::assert_that(is.numeric(connection_width), connection_width > 0)

  assertthat::assert_that(is.null(connection_ymax) |
                            assertthat::is.number(connection_ymax))
  if (!is.null(connection_ymax)) assertthat::assert_that(connection_ymax > 0)

  if (!is.null(gene_model)) {
    assertthat::assert_that(is.data.frame(gene_model),
                            assertthat::has_name(gene_model, "chromosome"),
                            assertthat::has_name(gene_model, "start"),
                            assertthat::has_name(gene_model, "end"),
                            assertthat::has_name(gene_model, "strand"),
                            assertthat::has_name(gene_model, "transcript"))
    assertthat::assert_that(class(gene_model$start) %in% c("integer", "numeric"))
    assertthat::assert_that(class(gene_model$end) %in% c("integer", "numeric"))
    assertthat::assert_that(is_color(gene_model_color))
  }

  if (!is.null(comparison_track)) {
    assertthat::assert_that(is.data.frame(comparison_track))
    assertthat::assert_that(assertthat::has_name(comparison_track, "Peak1"),
                            assertthat::has_name(comparison_track, "Peak2"),
                            assertthat::has_name(comparison_track, "coaccess"))
    assertthat::assert_that(is_color(comparison_connection_color,
                                     comparison_track))
    assertthat::assert_that(assertthat::is.number(comparison_coaccess_cutoff))
    assertthat::assert_that(comparison_coaccess_cutoff >=0)
    assertthat::assert_that(is_color(comparison_peak_color, comparison_track))
    assertthat::assert_that(is.null(comparison_ymax) |
                              assertthat::is.number(comparison_ymax))
    if (!is.null(comparison_ymax)) assertthat::assert_that(comparison_ymax > 0)
    assertthat::assert_that(is.numeric(comparison_connection_width),
                            comparison_connection_width > 0)
  }

  assertthat::assert_that(is.logical(include_axis_track),
                          is.logical(return_as_list))

  assertthat::assert_that(is.logical(collapseTranscripts) |
                            is.character(collapseTranscripts))
  if(is.character(collapseTranscripts)) {
    assertthat::assert_that(collapseTranscripts %in% c("gene", "longest",
                                                       "shortest", "meta"))
  }

  if (!is.null(viewpoint)) {
    assertthat::assert_that(is.character(viewpoint))
    assertthat::assert_that(is_color(viewpoint_color, connection_df))
    assertthat::assert_that(is_color(viewpoint_fill, connection_df))
    assertthat::assert_that(assertthat::is.number(viewpoint_alpha),
                            viewpoint_alpha > 0,
                            viewpoint_alpha <= 1)
    if (!is.null(comparison_track)) {
      assertthat::assert_that(is.logical(comparison_viewpoint))
    }
    viewpoint <- gsub(",", "", viewpoint)
    viewpoint <- gsub(":|-", "_", viewpoint)
  }

  # assign peak colors
  if (peak_color %in% names(connection_df)) {
    connection_df$peak_color <- get_colors(connection_df[,peak_color])
  } else {
    connection_df$peak_color <- peak_color
  }

  # if chr_1 already exists, subset by chromosome
  if("chr_1" %in% names(connection_df)) {
    connection_df <- connection_df[connection_df$chr_1 %in%
                                     c(chr, paste0("chr", chr),
                                       gsub("chr", "", chr)),]
  }
  sub <- generate_plotting_subset(connection_df, chr, minbp, maxbp)

  if (nrow(sub) == 0) stop("No peaks in the specified range, nothing to plot!")

  if (!peak_color %in% names(connection_df)) sub$peak_color <- peak_color
  anntrack <- make_peak_track(sub)

  sub <- sub[sub$coaccess >= coaccess_cutoff &
               sub$coaccess > 0 &
               !is.na(sub$coaccess),]

  if (!is.null(viewpoint)) {
    viewpoint_overs <- find_overlapping_coordinates(c(sub$Peak1, sub$Peak2),
                                                    viewpoint)
    sub <- sub[sub$Peak1 %in% viewpoint_overs | sub$Peak2 %in% viewpoint_overs,]
    if (nrow(sub) == 0) warning("No connections in viewpoint")
  }

  color_names <- NULL
  if (!nrow(sub) == 0) {
    if (connection_color %in% names(sub)) {
      color_levs <- levels(as.factor(sub[,connection_color]))
      color_names <- rep("temp", length(color_levs))
      names(color_names) <- color_levs
      new_connection_color <- get_colors(sub[,connection_color])
      for(n in color_levs) {
          color_names[n] <-
              new_connection_color[which(sub[,connection_color] == n)[1]]
      }
      connection_color <- new_connection_color
    }
    sub$color <- connection_color

    sub$width <- connection_width

    sub <- sub[,c("chr", "bp1", "bp2", "chr_2", "bp1_2", "bp2_2", "coaccess",
                  "width", "color")]
    names(sub) <- c("chrom1","start1","stop1","chrom2","start2","stop2",
                    "height", "width", "color")
  }
  else {
    warning("No connections above coaccess_cutoff")
  }

  if (connection_color_legend == FALSE) color_names <- NULL

  if (is.null(connection_ymax))
    connection_ymax <- signif_up(max(sub$height, na.rm=TRUE))
  if (is.na(connection_ymax) | connection_ymax <= 0) {
    connection_ymax <- 1
    warning("connection_ymax calc failed")
  }

  ctrack <- CustomTrack(plottingFunction=function(GdObject, prepare) {
    Gviz::displayPars(GdObject) <- list(ylim = c(0,connection_ymax))
    if(!prepare) {
      plotBedpe(sub, chrom = chr, chromstart=minbp, chromend=maxbp,
                connection_ymax, coaccess_cutoff,
                connection_width, alpha_by_coaccess,
                color_names)
     }
    return(invisible(GdObject))}, name="", fontsize.group=6,fontsize=6)

  component_list <- list(ctrack, anntrack)
  component_num <- 2
  size_track <- c(1,.3)

  if (include_axis_track) {
    atrack <- Gviz::GenomeAxisTrack(fontsize=6)
    size_track <- c(size_track, 0.2)
    component_num <- component_num + 1
    component_list[[component_num]] <- atrack
  }

  if(!is.null(comparison_track)) {
    if (comparison_peak_color %in% names(comparison_track)) {
      comparison_track$peak_color <-
        get_colors(comparison_track[,comparison_peak_color])
    } else {
      comparison_track$peak_color <- comparison_peak_color
    }

    sub2 <- generate_plotting_subset(comparison_track, chr, minbp, maxbp)
    if (!nrow(sub2) == 0) {
      if (!comparison_peak_color %in% names(comparison_track)) {
        sub2$peak_color <- comparison_peak_color
      }
    }
      anntrack2 <- make_peak_track(sub2)
    color_names2 <- NULL
    if (!nrow(sub2) == 0) {
      if (comparison_connection_color %in% names(sub2)) {
        color_levs <- levels(as.factor(sub2[,comparison_connection_color]))
        color_names2 <- rep("temp", length(color_levs))
        names(color_names2) <- color_levs
        new_connection_color <- get_colors(sub2[,comparison_connection_color])
        for(n in color_levs) {
          color_names2[n] <-
            new_connection_color[which(sub2[,comparison_connection_color] ==
                                         n)[1]]
        }
        comparison_connection_color <- new_connection_color
      }
      if (comparison_connection_color_legend == FALSE) color_names2 <- NULL

      sub2$color <- comparison_connection_color

      sub2 <- sub2[sub2$coaccess >= comparison_coaccess_cutoff &
                     !is.na(sub2$coaccess),]
      if (!nrow(sub2) == 0 &
          !is.null(viewpoint) &
          comparison_viewpoint == TRUE) {
        viewpoint_overs <- find_overlapping_coordinates(c(sub2$Peak1,
                                                          sub2$Peak2),
                                                        viewpoint)
        sub2 <- sub2[sub2$Peak1 %in% viewpoint_overs |
                       sub2$Peak2 %in% viewpoint_overs,]
        if (nrow(sub2) == 0) warning("No comparison connections in viewpoint")
      }
    }
    if (!nrow(sub2) == 0) {
      sub2$width <- comparison_connection_width
      sub2 <- sub2[,c("chr", "bp1", "bp2", "chr_2", "bp1_2", "bp2_2",
                      "coaccess", "width", "color")]
      names(sub2) <- c("chrom1","start1","stop1","chrom2","start2","stop2",
                       "height", "width", "color")
    }

    if (is.null(comparison_ymax))
      comparison_ymax <- signif_up(max(sub2$height, na.rm=TRUE))
    if (is.na(comparison_ymax) | comparison_ymax <= 0) {
      comparison_ymax <- 1
      warning("comparison_ymax calc failed")
    }

    ctrack2 <- CustomTrack(plottingFunction=function(GdObject, prepare) {
     Gviz::displayPars(GdObject) <- list(ylim = c(0,comparison_ymax))
      if(!prepare) {
        plotBedpe(sub2, chrom = chr, chromstart=minbp, chromend=maxbp,
                  comparison_ymax, comparison_coaccess_cutoff,
                  comparison_connection_width, alpha_by_coaccess, color_names2)
      }
      return(invisible(GdObject))}, name="", fontsize.group=6,fontsize=6)

    component_num <- component_num + 1
    component_list[[component_num]] <- ctrack2
    size_track <- c(size_track, 1)
    component_num <- component_num + 1
    component_list[[component_num]] <- anntrack2
    size_track <- c(size_track, .3)
  }

  if(!is.null(gene_model)) {
    gene_model <- gene_model[!is.na(gene_model$chromosome) & 
                             !is.na(gene_model$start) &
                             !is.na(gene_model$end) &
                             !is.na(gene_model$strand) &
                             !is.na(gene_model$transcript),]

    gene_model <-
      gene_model[gene_model$chromosome == chr &
                   ((gene_model$start > minbp & gene_model$start < maxbp) |
                    (gene_model$end > minbp & gene_model$end < maxbp) |
                    (gene_model$start < minbp & gene_model$end > maxbp)),]
    # Define gene model track
    grtrack <- make_gene_model_track(gene_model, chr,
                                     collapseTranscripts, gene_model_color,
                                     gene_model_shape)
    component_num <- component_num + 1
    component_list[[component_num]] <- grtrack
    size_track <- c(size_track, .3)
  }

  if(!is.null(viewpoint)) {
    view_chr <- split_peak_names(viewpoint)[,1]
    if (chr != view_chr) warning("Viewpoint not on correct chromosome")
    else {
      if (return_as_list) {
        message(paste0("In order to use return_as_list functionality along",  
                       "with a viewpoint, the final step of track ", 
                       "highlighting must be skipped. After your ",
                       "modifications, you will need to use the ", 
                       "Gviz::HighlightTrack function for final plotting. See ", 
                       "this link for details: ", 
                       "https://www.bioconductor.org/packages/devel/bioc/vignettes/Gviz/inst/doc/Gviz.html#61_Highlighting"))
        return(component_list)
      }
      view_start <- split_peak_names(viewpoint)[,2]
      view_end <- split_peak_names(viewpoint)[,3]
      ht1 <- Gviz::HighlightTrack(trackList = component_list,
                                  start = as.numeric(view_start),
                                  end = as.numeric(view_end),
                                  chromosome = chr,
                                  col = viewpoint_color,
                                  fill = viewpoint_fill,
                                  inBackground = FALSE,
                                  alpha = viewpoint_alpha)
      component_list <- ht1
    }
  }

  if (return_as_list) return(component_list)
  return(plotTracks(component_list, title.width = .5, showTitle = TRUE,
              from = minbp, to = maxbp, chromosome = chr, sizes = size_track,
             transcriptAnnotation = "symbol", background.title = "transparent",
             col.border.title="transparent", lwd.border.title = "transparent",
              col.axis = "black", fontsize.group = 6,#col.title="white",
             fontcolor.legend = "black"))
}

generate_plotting_subset <- function(connections, chr, minbp, maxbp) {
  connections$Peak1 <- as.character(connections$Peak1)
  connections$Peak2 <- as.character(connections$Peak2)

  pcolor_map <- data.frame(Peak1 = connections$Peak1,
                           peak_color = connections$peak_color)
  pcolor_map <- pcolor_map[!duplicated(pcolor_map),]
  connections$peak_color <- NULL

  if(sum(!c("chr_1", "chr_2", "bp1_1", "bp2_1", "bp2_1", "bp2_2") %in%
         names(connections)) != 0 ) {
    suppressWarnings(connections$chr <- NULL)
    suppressWarnings(connections$bp1 <- NULL)
    suppressWarnings(connections$bp2 <- NULL)
    suppressWarnings(connections$chr_2 <- NULL)
    suppressWarnings(connections$bp1_2 <- NULL)
    suppressWarnings(connections$bp2_2 <- NULL)

    connections <- cbind(connections, df_for_coords(connections$Peak1)[,c(1, 2, 3)])
    cons2 <- df_for_coords(connections$Peak2) #slow
    cons2$Peak <- NULL
    names(cons2) <- c("chr_2", "bp1_2", "bp2_2")
    connections <- cbind(connections, cons2) #slow
  } else {
    if(!grepl("chr", connections$chr_1[1])) {
      connections$chr_1 <- paste0("chr", connections$chr_1)
    }
    if(!grepl("chr", connections$chr_2[1])) {
      connections$chr_2 <- paste0("chr", connections$chr_2)
    }
    names(connections)[names(connections) == "chr_1"] <- "chr"
    names(connections)[names(connections) == "bp1_1"] <- "bp1"
    names(connections)[names(connections) == "bp2_1"] <- "bp2"
  }

  sub <- connections[connections$chr_2 == chr & connections$bp1 <= maxbp &
                       connections$bp2 <= maxbp & connections$bp1 >= minbp &
                       connections$bp2 >= minbp & connections$bp1_2 <= maxbp &
                       connections$bp2_2 <= maxbp & connections$bp1_2 >= minbp &
                       connections$bp2_2 >= minbp,]

  sub2 <- sub
  names(sub2)[unlist(lapply(c("Peak1", "Peak2", "coaccess", "chr",
                       "bp1", "bp2", "chr_2", "bp1_2", "bp2_2"),
    function(x) which(names(sub2) == x)))] <- c("Peak2", "Peak1", "coaccess",
                                               "chr_2", "bp1_2", "bp2_2",
                                               "chr", "bp1", "bp2")

  sub <- rbind(sub, sub2)
  sub <- sub[!duplicated(sub),]

  sub <- merge(sub, pcolor_map, all.x = TRUE)
  sub$peak_color <- as.character(sub$peak_color)
  sub$peak_color[is.na(sub$peak_color)] <- "black"

  return(sub)
}

make_peak_track <- function(df) {
  df <- df[!duplicated(df[,c("chr", "bp1", "bp2", "peak_color")]),]

  if (sum(duplicated(df[,c("chr", "bp1", "bp2")])) > 0)
    stop(paste("Multiple peak colors correspond to a single peak. Be sure that",
               "your peak_color column name assigns colors for Peak1 only",
               collapse = " "))
  gr <-  GenomicRanges::GRanges(as.character(df$chr),
                 IRanges::IRanges(as.numeric(as.character(df$bp1)),
                         as.numeric(as.character(df$bp2))))
  options(ucscChromosomeNames=FALSE)
  anntrack <- Gviz::AnnotationTrack(gr, name="Peaks", fill = df$peak_color,
                                    lwd = .0000001, fontsize.group=6,
                                    fontsize=6, cex.feature = 0.5)

  anntrack
}

make_gene_model_track <- function(txdb,
                                  chr,
                                  collapseTranscripts,
                                  gene_model_color, shape) {
  if(nrow(txdb) != 0) {
    txdb <- txdb[txdb$chromosome == chr,]
    txdb$transcript <- as.character(txdb$transcript)
    if("feature" %in% names(txdb)) {
      txdb$feature <- as.character(txdb$feature)
    }
    grtrack <- Gviz::GeneRegionTrack(txdb, chromosome = chr, geneSymbols=TRUE,
                                     name = "Gene Model", fill=gene_model_color,
                                     col= gene_model_color,  fontcolor="black",
                                     fontcolor.group="black", fontsize.group=6,
                                     fontsize=6,
                                     collapseTranscripts = collapseTranscripts,
                                     shape=shape)
  } else {
    grtrack <- GeneRegionTrack()
  }
  grtrack
}

plotBedpe <- function(bedpedata,
                      chrom,
                      chromstart,
                      chromend,
                      ymax,
                      coaccess_cutoff,
                      width,
                      alpha_by_coaccess,
                      color_names = NULL)
{ ###### All borrowed and modified from Sushi package.

  if (nrow(bedpedata) == 0) {
    warning("Nothing to plot")
    return()
  }

  bedpedata  <- bedpedata[,c("chrom1","start1","stop1","chrom2","start2",
                             "stop2", "height", "width", "color")]

  # normalize height
  maxheight <- ymax

  bedpedata$alpha <- .6
  if(alpha_by_coaccess) {
    bedpedata$alpha <- (bedpedata$height-coaccess_cutoff)/maxheight
  }
  bedpedata$height <- bedpedata$height/maxheight
  # remove any rows with 0 height
  bedpedata <- bedpedata[abs(bedpedata$height) > 0,]

  # reclass data
  if (any(class(bedpedata) == "data.table")) {
    for(i in c("start1", "stop1", "start2", "stop2")) {
      bedpedata[[i]] <- as.numeric(as.character((bedpedata[[i]])))
    }
  } else {
    for(i in c("start1", "stop1", "start2", "stop2")) {
      bedpedata[,i] <- as.numeric(as.character((bedpedata[,i])))
    }
  }

  # add position columns
  bedpedata$pos1 = apply(bedpedata[,c("start1","stop1")],1,mean)
  bedpedata$pos2 = apply(bedpedata[,c("start2","stop2")],1,mean)

  totalrange <- as.numeric(as.character(chromend)) -
    as.numeric(as.character(chromstart))
  if (nrow(bedpedata) == 0) warning("Nothing to plot")

  #legFactors <- sort(names(which(apply(legInfo, 2, any))))
  #boxSize <-  if(length(setdiff(legFactors, c("col", "cex")))==0) 0.1 else 0.3
  #pcols <- Gviz:::.getPlottingFeatures(GdObject)

  if (!is.null(color_names)) {
    boxSize <- .3
    spacing <- 0.2
    vspace <- .05
    for (i in seq_len(length(color_names))) {
      grid::grid.lines(unit(c(spacing,spacing + boxSize), "inches"),
                       c(1 - vspace*i, 1 - vspace*i),
                       gp=grid::gpar(col=color_names[i], lwd=width))
      grid::grid.text(x=unit(.1 + (boxSize + spacing), "inches"),
                      y=1 - vspace*i, just=c(0, 0.5),
                      label=names(color_names)[i])
    }
  }

  # plot the data
  grid::grid.function(function(x) list(x=x, y=(coaccess_cutoff/(ymax))),
                      gp=grid::gpar(col="black", lty="dashed", lwd=width)) #
  for (row in (seq_len(nrow(bedpedata)))) {
    x1     = bedpedata$pos1[row]
    x2     = bedpedata$pos2[row]
    height = bedpedata$height[row]
    width  = bedpedata$width[row]
    color  = bedpedata$color[row]
    alpha  = bedpedata$alpha[row]
    plotpair(x1,x2,height,totalrange,width,color, chromstart, alpha)
  }
}

# Define a function that plots a looping interaction on a graph
plotpair <- function(start, end, height, totalrange,
                     width, color, chromstart, alpha) {
  #scale values for plotting
  x1 = (min(start,end) - as.numeric(as.character(chromstart)))/totalrange
  x2 = (max(start,end) - as.numeric(as.character(chromstart)))/totalrange
  hx1 <- (x1 + x2)/2
  hy1 <- height/.725

  grid::grid.bezier(x = c(x1, hx1, hx1, x2), y = c(0, hy1, hy1, 0),
                    default.units = "npc",
                    gp=grid::gpar(col=color, lwd=width,
                                  alpha = (alpha*.9 + .1)))
}


get_colors <- function(type_list) {
  if (is_color(type_list)) {
    return(as.character(type_list))
  }
  type_list <- as.numeric(as.factor(type_list))
  n <- length(unique(type_list))
  if(n < 2) return(rep("black", times = length(type_list)))
  return(rainbow(n)[type_list])
}


#' Plot accessibility by pseudotime
#'
#' Make a barplot of chromatin accessibility across pseudotime
#'
#' @param cds_subset Subset of the CDS object you want to plot. The CDS must
#'   have a column in the pData table called "Pseudotime".
#' @param breaks Number of breaks along pseudotime. Controls the coarseness of
#'   the plot.
#'
#' @details This function plots each site in the CDS subset by cell pseudotime
#'   as a barplot. Cells are divided into bins by pseudotime (number determined
#'   by \code{breaks}) and the percent of cells in each bin that are accessible
#'   is represented by bar height. In addition, the black line represents the
#'   pseudotime-dependent average accessibility from a smoothed binomial
#'   regression.
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' plot_accessibility_in_pseudotime(input_cds_lin[c("chr18_38156577_38158261",
#'                                                  "chr18_48373358_48374180",
#'                                                  "chr18_60457956_60459080")])
#' }
#'
plot_accessibility_in_pseudotime <- function(cds_subset,
                                             breaks = 10) {
  assertthat::assert_that(is(cds_subset, "CellDataSet"))
  assertthat::assert_that(assertthat::is.count(breaks))
  assertthat::assert_that(breaks >=2)
  assertthat::assert_that(nrow(fData(cds_subset)) <= 30,
                          msg = paste("Too many sites to plot. Be sure you are",
                                      "passing only a subset of your CDS.",
                                      collapse = " "))

  min_expr = 0
  fData(cds_subset)$site_name <- row.names(fData(cds_subset))
  
  df <- as.data.frame(as.matrix(exprs(cds_subset)))
  df$f_id <- row.names(df)
  cds_exprs <- tidyr::gather(df, "Cell", "expression", -f_id)

  cds_exprs$expression <- as.numeric(cds_exprs$expression > 0)

  cds_exprs <- merge(cds_exprs, fData(cds_subset), by.x="f_id",
                     by.y="row.names")
  cds_exprs <- merge(cds_exprs, pData(cds_subset), by.x="Cell",
                     by.y="row.names")

  trend_formula <- "expression~ sm.ns(Pseudotime, df=3)"

  merged_df_with_vgam <- plyr::ddply(cds_exprs, plyr::.(f_id), function(x) {
    fit_res <- tryCatch({
      vg <- suppressWarnings(VGAM::vgam(formula = as.formula(trend_formula),
                                        family = cds_subset@expressionFamily,
                                        data = x, maxit=30, checkwz=FALSE))
      res <- predict(vg, type="response")
      res[res < min_expr] <- min_expr
      res
    }
    ,error = function(e) {
      print("Error! Curve fit failed!")
      print(e)
      res <- rep(NA, nrow(x))
      res
    })

    expectation <- as.numeric(fit_res)

    data.frame(Pseudotime=x$Pseudotime, expectation=expectation)
  })

  cds_exprs$br <- cut(cds_exprs$Pseudotime,breaks=breaks)

  df <- as.data.frame(with(cds_exprs, tapply(expression, list(br, f_id),
                                             mean)))
  df$Var1 <- row.names(df)
  mean.wt <- tidyr::gather(df, "Var2", "value", -Var1)
  
  # fix to avoid reshape v reshape2 incompatability
  names(mean.wt) <- c("Var1", "Var2", "value")
  mean.wt <- cbind(mean.wt, stringr::str_split_fixed(mean.wt$Var1, ",", 2))

  names(mean.wt) <- c("interval", "feature_label", "mean", "int_start",
                      "int_end")

  mean.wt$int_start <- as.numeric(as.character(gsub("\\(", "",
                                                    mean.wt$int_start)))
  merged_df_with_vgam$feature_label <- merged_df_with_vgam$f_id
  mean.wt$mean <- mean.wt$mean * 100
  merged_df_with_vgam$expectation <- merged_df_with_vgam$expectation * 100

  g_plot <- ggplot2::ggplot(data=mean.wt) +
    ggplot2::geom_bar(stat="identity",
                      ggplot2::aes_string(x = "int_start", y = "mean"),
                      color="#3C1642", fill= "#1DD3B0", size = .2) +
    ggplot2::geom_line(ggplot2::aes_string(x = "Pseudotime",
                                           y = "expectation"),
                       data=merged_df_with_vgam) +
    ggplot2::facet_wrap(~feature_label, nrow=NULL, ncol=1, scales="free_y") +
    ggplot2::ylab("Percent of cells accessible") +
    ggplot2::xlab("Pseudotime") +
    ggplot2::theme(text = ggplot2::element_text(size=6)) +
    ggplot2::theme(axis.line.x = ggplot2::element_line(size = .2),
                   axis.line.y = ggplot2::element_line(size = .2)) +
    monocle_theme_opts()

  return(g_plot)
}

# This function copied from
# https://stackoverflow.com/questions/37583715/
# round-up-values-to-a-specific-significant-figure-in-r

signif_up <- function(x) {
  num_string <- format(x, scientific=TRUE)

  n <- strsplit(num_string, "e")
  n1 <- vapply(n, function(x) as.numeric(x[1]), .1)
  n2 <- vapply(n, function(x) as.numeric(x[2]), .1)

  ceiling(n1) * 10^(n2)
}
